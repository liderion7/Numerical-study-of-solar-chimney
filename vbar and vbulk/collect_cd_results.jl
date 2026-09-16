# collect_cd_results.jl
# Collect Task-D outputs from all completed angle cases and build Vbulk(theta).
# Handles both standard task_d/task_d_result.csv and
# task_d_manual/task_d_manual_result.csv (provisional single-cycle) results.

using DelimitedFiles
using GLMakie
using Printf

base = @__DIR__
root = get(ENV, "RESULTS_ROOT", joinpath(base, "data"))

function read_parameter_file(path)
    d = Dict{String, String}()
    isfile(path) || return d
    for line in eachline(path)
        occursin('=', line) || continue
        key, value = strip.(split(line, '='; limit = 2))
        d[key] = value
    end
    d
end

function parse_or_nan(params, key)
    haskey(params, key) || return NaN
    try
        parse(Float64, params[key])
    catch
        NaN
    end
end

rows = Vector{Dict{String, Any}}()

# Cases live either directly under the results root (single-angle runs) or one
# level deeper inside a sweep directory (e.g. data/coarse_N64/<angle>degree).
function candidate_case_dirs(root)
    dirs = String[]
    for entry in sort(readdir(root; join = true))
        isdir(entry) || continue
        isfile(joinpath(entry, "case_parameters.txt")) && push!(dirs, entry)
        for sub in sort(readdir(entry; join = true))
            isdir(sub) || continue
            isfile(joinpath(sub, "case_parameters.txt")) && push!(dirs, sub)
        end
    end
    dirs
end

for entry in candidate_case_dirs(root)
    params = read_parameter_file(joinpath(entry, "case_parameters.txt"))

    # Production sweep only: exclude validation/low-Ra backup directories.
    ra = parse_or_nan(params, "Ra")
    (isfinite(ra) && abs(ra - 1.0e7) < 1e-9) || continue

    f_std = joinpath(entry, "task_d", "task_d_result.csv")
    f_man = joinpath(entry, "task_d_manual", "task_d_manual_result.csv")

    if isfile(f_std)
        data = readdlm(f_std, ',', '\n'; skipstart = 1)
        size(data, 1) >= 1 || continue
        r = data[1, :]
        # header: theta,Ra,Pr,drive_force_y,t_start,t_end,Vbulk,flow_area,Vabs,Uwind,regime,status
        push!(rows, Dict(
            "theta" => Float64(r[1]),
            "Vbulk" => Float64(r[7]),
            "Vabs" => Float64(r[9]),
            "Uwind" => Float64(r[10]),
            "t_start" => Float64(r[5]),
            "t_end" => Float64(r[6]),
            "regime" => string(r[11]),
            "status" => string(r[12]),
            "grid_n" => parse_or_nan(params, "N_wall"),
            "runtime" => parse_or_nan(params, "runtime_seconds"),
        ))
    elseif isfile(f_man)
        data = readdlm(f_man, ',', Float64, '\n'; skipstart = 1)
        size(data, 1) >= 1 || continue
        r = vec(data[1, :])
        # header: t_start,t_end,Vbulk,flow_area,Vabs,Uwind
        regime = "provisional_single_cycle"
        status = "provisional_single_cycle"
        sel = joinpath(entry, "task_c", "task_c_selected_interval.csv")
        if isfile(sel)
            seldata = readdlm(sel, ',', '\n'; skipstart = 1)
            if size(seldata, 1) >= 1 && size(seldata, 2) >= 6
                regime = string(seldata[1, 1])
                status = string(seldata[1, 6])
            end
        end
        push!(rows, Dict(
            "theta" => parse_or_nan(params, "theta"),
            "Vbulk" => Float64(r[3]),
            "Vabs" => Float64(r[5]),
            "Uwind" => Float64(r[6]),
            "t_start" => Float64(r[1]),
            "t_end" => Float64(r[2]),
            "regime" => regime,
            "status" => status,
            "grid_n" => parse_or_nan(params, "N_wall"),
            "runtime" => parse_or_nan(params, "runtime_seconds"),
        ))
    end
end

isempty(rows) && error("No completed Task-D results found under $root")

# Deduplicate by theta, preferring the standard task_d result over task_d_manual.
best = Dict{Float64, Dict{String, Any}}()
for r in rows
    θ = Float64(r["theta"])
    if !haskey(best, θ) || (r["status"] != "provisional_single_cycle" && best[θ]["status"] == "provisional_single_cycle")
        best[θ] = r
    end
end

θs = sort(collect(keys(best)))
n = length(θs)

outcsv = joinpath(root, "bulk_velocity_vs_angle.csv")
open(outcsv, "w") do io
    println(io, "theta,Vbulk,Vabs,Uwind,t_start,t_end,regime,status,grid_n,runtime")
    for θ in θs
        r = best[θ]
        println(io,
            r["theta"], ',',
            r["Vbulk"], ',',
            r["Vabs"], ',',
            r["Uwind"], ',',
            r["t_start"], ',',
            r["t_end"], ',',
            r["regime"], ',',
            r["status"], ',',
            r["grid_n"], ',',
            r["runtime"],
        )
    end
end

Vbulk = [Float64(best[θ]["Vbulk"]) for θ in θs]
Vbulk_abs = abs.(Vbulk)
imax = argmax(Vbulk_abs)
θopt = θs[imax]
Vopt = Vbulk_abs[imax]

fig = Figure(size = (900, 650))
ax = Axis(fig[1, 1];
    title = "Task D: |Vbulk| versus chimney angle",
    xlabel = "theta (degree)",
    ylabel = "|Vbulk|",
)
lines!(ax, θs, Vbulk_abs)
scatter!(ax, θs, Vbulk_abs)
scatter!(ax, [θopt], [Vopt]; markersize = 18)
save(joinpath(root, "Vbulk_vs_angle.png"), fig)

open(joinpath(root, "optimal_angle.txt"), "w") do io
    println(io, "theta_opt = ", θopt)
    println(io, "Vbulk_abs_max = ", Vopt)
    println(io, "Vbulk_at_opt = ", Vbulk[imax])
end

@printf("coarse optimum angle (by |Vbulk|) = %.6g deg\n", θopt)
@printf("max |Vbulk|                       = %.8e\n", Vopt)
@printf("Vbulk at optimum                  = %.8e\n", Vbulk[imax])
println("Wrote: ", outcsv)
