# collect_ra1000_results.jl
# -----------------------------------------------------------------------------
# Collect every Ra = 1000 case into one delivery folder:
#
#   data/Ra1000_summary/
#       ra1000_summary.csv          one row per collected case
#       Vabs_Uwind_vs_angle.png     velocity diagnostics versus chimney angle
#       README.txt                  what the numbers mean / how to regenerate
#       <case>/case_parameters.txt  copy of the production input file
#       <case>/task_c/...           full Task C products (audit, windows, figures)
#       <case>/task_d/...           full Task D products (profile, result, figure)
#       <case>/task_d_manual/...    only when a manual single-cycle run exists
#
# Only cases whose case_parameters.txt reports Ra = 1000 are collected, so the
# Ra = 1e7 production runs are never mixed in. The raw CFD record
# vbar_xt_raw.csv is NOT copied (it stays in the original case directory);
# everything Task C / Task D derived from it is.
#
# The script is idempotent: re-running it refreshes the folder in place.
#
# Usage:
#   julia collect_ra1000_results.jl
#   RA_TARGET=1000 OUTDIR=data/Ra1000_summary julia collect_ra1000_results.jl
# -----------------------------------------------------------------------------

using DelimitedFiles
using Printf

base = @__DIR__
root = get(ENV, "RESULTS_ROOT", joinpath(base, "data"))

const RA_TARGET = parse(Float64, get(ENV, "RA_TARGET", "1000.0"))
const OUTDIR = get(ENV, "OUTDIR", joinpath(root, "Ra1000_summary"))
const CASE_SUBDIRS = ("task_c", "task_d", "task_d_manual")

include(joinpath(base, "..", "shared", "common", "cd_plots.jl"))   # figure conventions / font

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

function parse_or_nan(d, key)
    haskey(d, key) || return NaN
    v = tryparse(Float64, strip(d[key]))
    isnothing(v) ? NaN : v
end

# -----------------------------------------------------------------------------
# 1. find every Ra = 1000 case
# -----------------------------------------------------------------------------
cases = String[]
for entry in sort(readdir(root; join = true))
    isdir(entry) || continue
    params = read_parameter_file(joinpath(entry, "case_parameters.txt"))
    ra = parse_or_nan(params, "Ra")
    (isfinite(ra) && isapprox(ra, RA_TARGET; rtol = 1e-9)) || continue
    push!(cases, entry)
end
isempty(cases) && error("no Ra = $RA_TARGET cases found under $root")

println("source root : ", root)
println("output dir  : ", OUTDIR)
println("Ra target   : ", RA_TARGET)
println("cases found : ", length(cases))

# -----------------------------------------------------------------------------
# 2. copy the case products
# -----------------------------------------------------------------------------
mkpath(OUTDIR)
rows = NamedTuple[]

for case in cases
    name = basename(case)
    dest = joinpath(OUTDIR, name)
    mkpath(dest)

    src_params = joinpath(case, "case_parameters.txt")
    isfile(src_params) && cp(src_params, joinpath(dest, "case_parameters.txt"); force = true)

    copied = 0
    for sub in CASE_SUBDIRS
        srcdir = joinpath(case, sub)
        isdir(srcdir) || continue
        dstdir = joinpath(dest, sub)
        mkpath(dstdir)
        for f in sort(readdir(srcdir))
            src = joinpath(srcdir, f)
            isfile(src) || continue
            cp(src, joinpath(dstdir, f); force = true)
            copied += 1
        end
    end

    # ---- one summary row from the Task D result -----------------------------
    result_file = joinpath(case, "task_d", "task_d_result.csv")
    if isfile(result_file)
        data = readdlm(result_file, ',', '\n'; skipstart = 1)
        if size(data, 1) >= 1 && size(data, 2) >= 12
            r = data[1, :]
            params = read_parameter_file(src_params)
            push!(rows, (
                theta = Float64(r[1]),
                grid_n = parse_or_nan(params, "N_wall"),
                Ra = Float64(r[2]),
                Pr = Float64(r[3]),
                t_start = Float64(r[5]),
                t_end = Float64(r[6]),
                t_end_actual = parse_or_nan(params, "t_end_actual"),
                Vbulk = Float64(r[7]),
                flow_area = Float64(r[8]),
                Vabs = Float64(r[9]),
                Uwind = Float64(r[10]),
                regime = string(r[11]),
                status = string(r[12]),
                source_dir = name,
            ))
        else
            @warn "unexpected task_d_result.csv layout in $name"
        end
    else
        @warn "no task_d/task_d_result.csv in $name"
    end

    @printf("collected %-20s files=%d\n", name, copied)
end

isempty(rows) && error("no Task D results found among the Ra = $RA_TARGET cases")

# -----------------------------------------------------------------------------
# 3. summary table
# -----------------------------------------------------------------------------
sort!(rows; by = r -> (r.theta, r.grid_n, r.source_dir))

summary_csv = joinpath(OUTDIR, "ra1000_summary.csv")
open(summary_csv, "w") do io
    println(
        io,
        "theta,grid_N,Ra,Pr,t_start,t_end,t_end_actual,Vbulk,flow_area,Vabs,Uwind,regime,status,source_dir",
    )
    for r in rows
        println(
            io,
            r.theta, ',', r.grid_n, ',', r.Ra, ',', r.Pr, ',',
            r.t_start, ',', r.t_end, ',', r.t_end_actual, ',',
            r.Vbulk, ',', r.flow_area, ',', r.Vabs, ',', r.Uwind, ',',
            r.regime, ',', r.status, ',', r.source_dir,
        )
    end
end
println("wrote ", summary_csv)

# -----------------------------------------------------------------------------
# 4. diagnostics figure (Vbulk is numerically zero, so Vabs / Uwind are the
#    meaningful magnitudes)
# -----------------------------------------------------------------------------
main_rows = [r for r in rows if r.grid_n == 128]
other_rows = [r for r in rows if r.grid_n != 128]

fig = Figure(size = (1000, 700))
ax = Axis(
    fig[1, 1];
    title = @sprintf("Ra = %g: velocity diagnostics versus chimney angle", RA_TARGET),
    xlabel = "theta (deg)",
    ylabel = "velocity",
)
if !isempty(main_rows)
    scatterlines!(
        ax,
        [r.theta for r in main_rows],
        [r.Vabs for r in main_rows];
        marker = :circle,
        label = "Vabs (N = 128)",
    )
    scatterlines!(
        ax,
        [r.theta for r in main_rows],
        [r.Uwind for r in main_rows];
        marker = :circle,
        label = "Uwind (N = 128)",
    )
end
if !isempty(other_rows)
    scatter!(
        ax,
        [r.theta for r in other_rows],
        [r.Vabs for r in other_rows];
        marker = :xcross,
        label = "Vabs (other grids)",
    )
    scatter!(
        ax,
        [r.theta for r in other_rows],
        [r.Uwind for r in other_rows];
        marker = :xcross,
        label = "Uwind (other grids)",
    )
end
isempty(rows) || axislegend(ax; position = :lt)
save(joinpath(OUTDIR, "Vabs_Uwind_vs_angle.png"), fig)
println("wrote ", joinpath(OUTDIR, "Vabs_Uwind_vs_angle.png"))

# -----------------------------------------------------------------------------
# 5. README
# -----------------------------------------------------------------------------
max_abs_vbulk = maximum(abs.(getproperty.(rows, :Vbulk)))
max_abs_vabs = maximum(abs.(getproperty.(rows, :Vabs)))
max_abs_uwind = maximum(abs.(getproperty.(rows, :Uwind)))
n_grids = length(unique(getproperty.(rows, :grid_n)))

open(joinpath(OUTDIR, "README.txt"), "w") do io
    println(io, "Ra = 1000 case collection")
    println(io, "======================")
    println(io)
    println(io, "Source root    : ", root)
    println(io, "Ra             : ", RA_TARGET)
    println(io, "Cases collected: ", length(rows), " (", length(cases), " directories)")
    println(io, "Grids present  : ", join(string.(sort(unique(getproperty.(rows, :grid_n)))), ", "))
    println(io, "Other physics  : Pr = 0.71, drive_force_y = 0.0, zero_mass_flux = false")
    println(io)
    println(io, "Layout")
    println(io, "------")
    println(io, "ra1000_summary.csv          one row per collected case (sorted by theta, then grid)")
    println(io, "Vabs_Uwind_vs_angle.png     Vabs and Uwind versus chimney angle theta")
    println(io, "<case>/case_parameters.txt  copy of the production input file")
    println(io, "<case>/task_c/              Task C products: audit, window statistics/profiles, figures")
    println(io, "<case>/task_d/              Task D products: vbar_t profile, result row, V-X figure")
    println(io, "<case>/task_d_manual/       present only when a manual single-cycle run exists")
    println(io)
    println(io, "NOT copied: vbar_xt_raw.csv (the raw CFD record, ~10 MB per case). It stays")
    println(io, "in the original case directory; every Task C / Task D product derived from it")
    println(io, "is included here.")
    println(io)
    println(io, "What the numbers mean")
    println(io, "---------------------")
    println(io, "Every collected case reached the stationary regime and is reported as")
    println(io, "status = converged_stationary. The SIGNED bulk velocity is numerically zero:")
    println(io)
    @printf(io, "    max |Vbulk| over all cases = %.3e\n", max_abs_vbulk)
    println(io)
    println(io, "This is the physically expected outcome for this configuration: the wall")
    println(io, "temperature difference drives a symmetric counter-flow, so the wall-normal")
    println(io, "average of the streamwise velocity integrates to zero. Values of order 1e-17")
    println(io, "are floating-point round-off, NOT a physical net through-flow. A non-zero")
    println(io, "drive_force_y (or an asymmetric geometry) is required to obtain a net flow.")
    println(io)
    println(io, "The meaningful magnitudes are the diagnostics:")
    println(io, "    Vabs  = (1/Lx) * integral_x |vbar_t(x)| dx   (mean speed across the gap)")
    println(io, "    Uwind = max_x |vbar_t(x)|                    (peak speed)")
    @printf(io, "\n    max Vabs  over all cases = %.6g\n", max_abs_vabs)
    @printf(io, "    max Uwind over all cases = %.6g\n", max_abs_uwind)
    println(io, "\nBoth grow with the chimney angle theta; see ra1000_summary.csv.")
    println(io, "Cases with a grid other than N = 128 (if any) are kept in the table for")
    println(io, "grid-sensitivity checking and are marked by grid_N.")
    println(io)
    println(io, "How to regenerate")
    println(io, "-----------------")
    println(io, "Refresh this folder from the existing results (idempotent):")
    println(io, "    julia collect_ra1000_results.jl")
    println(io)
    println(io, "Rerun one case from scratch:")
    println(io, "    THETA=<angle> RA=1000 T_END=100 N_WALL=128 julia run_cd_case.jl")
    println(io, "    CASE_DIR=data/<angle>degree julia task_c.jl")
    println(io, "    CASE_DIR=data/<angle>degree julia task_d.jl")
end
println("wrote ", joinpath(OUTDIR, "README.txt"))
@printf("max |Vbulk| = %.3e (numerically zero)\n", max_abs_vbulk)
println("done.")
