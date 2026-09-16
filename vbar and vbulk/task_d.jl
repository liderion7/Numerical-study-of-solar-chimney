# task_d.jl
# -----------------------------------------------------------------------------
# Task D: use the Task-C screened data to calculate the final time-averaged
# velocity profile and the signed bulk velocity.
#
# Definitions:
#   vbar(x,t)   = (1/Ly) integral_y v(x,y,t) dy
#   vbar_t(x)   = (1/T)  integral_t vbar(x,t) dt
#   V_bulk      = (1/Lx) integral_x vbar_t(x) dx
#
# Both integrations are trapezoidal and use the actual time/space coordinates.
# -----------------------------------------------------------------------------

using IncompressibleNavierStokes
using GLMakie
using Printf
using DelimitedFiles

base = @__DIR__
include(joinpath(base, "..", "shared", "common", "cd_tools.jl"))

# Default: the most complete production 90-degree case (main_task_cd/data).
const CASE_DIR = get(ENV, "CASE_DIR", joinpath(base, "data", "90degree_Ra1e7_T550"))
const TASKC_DIR = joinpath(CASE_DIR, "task_c")
const INPUT = get(ENV, "TASKC_FILE", joinpath(TASKC_DIR, "task_c_clean_vbar.csv"))
const PARAM_FILE = joinpath(CASE_DIR, "case_parameters.txt")
const OUTDIR = joinpath(CASE_DIR, "task_d")
mkpath(OUTDIR)

isfile(INPUT) || error("Task C clean file not found: $INPUT. Run task_c.jl first and ensure Task D ready = true.")

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

params = read_parameter_file(PARAM_FILE)
xmin = parse(Float64, get(params, "x_min", "0.0"))
xmax = parse(Float64, get(params, "x_max", "1.0"))
xlims = (xmin, xmax)

times, x, V = read_vbar_long_csv(INPUT)
audit = audit_task_c_data(times, x, V)
audit.all_finite || error("Task D input contains NaN/Inf")
audit.time_strictly_increasing || error("Task D time coordinate is not strictly increasing")
audit.x_strictly_increasing || error("Task D x coordinate is not strictly increasing")

# Read the Task-C decision (regime/status) if it exists.
regime = "unknown"
status = "unknown"
selected_file = joinpath(TASKC_DIR, "task_c_selected_interval.csv")
if isfile(selected_file)
    sel = readdlm(selected_file, ',', '\n'; skipstart = 1)
    if size(sel, 1) >= 1 && size(sel, 2) >= 6
        regime = string(sel[1, 1])
        status = string(sel[1, 6])
    end
end

# Task C already clipped the clean file to the approved interval.
t_start = times[1]
t_end = times[end]

# Step 1: true time-weighted average at every x.
vbar_t = time_weighted_profile(times, V, t_start, t_end)

# Step 2: true spatial integral/average across the x direction.
diag = task_d_velocity_diagnostics(x, vbar_t, xlims)

# Full wall profile for report/plot, including v=0 at both no-slip walls.
xfull, vfull = full_wall_profile(x, vbar_t, xlims)
flow_area = trapz1(xfull, vfull)  # nondimensional volume-flux proxy per unit depth

# Save final profile.
open(joinpath(OUTDIR, "vbar_t_profile.csv"), "w") do io
    println(io, "x,vbar_t")
    for i in eachindex(xfull)
        println(io, xfull[i], ',', vfull[i])
    end
end

# Save one-row summary for later angle sweep.
theta = parse(Float64, get(params, "theta", "NaN"))
Ra = parse(Float64, get(params, "Ra", "NaN"))
Pr = parse(Float64, get(params, "Pr", "NaN"))
drive_force_y = parse(Float64, get(params, "drive_force_y", "NaN"))

open(joinpath(OUTDIR, "task_d_result.csv"), "w") do io
    println(io, "theta,Ra,Pr,drive_force_y,t_start,t_end,Vbulk,flow_area,Vabs,Uwind,regime,status")
    println(
        io,
        theta, ',', Ra, ',', Pr, ',', drive_force_y, ',',
        t_start, ',', t_end, ',', diag.Vbulk, ',', flow_area, ',',
        diag.Vabs, ',', diag.Uwind, ',', regime, ',', status,
    )
end

open(joinpath(OUTDIR, "task_d_summary.txt"), "w") do io
    println(io, "Task D summary")
    println(io, "==============")
    println(io, "theta = ", theta)
    println(io, "Ra = ", Ra)
    println(io, "Pr = ", Pr)
    println(io, "drive_force_y = ", drive_force_y)
    println(io, "approved averaging interval = [", t_start, ", ", t_end, "]")
    println(io, "regime = ", regime)
    println(io, "status = ", status)
    println(io)
    println(io, "Primary Task-D result")
    println(io, "Vbulk = ", diag.Vbulk)
    println(io, "flow area integral = ", flow_area)
    println(io)
    println(io, "Additional diagnostics (not replacements for Vbulk)")
    println(io, "Vabs = ", diag.Vabs)
    println(io, "Uwind = ", diag.Uwind)
    println(io)
    if abs(diag.Vbulk) < 1e-6 && abs(drive_force_y) < 1e-14
        println(io, "NOTE: Vbulk is approximately zero while drive_force_y is zero.")
        println(io, "This can be physically consistent for symmetric natural-convection counterflow.")
        println(io, "A non-zero pressure-gradient/external drive is required if the final chimney task expects net through-flow.")
    end
    if status == "provisional_single_cycle"
        println(io, "NOTE: this result is a provisional single-cycle estimate.")
        println(io, "It is NOT a fully statistically converged value.")
    end
end

# Report-ready V-X plot.
# The shared module is loaded here, next to the plot it serves, so that the
# Task C / Task D axis convention (horizontal = position x running 1 to 0,
# vertical = velocity) has a single definition.
include(joinpath(base, "..", "shared", "common", "cd_plots.jl"))

plot_task_d_profile(
    OUTDIR,
    "V_X_time_averaged.png",
    xfull,
    vfull;
    title = "Task D: time-averaged V-X profile",
)

println("------------------------------------------------------------")
@printf("Task D interval: [%.6g, %.6g]\n", t_start, t_end)
@printf("Vbulk            = % .8e\n", diag.Vbulk)
@printf("flow area        = % .8e\n", flow_area)
@printf("Vabs diagnostic  = % .8e\n", diag.Vabs)
@printf("Uwind diagnostic = % .8e\n", diag.Uwind)
println("Outputs: ", OUTDIR)
println("------------------------------------------------------------")
