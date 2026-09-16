# run_cd_case.jl
# -----------------------------------------------------------------------------
# Production CFD data generator for Task C / Task D.
#
# Coordinate convention:
#   x : wall-normal, hot wall -> cold wall
#   y : vertical / streamwise, periodic
#   theta = 90 deg : chimney parallel to y
#   T = T(x), decreasing from +0.5 to -0.5
#   v = u_y, and the Task-C profile is vbar(x,t) = <u_y>_y
#
# This runner deliberately does NOT impose zero mass flux. Task D must measure
# the bulk velocity produced by the physics, not a prescribed zero value.
#
# Optional full-field snapshots (for report figures) are controlled by
# SNAPSHOT_TIMES, e.g. SNAPSHOT_TIMES="20,60,100". They reuse the generic
# field writer from the validation module and do not affect the solution.
# -----------------------------------------------------------------------------

using CUDA
using CUDSS
using IncompressibleNavierStokes
using Printf

base = @__DIR__
include(joinpath(base, "..", "shared", "common", "cd_grid.jl"))
include(joinpath(base, "..", "shared", "common", "cd_solver.jl"))
include(joinpath(base, "..", "shared", "common", "cd_tools.jl"))
include(joinpath(base, "..", "shared", "common", "cd_validation.jl"))   # generic field snapshots

# Parameters can be overridden from the shell, e.g.
#   THETA=90 T_END=400 DRIVE_FORCE_Y=0.02 julia run_cd_case.jl
const THETA = parse(Float64, get(ENV, "THETA", "90.0"))
const T_END = parse(Float64, get(ENV, "T_END", "100.0"))
const RA = parse(Float64, get(ENV, "RA", "1000.0"))
const PR = parse(Float64, get(ENV, "PR", "0.71"))
const DRIVE_FORCE_Y = parse(Float64, get(ENV, "DRIVE_FORCE_Y", "0.0"))

const N_WALL = parse(Int, get(ENV, "N_WALL", "128"))
const GAP = parse(Float64, get(ENV, "GAP", "1.0"))
const PERIODIC_LENGTH = parse(Float64, get(ENV, "PERIODIC_LENGTH", "2.0"))
const CLUSTERING = parse(Float64, get(ENV, "CLUSTERING", "1.2"))
const SAMPLE_DT = parse(Float64, get(ENV, "SAMPLE_DT", "0.2"))
const MONITOR_DT = parse(Float64, get(ENV, "MONITOR_DT", "5.0"))

# Optional field snapshots, e.g. SNAPSHOT_TIMES="20,60,100". Empty = disabled.
const SNAPSHOT_TIMES = let s = strip(get(ENV, "SNAPSHOT_TIMES", ""))
    isempty(s) ? Float64[] : [parse(Float64, strip(tok)) for tok in split(s, ',')]
end

# Output root. Overridable so that visualisation re-runs (different grid /
# T_END) never overwrite the production data used by Task C / Task D.
# Layout: production cases -> main_task_cd/data. A low-Ra validation run must
# point elsewhere, e.g. $env:RESULTS_ROOT = "..\validation_lowRa\data".
const RESULTS_ROOT = get(ENV, "RESULTS_ROOT", joinpath(base, "data"))
label = isapprox(THETA, round(Int, THETA); atol = 1e-10) ? string(round(Int, THETA)) : replace(string(THETA), "." => "p")
OUTDIR = joinpath(RESULTS_ROOT, "$(label)degree")
mkpath(OUTDIR)
raw_path = joinpath(OUTDIR, "vbar_xt_raw.csv")
monitor_path = joinpath(OUTDIR, "runtime_monitor.csv")

CUDA.functional() || error("CUDA is not functional")
backend = CUDA.CUDABackend()
println("GPU: ", CUDA.name(CUDA.device()))
println("CPU threads: ", Sys.CPU_THREADS)

coeff = cd_dimensionless_coefficients(RA, PR)

setup = build_cd_setup(
    N_WALL;
    gap = GAP,
    periodic_length = PERIODIC_LENGTH,
    clustering = CLUSTERING,
    backend = backend,
)

psolver = default_psolver(setup)
start = build_cd_initial_state(
    setup,
    psolver;
    gap = GAP,
    periodic_length = PERIODIC_LENGTH,
)

rec = make_task_c_recorder(
    setup;
    sample_interval = SAMPLE_DT,
    start_time = 0.0,
    output_path = raw_path,
)

monitor = make_runtime_monitor(
    setup;
    theta = THETA,
    t_end = T_END,
    report_interval = MONITOR_DT,
    output_path = monitor_path,
    xlims = (Float64(setup.xlims[1][1]), Float64(setup.xlims[1][2])),
)

# Generic full-field snapshot recorder (temperature + velocity on the pressure
# grid); a no-op when SNAPSHOT_TIMES is empty.
field_rec = make_huang_field_recorder(
    setup;
    output_dir = OUTDIR,
    snapshot_times = SNAPSHOT_TIMES,
)

println("------------------------------------------------------------")
@printf("theta             = %.6g deg\n", THETA)
@printf("Ra                = %.6e\n", RA)
@printf("Pr                = %.6g\n", PR)
@printf("nu*               = %.6e\n", coeff.viscosity)
@printf("kappa*            = %.6e\n", coeff.conductivity)
@printf("gap Lx            = %.6g\n", GAP)
@printf("periodic Ly       = %.6g\n", PERIODIC_LENGTH)
@printf("drive_force_y     = %.6e\n", DRIVE_FORCE_Y)
@printf("sample interval   = %.6g\n", SAMPLE_DT)
@printf("monitor interval  = %.6g\n", MONITOR_DT)
@printf("grid Nx           = %d\n", N_WALL)
@printf("grid Ny           = %d\n", length(setup.xp[2]))
@printf("total cells (x*y) = %d\n", prod(setup.N))
@printf("field snapshots   = %s\n",
        isempty(SNAPSHOT_TIMES) ? "off" : join(string.(SNAPSHOT_TIMES), ", "))
@printf("output directory  = %s\n", OUTDIR)
println("zero mass flux    = OFF (required for Task D bulk velocity)")
println("------------------------------------------------------------")

t0 = time()
state, outputs = solve_unsteady(;
    force! = boussinesq_cd!,
    setup,
    start,
    tlims = (0.0, T_END),
    psolver,
    params = (;
        viscosity = Float64(coeff.viscosity),
        conductivity = Float64(coeff.conductivity),
        gravity = 1.0,
        theta = Float64(THETA),
        drive_force_y = Float64(DRIVE_FORCE_Y),
        dodissipation = false,
    ),
    processors = (;
        rec,
        monitor,
        field_rec,
        log = timelogger(; nupdate = 2000),
    ),
)

elapsed = time() - t0
println("Finished: t = ", state.t, " in ", round(elapsed; digits = 1), " s")

isnothing(rec.xref[]) && error("Task-C recorder did not collect any data")
rec.close_io!()
monitor.close_io!()
write_vbar_long_csv(raw_path, rec.times, rec.xref[], rec.profiles)
println("Wrote: ", raw_path)
println("Monitor: ", monitor_path)
for p in field_rec.written
    println("Snapshot: ", p)
end

# Save the parameters needed by Task C/D and the report.
open(joinpath(OUTDIR, "case_parameters.txt"), "w") do io
    println(io, "theta = ", THETA)
    println(io, "Ra = ", RA)
    println(io, "Pr = ", PR)
    println(io, "viscosity = ", coeff.viscosity)
    println(io, "conductivity = ", coeff.conductivity)
    println(io, "gap_x = ", GAP)
    println(io, "periodic_length_y = ", PERIODIC_LENGTH)
    println(io, "clustering = ", CLUSTERING)
    println(io, "N_wall = ", N_WALL)
    println(io, "drive_force_y = ", DRIVE_FORCE_Y)
    println(io, "sample_interval_target = ", SAMPLE_DT)
    println(io, "x_min = ", setup.xlims[1][1])
    println(io, "x_max = ", setup.xlims[1][2])
    println(io, "y_min = ", setup.xlims[2][1])
    println(io, "y_max = ", setup.xlims[2][2])
    println(io, "zero_mass_flux = false")
    println(io, "snapshot_times = ", isempty(SNAPSHOT_TIMES) ? "none" : join(string.(SNAPSHOT_TIMES), ","))
    println(io, "runtime_seconds = ", elapsed)
    println(io, "t_end_actual = ", state.t)
end

println("Next: run task_c.jl on this case, then task_d.jl if Task C passes.")
