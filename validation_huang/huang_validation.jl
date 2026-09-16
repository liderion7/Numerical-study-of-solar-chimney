# huang_validation.jl
# -----------------------------------------------------------------------------
# Direct benchmark runner for Huang et al. (2022), Table 1:
#   plate BC = no-slip (NS)
#   sidewall BC = periodic (PD)
#   Pr = 4.3
#   aspect ratio Gamma = L/H = 2
#
# Put this file in the project root, and put cd_validation.jl in common/.
# The existing files common/cd_grid.jl and common/cd_solver.jl are reused.
#
# Recommended first run (PowerShell):
#   $env:RA="1e6"
#   $env:N_WALL="64"
#   julia huang_validation.jl
#
# Grid-convergence follow-up:
#   $env:N_WALL="128"
#   julia huang_validation.jl
#
# Default timing follows the paper logic:
#   0 -> 200      transient / spin-up
#   200 -> 1200   statistics (~1000 free-fall time units)
# -----------------------------------------------------------------------------

using CUDA
using CUDSS
using IncompressibleNavierStokes
using Printf

base = @__DIR__
include(joinpath(base, "..", "shared", "common", "cd_grid.jl"))
include(joinpath(base, "..", "shared", "common", "cd_solver.jl"))
include(joinpath(base, "..", "shared", "common", "cd_validation.jl"))

# ------------------------------- parameters ----------------------------------
const RA = parse(Float64, get(ENV, "RA", "1.0e6"))
const PR = parse(Float64, get(ENV, "PR", "4.3"))

# Huang-compatible geometry. Keep these values for direct comparison.
const GAP = parse(Float64, get(ENV, "GAP", "1.0"))
const PERIODIC_LENGTH = parse(Float64, get(ENV, "PERIODIC_LENGTH", "2.0"))
const THETA = 0.0
const DRIVE_FORCE_Y = 0.0

# Current finite-volume grid. This is NOT Huang's spectral-element N=7.
const N_WALL = parse(Int, get(ENV, "N_WALL", "64"))
const NY_ENV = get(ENV, "NY", "")
const NY = isempty(strip(NY_ENV)) ? nothing : parse(Int, NY_ENV)
const CLUSTERING = parse(Float64, get(ENV, "CLUSTERING", "1.2"))

# Time settings.
const AVG_START = parse(Float64, get(ENV, "AVG_START", "200.0"))
const T_END = parse(Float64, get(ENV, "T_END", "1200.0"))
const SAMPLE_DT = parse(Float64, get(ENV, "SAMPLE_DT", "0.5"))
const MONITOR_DT = parse(Float64, get(ENV, "MONITOR_DT", "20.0"))

# Initial condition. "perturbed" is robust; use "zero" for the exact paper text.
const INIT_MODE = get(ENV, "INIT_MODE", "perturbed")
const INIT_EPS = parse(Float64, get(ENV, "INIT_EPS", "1e-6"))
# Number of horizontal perturbation wavelengths across the periodic length.
# Selects the convection-roll pattern: 1 = original setup (wavelength L).
const PERT_WAVES = parse(Int, get(ENV, "PERT_WAVES", "1"))

# ------------------------------ safety checks --------------------------------
RA > 0 || error("RA must be positive")
N_WALL > 0 || error("N_WALL must be positive")
T_END > AVG_START || error("T_END must be greater than AVG_START")
SAMPLE_DT > 0 || error("SAMPLE_DT must be positive")
MONITOR_DT > 0 || error("MONITOR_DT must be positive")
PERT_WAVES >= 1 || error("PERT_WAVES must be >= 1")

isapprox(PR, 4.3; rtol = 0, atol = 1e-12) || error(
    "Direct Huang comparison requires Pr = 4.3. Current PR=$PR",
)
isapprox(GAP, 1.0; rtol = 0, atol = 1e-12) || error(
    "Direct Huang comparison requires nondimensional cell height GAP=1.0. Current GAP=$GAP",
)
isapprox(PERIODIC_LENGTH, 2.0; rtol = 0, atol = 1e-12) || error(
    "Direct Huang comparison requires aspect ratio Gamma=L/H=2, so PERIODIC_LENGTH=2.0 when GAP=1. Current value=$PERIODIC_LENGTH",
)

# ------------------------------- output path ---------------------------------
function safe_tag(x::Real)
    s = @sprintf("%.4g", Float64(x))
    replace(s, "." => "p", "+" => "", "-" => "m")
end

ra_tag = safe_tag(RA)
pr_tag = safe_tag(PR)
init_tag = lowercase(replace(strip(INIT_MODE), " " => "_"))
case_name = "Huang_NS_PD_Ra$(ra_tag)_Pr$(pr_tag)_N$(N_WALL)_$(init_tag)_W$(PERT_WAVES)"
OUTDIR = joinpath(base, "data", case_name)
mkpath(OUTDIR)

stats_path = joinpath(OUTDIR, "huang_statistics_timeseries.csv")

# --------------------------------- setup --------------------------------------
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
    ny = NY,
    backend = backend,
)

psolver = default_psolver(setup)
start = build_huang_initial_state(
    setup,
    psolver;
    gap = GAP,
    periodic_length = PERIODIC_LENGTH,
    mode = INIT_MODE,
    epsilon = INIT_EPS,
    perturbation_waves = PERT_WAVES,
)

stats_rec = make_huang_statistics_recorder(
    setup;
    Ra = RA,
    Pr = PR,
    H = GAP,
    sample_interval = SAMPLE_DT,
    start_time = AVG_START,
    output_path = stats_path,
)

monitor = make_huang_runtime_monitor(
    setup;
    Ra = RA,
    Pr = PR,
    H = GAP,
    t_end = T_END,
    report_interval = MONITOR_DT,
)

# Optional full-field snapshots, e.g. SNAPSHOT_TIMES="100,200,400". Empty = off.
const SNAPSHOT_TIMES = let s = strip(get(ENV, "SNAPSHOT_TIMES", ""))
    isempty(s) ? Float64[] : [parse(Float64, strip(tok)) for tok in split(s, ',')]
end

field_rec = make_huang_field_recorder(
    setup;
    output_dir = OUTDIR,
    snapshot_times = SNAPSHOT_TIMES,
)

ref = huang_reference_ns_pd(RA)

println("------------------------------------------------------------")
println("Huang et al. (2022) validation: NS plates / PD sidewalls")
println("------------------------------------------------------------")
@printf("Ra                = %.6e\n", RA)
@printf("Pr                = %.6g\n", PR)
@printf("theta             = %.1f deg (buoyancy along code +x)\n", THETA)
@printf("nu*               = %.8e\n", coeff.viscosity)
@printf("kappa*            = %.8e\n", coeff.conductivity)
@printf("H = gap           = %.6g\n", GAP)
@printf("L periodic        = %.6g\n", PERIODIC_LENGTH)
@printf("aspect ratio L/H  = %.6g\n", PERIODIC_LENGTH / GAP)
@printf("N_wall            = %d\n", N_WALL)
@printf("Ny pressure grid  = %d\n", length(setup.xp[2]))
@printf("clustering        = %.6g\n", CLUSTERING)
@printf("transient to      = %.6g\n", AVG_START)
@printf("T_end             = %.6g\n", T_END)
@printf("target averaging  = %.6g\n", T_END - AVG_START)
@printf("sample interval   = %.6g\n", SAMPLE_DT)
println("initial mode      = ", INIT_MODE)
@printf("initial epsilon   = %.3e\n", INIT_EPS)
@printf("perturbation waves= %d\n", PERT_WAVES)
println("output directory  = ", OUTDIR)
@printf("field snapshots   = %s\n",
        isempty(SNAPSHOT_TIMES) ? "off" : join(string.(SNAPSHOT_TIMES), ", "))
println("coordinate map    = code x -> paper vertical y; code y -> paper horizontal x")
if ref === nothing
    println("reference row     = none stored for this Ra")
else
    @printf("Huang target      = Nu %.3f, Re %.3f, Rex %.3f, Rey %.3f, tavg %.0f\n",
            ref.Nu, ref.Re, ref.Rex, ref.Rey, ref.tavg)
end
println("------------------------------------------------------------")

# --------------------------------- solve --------------------------------------
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
        theta = THETA,
        drive_force_y = DRIVE_FORCE_Y,
        dodissipation = false,
    ),
    processors = (;
        stats_rec,
        field_rec,
        monitor,
        log = timelogger(; nupdate = 2000),
    ),
)

elapsed = time() - t0
stats_rec.close_io!()

println("Finished CFD: t = ", state.t, " in ", round(elapsed; digits = 1), " s")
println("Statistics time series: ", stats_path)

# ------------------------------- final stats ---------------------------------
stats = finalize_huang_statistics(stats_rec; Ra = RA, Pr = PR, H = GAP)
written = write_huang_summary(
    OUTDIR,
    stats;
    Ra = RA,
    Pr = PR,
    N_wall = N_WALL,
    Ny = length(setup.xp[2]),
    init_mode = INIT_MODE,
    init_epsilon = INIT_EPS,
)

println()
println("==================== VALIDATION RESULT ====================")
@printf("Nu  = %.6f\n", stats.Nu)
@printf("Re  = %.6f\n", stats.Re)
@printf("Rex = %.6f\n", stats.Rex)
@printf("Rey = %.6f\n", stats.Rey)
@printf("sqrt(Rex^2 + Rey^2) = %.6f\n", stats.Re_identity)
@printf("actual t_avg = %.6f (%d samples)\n", stats.tavg, stats.nsamples)

if written.ref !== nothing
    r = written.ref
    @printf("\n%-8s %12s %12s %12s\n", "quantity", "CFD", "Huang", "abs err %")
    @printf("%-8s %12.5f %12.5f %12.3f\n", "Nu",  stats.Nu,  r.Nu,  abs(100*(stats.Nu-r.Nu)/r.Nu))
    @printf("%-8s %12.5f %12.5f %12.3f\n", "Re",  stats.Re,  r.Re,  abs(100*(stats.Re-r.Re)/r.Re))
    @printf("%-8s %12.5f %12.5f %12.3f\n", "Rex", stats.Rex, r.Rex, abs(100*(stats.Rex-r.Rex)/r.Rex))
    @printf("%-8s %12.5f %12.5f %12.3f\n", "Rey", stats.Rey, r.Rey, abs(100*(stats.Rey-r.Rey)/r.Rey))
end

println("===========================================================")
println("Summary: ", written.txt)
println("CSV:     ", written.csv)

# Save run metadata separately for reproducibility.
open(joinpath(OUTDIR, "huang_case_parameters.txt"), "w") do io
    println(io, "Ra = ", RA)
    println(io, "Pr = ", PR)
    println(io, "theta = ", THETA)
    println(io, "viscosity = ", coeff.viscosity)
    println(io, "conductivity = ", coeff.conductivity)
    println(io, "gap_H = ", GAP)
    println(io, "periodic_length_L = ", PERIODIC_LENGTH)
    println(io, "aspect_ratio = ", PERIODIC_LENGTH / GAP)
    println(io, "plate_BC = NS")
    println(io, "sidewall_BC = PD")
    println(io, "N_wall = ", N_WALL)
    println(io, "Ny_pressure = ", length(setup.xp[2]))
    println(io, "clustering = ", CLUSTERING)
    println(io, "avg_start = ", AVG_START)
    println(io, "t_end_requested = ", T_END)
    println(io, "t_end_actual = ", state.t)
    println(io, "sample_dt = ", SAMPLE_DT)
    println(io, "initial_mode = ", INIT_MODE)
    println(io, "initial_epsilon = ", INIT_EPS)
    println(io, "runtime_seconds = ", elapsed)
end
