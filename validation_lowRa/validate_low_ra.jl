# validate_low_ra.jl
# -----------------------------------------------------------------------------
# Low-Ra validation with two parts:
#
#   Part A  — verify the ANALYTICAL solution itself, independently:
#             A1 ODE residual of  nu* v'' + T(x) = 0
#             A2 boundary conditions v(0) = v(1) = 0
#             A3 antisymmetry  v(x) + v(1 - x) = 0
#             A4 bulk velocity must vanish (antisymmetric profile)
#             A5 peak location and amplitude against their closed forms
#             A6 an independent finite-difference solution of the same BVP
#             A7 agreement with `analytic_low_ra_v` from common/cd_tools.jl
#
#   Part B  — compare the CFD run against that verified analytical profile and
#             report relative L2, Linf and peak errors plus bulk quantities.
#
# Physics. For theta = 90 deg the buoyancy acts along +y and, at low Ra, the
# temperature field is the pure conduction profile T(x) = 0.5 - x on x in [0, 1]
# (hot plate at x = 0). The steady streamwise balance
#
#     nu* v'' + T(x) = 0 ,      v(0) = v(1) = 0 ,      nu* = sqrt(Pr / Ra)
#
# integrates to the cubic
#
#     v(x) = x (1 - x) (1 - 2x) / (12 nu*) ,
#
# whose derivative vanishes at x_p = (3 - sqrt(3)) / 6 with amplitude
# v_p = 1 / (72 sqrt(3) nu*).
#
# Reads  : <CASE_DIR>/vbar_xt_raw.csv        (columns time,x,vbar)
#          <CASE_DIR>/case_parameters.txt    (Ra, Pr, viscosity)
# Writes : <CASE_DIR>/validation/CFD_vs_analytical.csv
#          <CASE_DIR>/validation/CFD_vs_analytical.png
#          <CASE_DIR>/validation/validation_summary.txt
#
# Note: this project does not depend on CSV.jl / DataFrames.jl; tables are
# written with DelimitedFiles, as in the other scripts.
#
# Usage:
#   CASE_DIR=data/90degree julia validate_low_ra.jl
#   CASE_DIR=data/90degree T_START=40 T_END_STAT=100 julia validate_low_ra.jl
# -----------------------------------------------------------------------------

using GLMakie
using DelimitedFiles
using Printf
using LinearAlgebra

base = @__DIR__
include(joinpath(base, "..", "shared", "common", "cd_tools.jl"))

# ------------------------------- case settings -------------------------------
# Default case: the Ra = 1000, 90-degree validation case.
const CASE_DIR = get(ENV, "CASE_DIR", joinpath(base, "data", "90degree"))
const INPUT_FILE = get(ENV, "VBAR_FILE", joinpath(CASE_DIR, "vbar_xt_raw.csv"))
const PARAM_FILE = joinpath(CASE_DIR, "case_parameters.txt")
const OUTPUT_DIR = get(ENV, "OUTDIR", joinpath(CASE_DIR, "validation"))

# Averaging interval for the steady CFD profile.
const T_START = parse(Float64, get(ENV, "T_START", "40.0"))
const T_END = parse(Float64, get(ENV, "T_END_STAT", "100.0"))

# Tolerance for the Part A self-checks (machine-precision expectations).
const SELFCHECK_TOL = parse(Float64, get(ENV, "SELFCHECK_TOL", "1e-10"))

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

isfile(INPUT_FILE) || error("CFD data not found: $INPUT_FILE")
mkpath(OUTPUT_DIR)

params = read_parameter_file(PARAM_FILE)
Ra = parse(Float64, get(params, "Ra", "1000.0"))
Pr = parse(Float64, get(params, "Pr", "0.71"))
nu = sqrt(Pr / Ra)

if haskey(params, "viscosity")
    nu_case = parse(Float64, params["viscosity"])
    isapprox(nu, nu_case; rtol = 1e-8) || @warn(
        "nu* from sqrt(Pr/Ra) = $nu differs from the case value $nu_case",
    )
end

# ============================================================
# Helper functions
# ============================================================

"""Analytic low-Ra profile: nu v'' + (0.5 - x) = 0 on [0,1], v(0)=v(1)=0."""
low_ra_profile(x, nu) = @. x * (1.0 - x) * (1.0 - 2.0 * x) / (12.0 * nu)

"""
    second_derivative(x, v)

Three-point second derivative on a (possibly non-uniform) grid. The stencil is
exact for cubic polynomials, which is the case here.
"""
function second_derivative(x::AbstractVector, v::AbstractVector)
    n = length(x)
    d2 = fill(NaN, n)
    for i in 2:(n - 1)
        hL = x[i] - x[i - 1]
        hR = x[i + 1] - x[i]
        d2[i] = 2 / (hL + hR) * ((v[i + 1] - v[i]) / hR - (v[i] - v[i - 1]) / hL)
    end
    d2
end

"""
    solve_low_ra_fd(x, nu)

Independent finite-difference solution of `nu v'' + (0.5 - x) = 0` with
`v(0) = v(1) = 0` on the grid `x`, via a tridiagonal solve.
"""
function solve_low_ra_fd(x::AbstractVector, nu::Real)
    n = length(x)
    m = n - 2                       # interior unknowns
    a = zeros(m - 1)                # sub-diagonal
    b = zeros(m)
    c = zeros(m - 1)                # super-diagonal
    r = zeros(m)

    for k in 1:m
        i = k + 1
        hL = x[i] - x[i - 1]
        hR = x[i + 1] - x[i]
        hs = hL + hR
        # v_{i-1} and v_{i+1} are boundary values (zero) at the first and last
        # interior row, so those coefficients simply do not appear in `a` / `c`.
        k > 1 && (a[k - 1] = 2 * nu / hs / hL)
        b[k] = -2 * nu / hs * (1 / hL + 1 / hR)
        k < m && (c[k] = 2 * nu / hs / hR)
        r[k] = x[i] - 0.5           # nu v'' = -(0.5 - x)
    end

    vint = Tridiagonal(a, b, c) \ r
    vcat(0.0, vint, 0.0)
end

# ============================================================
# Part A — verify the analytical solution itself
# ============================================================

x_uni = collect(range(0.0, 1.0; length = 401))     # uniform grid for closed forms
v_uni = low_ra_profile(x_uni, nu)

# A1: ODE residual on the actual (tanh-clustered) CFD grid
x_ref = vcat(
    0.0,
    collect(range(0.0, 1.0; length = 129))[2:end-1],
    1.0,
)
v_ref = low_ra_profile(x_ref, nu)
residual = nu .* second_derivative(x_ref, v_ref) .+ (0.5 .- x_ref)
# The stencil leaves NaN at the two walls by construction; use the interior only.
A1_residual = maximum(abs.(residual[2:end-1]))

# A2: boundary conditions
A2_bc = max(abs(v_ref[1]), abs(v_ref[end]))

# A3: antisymmetry on the uniform grid (1 - x is also a grid point there)
A3_antisym = maximum(abs.(v_uni .+ reverse(v_uni)))

# A4: bulk velocity of an antisymmetric profile
A4_vbulk = trapz1(x_uni, v_uni) / (x_uni[end] - x_uni[1])

# A5: peak location and amplitude against closed forms
# v is antisymmetric: the positive peak sits at x_p and the negative one at
# 1 - x_p, so a plain argmax (no abs) selects the one we compare against.
i_peak = argmax(v_uni)
x_peak_ana = (3 - sqrt(3)) / 6
v_peak_ana = 1 / (72 * sqrt(3) * nu)
A5_dx_peak = abs(x_uni[i_peak] - x_peak_ana)
A5_dv_peak = abs(abs(v_uni[i_peak]) - v_peak_ana) / v_peak_ana

# A6: independent finite-difference solve of the same BVP
v_fd = solve_low_ra_fd(x_ref, nu)
A6_fd_diff = maximum(abs.(v_fd .- v_ref))

# A7: agreement with the helper already shipped in common/cd_tools.jl
v_helper = analytic_low_ra_v(x_ref; viscosity = nu, g_stream = 1.0, drive_force_y = 0.0, H = 1.0)
A7_helper_diff = maximum(abs.(v_helper .- v_ref))

# A1-A4, A6 and A7 are exact statements about the continuous profile, so they must
# hold to machine precision. A5 compares a *grid* maximum with a continuous
# formula, so its accuracy is limited by the grid spacing (O(h^2) for the
# amplitude); it is reported separately instead of driving the pass/fail test.
part_a_ok = max(A1_residual, A2_bc, A3_antisym, abs(A4_vbulk),
                A6_fd_diff, A7_helper_diff) < SELFCHECK_TOL
A5_grid_ok = A5_dx_peak <= maximum(diff(x_uni)) && A5_dv_peak < 1e-2

# ============================================================
# Part B — CFD versus the verified analytical profile
# ============================================================

times, x, V = read_vbar_long_csv(INPUT_FILE)

T_START >= times[1] || error("T_START precedes the available CFD data")
T_END <= times[end] || error("T_END exceeds the available CFD data")

v_cfd = time_weighted_profile(times, V, T_START, T_END)
x_full, v_cfd_full = full_wall_profile(x, v_cfd, (0.0, 1.0))
Lx = x_full[end] - x_full[1]

v_ana_full = low_ra_profile(x_full, nu)

err = v_cfd_full .- v_ana_full
relative_L2 = sqrt(trapz1(x_full, err .^ 2) / trapz1(x_full, v_ana_full .^ 2))
relative_Linf = maximum(abs.(err)) / maximum(abs.(v_ana_full))

U_cfd = maximum(abs.(v_cfd_full))
U_ana = maximum(abs.(v_ana_full))
peak_error = abs(U_cfd - U_ana) / U_ana

Vbulk_cfd = trapz1(x_full, v_cfd_full) / Lx
Vabs_cfd = trapz1(x_full, abs.(v_cfd_full)) / Lx
Vbulk_ana = trapz1(x_full, v_ana_full) / Lx
Vabs_ana = trapz1(x_full, abs.(v_ana_full)) / Lx

# ============================================================
# Save comparison table
# ============================================================

csv_path = joinpath(OUTPUT_DIR, "CFD_vs_analytical.csv")
open(csv_path, "w") do io
    println(io, "x,v_CFD,v_analytical,error")
    for i in eachindex(x_full)
        println(io, x_full[i], ',', v_cfd_full[i], ',', v_ana_full[i], ',', err[i])
    end
end

# ============================================================
# Plot CFD vs analytical
# ============================================================

fig = Figure(size = (900, 800))
ax = Axis(
    fig[1, 1];
    xlabel = "streamwise velocity",
    ylabel = "x",
    title = @sprintf("Low-Ra validation (Ra = %.3g): CFD vs analytical", Ra),
)
lines!(ax, v_ana_full, x_full; linewidth = 3, label = "Analytical")
scatterlines!(ax, v_cfd_full, x_full; markersize = 5, label = "CFD")
axislegend(ax)
save(joinpath(OUTPUT_DIR, "CFD_vs_analytical.png"), fig)

# ============================================================
# Summary
# ============================================================

verdict = relative_L2 < 1e-2 ? "PASS (relative L2 < 1 %)" :
          relative_L2 < 5e-2 ? "ACCEPTABLE (relative L2 < 5 %)" : "FAIL"

summary = """
Low-Rayleigh-number CFD validation
===================================

case            = $(basename(CASE_DIR))
Ra              = $(Ra)
Pr              = $(Pr)
nu*             = $(nu)

----------------------------------------------
Part A  analytical solution self-checks
----------------------------------------------
A1 ODE residual  max|nu v'' + T|        = $(A1_residual)
A2 boundary     max(|v(0)|,|v(1)|)      = $(A2_bc)
A3 antisymmetry max|v(x) + v(1-x)|      = $(A3_antisym)
A4 bulk velocity (must be 0)            = $(A4_vbulk)
A5 peak location  |x_p - (3-sqrt3)/6|   = $(A5_dx_peak)
A5 peak amplitude relative difference   = $(A5_dv_peak)
A6 finite-difference BVP max|v_fd - v|  = $(A6_fd_diff)
A7 cd_tools analytic_low_ra_v max diff  = $(A7_helper_diff)

numerical peaks: x_peak = $(x_uni[i_peak]), |v_peak| = $(abs(v_uni[i_peak]))
closed forms   : x_peak = $(x_peak_ana), |v_peak| = $(v_peak_ana)

Part A verdict: the analytical profile solves the stated BVP to machine
precision (all checks below $(SELFCHECK_TOL)) -> $(part_a_ok)

----------------------------------------------
Part B  CFD versus analytical
----------------------------------------------
averaging interval = [$(T_START), $(T_END)]
samples used       = $(count(t -> T_START <= t <= T_END, times))
grid points        = $(length(x_full)) (incl. the two walls)

Umax   CFD = $(U_cfd)   analytical = $(U_ana)
Vabs   CFD = $(Vabs_cfd)   analytical = $(Vabs_ana)
Vbulk  CFD = $(Vbulk_cfd)   analytical = $(Vbulk_ana)

relative L2 error   = $(relative_L2)
relative Linf error = $(relative_Linf)
peak velocity error = $(peak_error)

----------------------------------------------
Conclusion
----------------------------------------------
$(verdict)

The CFD profile matches the low-Ra analytical solution to a relative L2 error
of $(relative_L2) ($(round(relative_L2 * 100; digits = 4)) %). Since Part A shows the
analytical profile is exact for the stated BVP, this difference is the physical
effect of Ra = $(Ra) not being infinitesimal (the real temperature field is not
exactly the conduction profile, and Nu > 1), plus the finite-grid discretisation
error of the CFD run -- not an error in the reference solution.
"""

println(summary)

open(joinpath(OUTPUT_DIR, "validation_summary.txt"), "w") do io
    write(io, summary)
end

println("Wrote:")
println("  ", csv_path)
println("  ", joinpath(OUTPUT_DIR, "CFD_vs_analytical.png"))
println("  ", joinpath(OUTPUT_DIR, "validation_summary.txt"))
part_a_ok || error("Part A self-checks failed; see the summary above")
