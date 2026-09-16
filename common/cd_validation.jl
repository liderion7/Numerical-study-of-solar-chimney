# common/cd_validation.jl
# -----------------------------------------------------------------------------
# Huang et al. (2022) validation utilities for the existing Task C/D solver.
#
# This file is intentionally independent of Task C / Task D post-processing.
# It evaluates the full 2-D velocity and temperature fields needed for the
# published RBC benchmark quantities:
#
#   Re   = H/nu * sqrt(<u^2 + v^2>_{V,t})
#   Re_x = H/nu * sqrt(<u^2>_{V,t})
#   Re_y = H/nu * sqrt(<v^2>_{V,t})
#   Nu   = 1 + sqrt(Ra*Pr) <v*T>_{V,t}
#
# Coordinate mapping for THIS project when theta = 0 deg:
#
#   code x  = Huang vertical y  (hot plate -> cold plate)
#   code y  = Huang horizontal x (periodic direction)
#   code u_x = Huang vertical velocity v
#   code u_y = Huang horizontal velocity u
#
# Therefore:
#   Huang Re_x uses <code u_y^2>
#   Huang Re_y uses <code u_x^2>
#   Nu uses <code u_x * T>
# -----------------------------------------------------------------------------

using IncompressibleNavierStokes
using Statistics
using DelimitedFiles
using Observables
using Printf

# Published Huang et al. (2022), Table 1, NS plates + periodic sidewalls (PD).
# These are the cases directly compatible with the current project geometry.
const HUANG_NS_PD_REFERENCE = Dict(
    1.0e6 => (Nu = 7.9,  Re = 67.0,   Rex = 51.3,   Rey = 43.1,   tavg = 1000.0),
    1.5e6 => (Nu = 8.8,  Re = 84.6,   Rex = 65.2,   Rey = 53.9,   tavg = 1000.0),
    2.0e6 => (Nu = 9.6,  Re = 99.7,   Rex = 77.3,   Rey = 63.0,   tavg = 1000.0),
    3.0e6 => (Nu = 10.4, Re = 127.6,  Rex = 99.7,   Rey = 79.6,   tavg = 1000.0),
    5.0e6 => (Nu = 11.7, Re = 171.6,  Rex = 134.5,  Rey = 106.6,  tavg = 1000.0),
    1.0e7 => (Nu = 14.3, Re = 257.6,  Rex = 202.4,  Rey = 159.3,  tavg = 800.0),
    2.0e7 => (Nu = 17.4, Re = 384.4,  Rex = 302.9,  Rey = 236.7,  tavg = 800.0),
    5.0e7 => (Nu = 22.4, Re = 634.5,  Rex = 501.3,  Rey = 389.0,  tavg = 800.0),
    1.0e8 => (Nu = 26.6, Re = 958.9,  Rex = 759.1,  Rey = 585.9,  tavg = 1000.0),
    2.0e8 => (Nu = 32.1, Re = 1446.2, Rex = 1118.1, Rey = 917.2,  tavg = 800.0),
    5.0e8 => (Nu = 41.3, Re = 2619.3, Rex = 2026.7, Rey = 1659.3, tavg = 800.0),
    1.0e9 => (Nu = 50.9, Re = 3575.0, Rex = 2820.5, Rey = 2196.7, tavg = 500.0),
)

"""Return the Huang NS/PD reference row for `Ra`, allowing tiny float roundoff."""
function huang_reference_ns_pd(Ra::Real)
    for (key, val) in HUANG_NS_PD_REFERENCE
        if isapprox(Float64(Ra), key; rtol = 1e-12, atol = 0.0)
            return val
        end
    end
    return nothing
end

"""
    build_huang_initial_state(setup, psolver; ...)

Initial condition for Huang validation.

Modes:
- `"zero"`      : u=0, T=0 in the interior, closest to the paper statement.
- `"conduction"`: u=0, T=0.5-x/H.
- `"perturbed"` : conduction profile plus a tiny 2-D perturbation.

`perturbed` is the recommended practical default because an exactly symmetric
RBC state can otherwise remain numerically symmetric for a long time.
`perturbation_waves` sets how many horizontal perturbation wavelengths fit in
the periodic length, which selects the convection-roll pattern (1 = original).
"""
function build_huang_initial_state(
    setup,
    psolver;
    gap::Real = 1.0,
    periodic_length::Real = 2.0,
    mode::AbstractString = "perturbed",
    epsilon::Real = 1e-6,
    perturbation_waves::Integer = 1,
)
    gap > 0 || throw(ArgumentError("gap must be positive"))
    periodic_length > 0 || throw(ArgumentError("periodic_length must be positive"))
    epsilon >= 0 || throw(ArgumentError("epsilon must be non-negative"))
    perturbation_waves >= 1 || throw(ArgumentError("perturbation_waves must be >= 1"))

    mode_l = lowercase(strip(mode))

    tempfun = if mode_l == "zero" || mode_l == "huang"
        # Stationary isothermal state, u = 0 and theta = 0. `"huang"` is the
        # name used for the initial condition stated in the benchmark text; it
        # builds exactly the same field as `"zero"`.
        (x, y) -> zero(x)
    elseif mode_l == "conduction"
        (x, y) -> 0.5 - x / gap
    elseif mode_l == "perturbed"
        (x, y) -> begin
            base = 0.5 - x / gap
            # Zero perturbation at both plates; `perturbation_waves` horizontal
            # wavelengths across the periodic length (1 = the original setup).
            base +
            epsilon * sinpi(x / gap) *
            sin(2π * perturbation_waves * y / periodic_length)
        end
    else
        throw(ArgumentError("unknown Huang INIT_MODE='$mode'; use huang/zero, conduction, or perturbed"))
    end

    (
        u = velocityfield(
            setup,
            (dim, x, y) -> zero(x);
            psolver,
        ),
        temp = temperaturefield(setup, tempfun),
    )
end

"""Extract scalar temperature on pressure/cell-centre points."""
function _huang_temperature_at_p(temp, setup)
    (; Ip) = setup
    ix, iy = Ip.indices
    Traw = Array(temp)

    if ndims(Traw) == 2
        return Float64.(Traw[ix, iy])
    elseif ndims(Traw) == 3 && size(Traw, 3) == 1
        return Float64.(Traw[ix, iy, 1])
    else
        error("Unexpected temperature array size $(size(Traw)); expected a 2-D scalar field")
    end
end

"""
    huang_instantaneous_stats(u, temp, setup; Ra, Pr, H=1.0)

Return full-field instantaneous volume statistics in Huang's coordinate names.
The volume average uses the actual non-uniform finite-volume cell widths.
"""
function huang_instantaneous_stats(
    u,
    temp,
    setup;
    Ra::Real,
    Pr::Real,
    H::Real = 1.0,
)
    Ra > 0 || throw(ArgumentError("Ra must be positive"))
    Pr > 0 || throw(ArgumentError("Pr must be positive"))
    H > 0 || throw(ArgumentError("H must be positive"))

    (; Δ, Ip) = setup
    ix, iy = Ip.indices

    # Staggered velocity -> pressure / cell-centre points.
    up = Array(IncompressibleNavierStokes.interpolate_u_p(u, setup))
    ux_code = Float64.(up[ix, iy, 1])  # Huang vertical velocity v
    uy_code = Float64.(up[ix, iy, 2])  # Huang horizontal velocity u
    T = _huang_temperature_at_p(temp, setup)

    size(ux_code) == size(uy_code) == size(T) ||
        error("Velocity/temperature pressure-grid sizes do not match: " *
              "ux=$(size(ux_code)), uy=$(size(uy_code)), T=$(size(T))")

    dx = Float64.(Array(Δ[1])[ix])
    dy = Float64.(Array(Δ[2])[iy])
    W = reshape(dx, :, 1) .* reshape(dy, 1, :)
    volume = sum(W)
    volume > 0 || error("non-positive computational volume")

    mean_vertical2 = sum((ux_code .^ 2) .* W) / volume
    mean_horizontal2 = sum((uy_code .^ 2) .* W) / volume
    mean_vT = sum((ux_code .* T) .* W) / volume

    nu_star = sqrt(Pr / Ra)

    # Paper coordinate convention after mapping described at top of this file.
    Rex_inst = H / nu_star * sqrt(max(mean_horizontal2, 0.0))
    Rey_inst = H / nu_star * sqrt(max(mean_vertical2, 0.0))
    Re_inst = H / nu_star * sqrt(max(mean_horizontal2 + mean_vertical2, 0.0))

    # For the benchmark H=1, ΔT=1 and free-fall nondimensionalization.
    Nu_inst = 1 + sqrt(Ra * Pr) * mean_vT

    (
        ;
        mean_horizontal2,
        mean_vertical2,
        mean_vT,
        Re_inst,
        Rex_inst,
        Rey_inst,
        Nu_inst,
    )
end

"""
    make_huang_statistics_recorder(setup; ...)

Sample full-field statistics at physical-time intervals and stream them to CSV.
Only samples with t >= `start_time` are recorded, so the transient can be
excluded without storing unnecessary full-field statistics.
"""
function make_huang_statistics_recorder(
    setup;
    Ra::Real,
    Pr::Real,
    H::Real = 1.0,
    sample_interval::Real = 0.5,
    start_time::Real = 200.0,
    output_path::Union{Nothing, AbstractString} = nothing,
)
    sample_interval > 0 || throw(ArgumentError("sample_interval must be positive"))

    times = Float64[]
    mean_horizontal2 = Float64[]
    mean_vertical2 = Float64[]
    mean_vT = Float64[]
    next_sample = Ref(Float64(start_time))
    io_ref = Ref{Any}(nothing)

    if output_path !== nothing
        mkpath(dirname(output_path))
        io_ref[] = open(output_path, "w")
        println(io_ref[], "time,mean_horizontal2,mean_vertical2,mean_vT,Re_inst,Rex_inst,Rey_inst,Nu_inst")
        flush(io_ref[])
    end

    close_io! = () -> begin
        io_ref[] !== nothing && close(io_ref[])
        io_ref[] = nothing
    end

    initialize = function (state_obs)
        on(state_obs) do (; u, temp, t)
            tf = Float64(t)
            tf + 100eps(tf + 1) < next_sample[] && return

            s = huang_instantaneous_stats(u, temp, setup; Ra, Pr, H)

            push!(times, tf)
            push!(mean_horizontal2, s.mean_horizontal2)
            push!(mean_vertical2, s.mean_vertical2)
            push!(mean_vT, s.mean_vT)

            if io_ref[] !== nothing
                println(
                    io_ref[],
                    repr(tf), ',',
                    repr(s.mean_horizontal2), ',',
                    repr(s.mean_vertical2), ',',
                    repr(s.mean_vT), ',',
                    repr(s.Re_inst), ',',
                    repr(s.Rex_inst), ',',
                    repr(s.Rey_inst), ',',
                    repr(s.Nu_inst),
                )
                flush(io_ref[])
            end

            while next_sample[] <= tf + 100eps(tf + 1)
                next_sample[] += sample_interval
            end
        end
        nothing
    end

    (;
        initialize,
        finalize = (initialized, state) -> (close_io!(); initialized),
        times,
        mean_horizontal2,
        mean_vertical2,
        mean_vT,
        close_io!,
    )
end

"""Lightweight progress monitor for long Huang benchmark runs."""
function make_huang_runtime_monitor(
    setup;
    Ra::Real,
    Pr::Real,
    H::Real = 1.0,
    t_end::Real,
    report_interval::Real = 20.0,
)
    report_interval > 0 || throw(ArgumentError("report_interval must be positive"))

    next_report = Ref(0.0)
    wall_start = time()

    initialize = function (state_obs)
        on(state_obs) do (; u, temp, t, n)
            tf = Float64(t)
            tf + 100eps(tf + 1) < next_report[] && return

            s = huang_instantaneous_stats(u, temp, setup; Ra, Pr, H)
            elapsed = time() - wall_start
            eta = tf > 0 ? max(0.0, (t_end - tf) * elapsed / tf) : NaN

            @printf(
                "[Huang] t=%8.3f / %.1f  step=%d  Re_inst=%8.3f  Rex=%8.3f  Rey=%8.3f  Nu_inst=%8.3f  wall=%.1fs  ETA=%.1fs\n",
                tf,
                t_end,
                Int(n),
                s.Re_inst,
                s.Rex_inst,
                s.Rey_inst,
                s.Nu_inst,
                elapsed,
                eta,
            )

            while next_report[] <= tf + 100eps(tf + 1)
                next_report[] += report_interval
            end
        end
        nothing
    end

    (; initialize, finalize = (initialized, state) -> initialized)
end

"""Trapezoidal time average for non-uniform sample times."""
function _huang_time_average(times::AbstractVector, values::AbstractVector)
    length(times) == length(values) || throw(DimensionMismatch("times/values length mismatch"))
    length(times) >= 2 || error("Need at least two validation samples for time averaging")
    all(diff(times) .> 0) || error("Validation sample times must be strictly increasing")

    duration = times[end] - times[1]
    duration > 0 || error("Validation averaging duration must be positive")

    integral = sum(
        (times[2:end] .- times[1:end-1]) .* (values[2:end] .+ values[1:end-1]) ./ 2,
    )
    integral / duration
end

"""
    finalize_huang_statistics(rec; Ra, Pr, H=1.0)

Compute the authoritative time+volume averaged Nu, Re, Rex and Rey. Note that
Re is formed AFTER time-averaging the squared velocities, matching the paper.
"""
function finalize_huang_statistics(rec; Ra::Real, Pr::Real, H::Real = 1.0)
    length(rec.times) >= 2 || error(
        "Too few Huang validation samples ($(length(rec.times))). " *
        "Increase T_END or reduce AVG_START/SAMPLE_DT.",
    )

    mh2 = _huang_time_average(rec.times, rec.mean_horizontal2)
    mv2 = _huang_time_average(rec.times, rec.mean_vertical2)
    mvT = _huang_time_average(rec.times, rec.mean_vT)

    nu_star = sqrt(Pr / Ra)
    Rex = H / nu_star * sqrt(max(mh2, 0.0))
    Rey = H / nu_star * sqrt(max(mv2, 0.0))
    Re = H / nu_star * sqrt(max(mh2 + mv2, 0.0))
    Nu = 1 + sqrt(Ra * Pr) * mvT

    (
        Nu,
        Re,
        Rex,
        Rey,
        Re_identity = sqrt(Rex^2 + Rey^2),
        mean_horizontal2 = mh2,
        mean_vertical2 = mv2,
        mean_vT = mvT,
        t_start = rec.times[1],
        t_end = rec.times[end],
        tavg = rec.times[end] - rec.times[1],
        nsamples = length(rec.times),
    )
end

_huang_relerr(value, ref) = 100 * (value - ref) / ref
_huang_abserr_pct(value, ref) = abs(_huang_relerr(value, ref))

"""Write a human-readable benchmark summary and a one-row CSV."""
function write_huang_summary(
    outdir,
    stats;
    Ra::Real,
    Pr::Real,
    N_wall::Integer,
    Ny::Integer,
    init_mode::AbstractString,
    init_epsilon::Real,
)
    mkpath(outdir)
    ref = huang_reference_ns_pd(Ra)

    txt = joinpath(outdir, "huang_validation_summary.txt")
    open(txt, "w") do io
        println(io, "Huang et al. (2022) RBC validation")
        println(io, "====================================")
        println(io, "Benchmark subset = NS plates / periodic sidewalls (NS/PD)")
        println(io, "Project mapping  = code x -> paper vertical y; code y -> paper horizontal x")
        println(io, "theta            = 0 deg")
        println(io, "Ra               = ", Ra)
        println(io, "Pr               = ", Pr)
        println(io, "aspect ratio     = 2")
        println(io, "N_wall           = ", N_wall)
        println(io, "Ny               = ", Ny)
        println(io, "initial mode     = ", init_mode)
        println(io, "initial epsilon  = ", init_epsilon)
        println(io, "sample count     = ", stats.nsamples)
        println(io, "average t_start  = ", stats.t_start)
        println(io, "average t_end    = ", stats.t_end)
        println(io, "average duration = ", stats.tavg)
        println(io)
        @printf(io, "CFD Nu  = %.8f\n", stats.Nu)
        @printf(io, "CFD Re  = %.8f\n", stats.Re)
        @printf(io, "CFD Rex = %.8f\n", stats.Rex)
        @printf(io, "CFD Rey = %.8f\n", stats.Rey)
        @printf(io, "sqrt(Rex^2+Rey^2) = %.8f\n", stats.Re_identity)
        @printf(io, "Re identity difference = %.3e\n", stats.Re - stats.Re_identity)

        if ref === nothing
            println(io)
            println(io, "No exact Huang NS/PD reference row is stored for this Ra.")
        else
            println(io)
            println(io, "Quantity        CFD             Huang           signed err %     abs err %")
            println(io, "----------------------------------------------------------------------------")
            @printf(io, "Nu        %14.6f  %14.6f  %14.4f  %12.4f\n",
                    stats.Nu, ref.Nu, _huang_relerr(stats.Nu, ref.Nu), _huang_abserr_pct(stats.Nu, ref.Nu))
            @printf(io, "Re        %14.6f  %14.6f  %14.4f  %12.4f\n",
                    stats.Re, ref.Re, _huang_relerr(stats.Re, ref.Re), _huang_abserr_pct(stats.Re, ref.Re))
            @printf(io, "Rex       %14.6f  %14.6f  %14.4f  %12.4f\n",
                    stats.Rex, ref.Rex, _huang_relerr(stats.Rex, ref.Rex), _huang_abserr_pct(stats.Rex, ref.Rex))
            @printf(io, "Rey       %14.6f  %14.6f  %14.4f  %12.4f\n",
                    stats.Rey, ref.Rey, _huang_relerr(stats.Rey, ref.Rey), _huang_abserr_pct(stats.Rey, ref.Rey))
            println(io)
            println(io, "Huang reported averaging time for this row = ", ref.tavg)
        end
    end

    csv = joinpath(outdir, "huang_validation_comparison.csv")
    open(csv, "w") do io
        println(io, "Ra,Pr,N_wall,Ny,tavg,Nu_cfd,Nu_ref,Nu_abs_err_pct,Re_cfd,Re_ref,Re_abs_err_pct,Rex_cfd,Rex_ref,Rex_abs_err_pct,Rey_cfd,Rey_ref,Rey_abs_err_pct")
        if ref === nothing
            println(io,
                Ra, ',', Pr, ',', N_wall, ',', Ny, ',', stats.tavg, ',',
                stats.Nu, ",NaN,NaN,",
                stats.Re, ",NaN,NaN,",
                stats.Rex, ",NaN,NaN,",
                stats.Rey, ",NaN,NaN",
            )
        else
            println(io,
                Ra, ',', Pr, ',', N_wall, ',', Ny, ',', stats.tavg, ',',
                stats.Nu, ',', ref.Nu, ',', _huang_abserr_pct(stats.Nu, ref.Nu), ',',
                stats.Re, ',', ref.Re, ',', _huang_abserr_pct(stats.Re, ref.Re), ',',
                stats.Rex, ',', ref.Rex, ',', _huang_abserr_pct(stats.Rex, ref.Rex), ',',
                stats.Rey, ',', ref.Rey, ',', _huang_abserr_pct(stats.Rey, ref.Rey),
            )
        end
    end

    (; txt, csv, ref)
end

"""
    write_huang_field_snapshot(path, u, temp, setup; t)

Write one full-field snapshot as a long-table CSV with columns
`ix,iy,x,y,T,ux,uy`, sampled on the pressure grid. `ux` is the wall-normal
(paper-vertical) velocity component and `uy` the periodic (paper-horizontal)
one, following the mapping documented at the top of this file. The first line
records the physical time as a comment so the snapshot can be identified.
"""
function write_huang_field_snapshot(path, u, temp, setup; t::Real)
    (; Ip) = setup
    ix, iy = Ip.indices

    up = Array(IncompressibleNavierStokes.interpolate_u_p(u, setup))
    ux = Float64.(up[ix, iy, 1])
    uy = Float64.(up[ix, iy, 2])
    T = _huang_temperature_at_p(temp, setup)

    size(ux) == size(uy) == size(T) ||
        error("Velocity/temperature pressure-grid sizes do not match")

    x = Float64.(Array(setup.xp[1])[ix])
    y = Float64.(Array(setup.xp[2])[iy])

    mkpath(dirname(path))
    open(path, "w") do io
        @printf(io, "# t = %.8f\n", Float64(t))
        println(io, "ix,iy,x,y,T,ux,uy")
        for (j, yv) in enumerate(y)
            for (i, xv) in enumerate(x)
                println(
                    io,
                    i, ',', j, ',',
                    xv, ',', yv, ',',
                    T[i, j], ',', ux[i, j], ',', uy[i, j],
                )
            end
        end
    end
    path
end

"""
    make_huang_field_recorder(setup; output_dir, snapshot_times, prefix="field")

Save full-field snapshots at the requested physical times so that the
convection-roll pattern can be inspected after the run. With an empty
`snapshot_times` this processor is a no-op, so it can always be handed to
`solve_unsteady` without changing the run.
"""
function make_huang_field_recorder(
    setup;
    output_dir::AbstractString,
    snapshot_times::AbstractVector{<:Real} = Float64[],
    prefix::AbstractString = "field",
)
    times = sort(Float64.(collect(snapshot_times)))
    next_idx = Ref(1)
    written = String[]

    isempty(times) || mkpath(output_dir)

    initialize = function (state_obs)
        on(state_obs) do (; u, temp, t)
            tf = Float64(t)
            while next_idx[] <= length(times) && tf + 100eps(tf + 1) >= times[next_idx[]]
                path = joinpath(
                    output_dir,
                    @sprintf("%s_t%.1f.csv", prefix, times[next_idx[]]),
                )
                write_huang_field_snapshot(path, u, temp, setup; t = tf)
                push!(written, path)
                next_idx[] += 1
            end
        end
        nothing
    end

    (;
        initialize,
        finalize = (initialized, state) -> initialized,
        written,
        snapshot_times = times,
    )
end
