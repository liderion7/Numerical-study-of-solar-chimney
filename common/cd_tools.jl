using DelimitedFiles
using Statistics
using LinearAlgebra
using Observables

# -----------------------------------------------------------------------------
# Core numerical integration helpers
# -----------------------------------------------------------------------------

"""Trapezoidal integral for a one-dimensional non-uniform grid."""
function trapz1(x::AbstractVector, f::AbstractVector)
    length(x) == length(f) || throw(DimensionMismatch("x and f must have the same length"))
    length(x) >= 2 || throw(ArgumentError("at least two points are required"))
    sum((x[2:end] .- x[1:end-1]) .* (f[2:end] .+ f[1:end-1]) ./ 2)
end

"""
    full_wall_profile(x, v, xlims)

Add the two no-slip wall values v=0 to a profile stored at cell centres.
"""
function full_wall_profile(x::AbstractVector, v::AbstractVector, xlims)
    xmin, xmax = Float64.(xlims)
    xfull = vcat(xmin, Float64.(x), xmax)
    vfull = vcat(0.0, Float64.(v), 0.0)
    xfull, vfull
end

# -----------------------------------------------------------------------------
# Task-C profile extraction from the CFD state
# -----------------------------------------------------------------------------

"""
    task_c_vbar_profile(u, setup)

Return `(x, vbar)` where

    vbar(x,t) = (1/Ly) * integral v(x,y,t) dy,

and `v` is the y-component of velocity. This is the exact project convention:

- y is vertical / streamwise,
- x is wall-normal,
- theta=90° is parallel to y,
- T=T(x), hence the streamwise velocity profile is ultimately a function of x.

The averaging uses the true y-cell widths and therefore remains correct if the
grid is changed later.
"""
function task_c_vbar_profile(u, setup)
    (; xp, Δ, Ip) = setup
    ix, iy = Ip.indices

    # Staggered velocity -> pressure / cell-centre points.
    up = Array(IncompressibleNavierStokes.interpolate_u_p(u, setup))

    # v = y-component, shape (Nx, Ny).
    v = up[ix, iy, 2]

    x = Array(xp[1])[ix]
    Δy = Array(Δ[2])[iy]

    Ly_weight = sum(Δy)
    vbar = vec(sum(v .* reshape(Δy, 1, :); dims = 2)) ./ Ly_weight

    Float64.(x), Float64.(vbar)
end

"""
    instantaneous_bulk_velocity(x, vbar, xlims)

Signed cross-section bulk velocity based on the Task-D definition.
"""
function instantaneous_bulk_velocity(x, vbar, xlims)
    xfull, vfull = full_wall_profile(x, vbar, xlims)
    Lx = xfull[end] - xfull[1]
    trapz1(xfull, vfull) / Lx
end

"""RMS amplitude of the one-dimensional profile; useful for Task-C diagnostics."""
function profile_rms(x, vbar, xlims)
    xfull, vfull = full_wall_profile(x, vbar, xlims)
    Lx = xfull[end] - xfull[1]
    sqrt(trapz1(xfull, vfull .^ 2) / Lx)
end

# -----------------------------------------------------------------------------
# Raw CFD recorder
# -----------------------------------------------------------------------------

"""
    make_task_c_recorder(setup; sample_interval=0.2, start_time=0.0)

Processor that stores the complete Task-C data set `vbar(x,t)`. Sampling is
controlled by physical time, not solver step count, so adaptive time stepping
does not change the intended output cadence.
"""
function make_task_c_recorder(
    setup;
    sample_interval::Real = 0.2,
    start_time::Real = 0.0,
    output_path::Union{Nothing, AbstractString} = nothing,
)
    sample_interval > 0 || throw(ArgumentError("sample_interval must be positive"))

    times = Float64[]
    profiles = Vector{Vector{Float64}}()
    xref = Ref{Union{Nothing, Vector{Float64}}}(nothing)
    next_sample = Ref(Float64(start_time))
    io_ref = Ref{Any}(nothing)

    if output_path !== nothing
        mkpath(dirname(output_path))
        io_ref[] = open(output_path, "w")
        println(io_ref[], "time,x,vbar")
        flush(io_ref[])
    end

    close_io! = () -> begin
        io_ref[] !== nothing && close(io_ref[])
        io_ref[] = nothing
    end

    initialize = function (state_obs)
        on(state_obs) do (; u, t)
            tf = Float64(t)
            tf + 100eps(tf + 1) < next_sample[] && return

            x, vbar = task_c_vbar_profile(u, setup)
            isnothing(xref[]) && (xref[] = copy(x))

            push!(times, tf)
            push!(profiles, copy(vbar))

            # Stream to disk continuously so a crash does not lose all raw data.
            if io_ref[] !== nothing
                for j in eachindex(x)
                    println(
                        io_ref[],
                        repr(tf), ',',
                        repr(Float64(x[j])), ',',
                        repr(Float64(vbar[j])),
                    )
                end
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
        profiles,
        xref,
        close_io!,
    )
end

"""
    make_runtime_monitor(setup; theta, t_end, report_interval=5.0,
                         output_path=nothing, xlims=(0.0, 1.0))

Processor that prints a one-line status every `report_interval` simulation time
units and continuously appends `runtime_monitor.csv`. The classification is a
lightweight online screening only; the authoritative decision is Task C.
"""
function make_runtime_monitor(
    setup;
    theta::Real = NaN,
    t_end::Real = NaN,
    report_interval::Real = 5.0,
    output_path::Union{Nothing, AbstractString} = nothing,
    xlims = (0.0, 1.0),
)
    report_interval > 0 || throw(ArgumentError("report_interval must be positive"))

    next_report = Ref(Float64(0.0))
    last_t = Ref(Float64(0.0))
    last_n = Ref(0)
    start_wall = Ref(time())
    rms_times = Float64[]
    rms_values = Float64[]
    io_ref = Ref{Any}(nothing)

    if output_path !== nothing
        mkpath(dirname(output_path))
        io_ref[] = open(output_path, "w")
        println(io_ref[], "wall_time,theta,simulation_time,step,dt,umax,profile_rms,Vbulk_inst")
        flush(io_ref[])
    end

    close_io! = () -> begin
        io_ref[] !== nothing && close(io_ref[])
        io_ref[] = nothing
    end

    function light_status()
        t_now = rms_times[end]
        t_now < 40.0 && return "transient"
        peaks = detect_cycle_peaks(rms_times, rms_values; min_separation = 20.0)
        if length(peaks) >= 3
            pt = rms_times[peaks]
            periods = diff(pt)
            pmean = mean(periods)
            cv = std(periods) / (abs(pmean) + eps(Float64))
            cv <= 0.10 && return "periodic-candidate"
        end
        if length(rms_values) >= 4
            recent = rms_values[(end - 2):end]
            rel = abs(recent[end] - recent[1]) / (abs(recent[end]) + eps(Float64))
            rel < 0.05 && return "stationary-candidate"
        end
        "transient"
    end

    initialize = function (state_obs)
        on(state_obs) do (; u, temp, t, n)
            tf = Float64(t)
            tf + 100eps(tf + 1) < next_report[] && return

            x, vbar = task_c_vbar_profile(u, setup)
            rms = profile_rms(x, vbar, xlims)
            vbulk = instantaneous_bulk_velocity(x, vbar, xlims)
            umax = maximum(abs, Array(u))

            n_int = Int(n)
            wall_now = time()
            elapsed = wall_now - start_wall[]
            dt_mean = last_n[] == 0 ? NaN : (tf - last_t[]) / (n_int - last_n[])
            eta = isfinite(t_end) && tf > 0 ? (t_end - tf) * elapsed / tf : NaN

            push!(rms_times, tf)
            push!(rms_values, rms)
            status = light_status()
            peak_count = length(detect_cycle_peaks(rms_times, rms_values; min_separation = 20.0))

            println(
                "[theta=", theta, " deg] t = ", round(tf; digits = 3), " / ", t_end,
                " step = ", n_int,
                " dt = ", round(dt_mean; sigdigits = 3),
                " wall = ", round(elapsed; digits = 1), " s",
                " ETA = ", round(eta; digits = 1), " s",
                " umax = ", round(umax; sigdigits = 3),
                " profile_RMS = ", round(rms; sigdigits = 3),
                " Vbulk_inst = ", round(vbulk; sigdigits = 3),
                " peak_count = ", peak_count,
                " Task-C = ", status,
            )

            if io_ref[] !== nothing
                println(io_ref[], wall_now, ',', theta, ',', tf, ',', n_int, ',', dt_mean, ',', umax, ',', rms, ',', vbulk)
                flush(io_ref[])
            end

            last_t[] = tf
            last_n[] = n_int
            while next_report[] <= tf + 100eps(tf + 1)
                next_report[] += report_interval
            end
        end
        nothing
    end

    (;
        initialize,
        finalize = (initialized, state) -> (close_io!(); initialized),
        close_io!,
    )
end

"""Write Task-C raw data in clean long format: time,x,vbar."""
function write_vbar_long_csv(path, times, x, profiles)
    length(times) == length(profiles) || throw(DimensionMismatch("times/profiles length mismatch"))
    isempty(times) && error("no Task-C samples were recorded")

    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "time,x,vbar")
        for (k, t) in pairs(times)
            p = profiles[k]
            length(p) == length(x) || throw(DimensionMismatch("profile length mismatch"))
            for j in eachindex(x)
                println(io, repr(Float64(t)), ',', repr(Float64(x[j])), ',', repr(Float64(p[j])))
            end
        end
    end
    path
end

"""
    read_vbar_long_csv(path)

Read `time,x,vbar` and return a rectangular data set `(times, x, V)` where
`V[i,j] = vbar(x[j], times[i])`.
"""
function read_vbar_long_csv(path)
    raw = readdlm(path, ',', Float64, '\n'; skipstart = 1)
    size(raw, 2) == 3 || error("expected exactly three columns: time,x,vbar")

    tcol = raw[:, 1]
    xcol = raw[:, 2]
    vcol = raw[:, 3]

    times = unique(tcol)
    x = unique(xcol)

    nt = length(times)
    nx = length(x)
    nt * nx == length(vcol) || error("data are incomplete or contain inconsistent time/x blocks")

    V = Matrix{Float64}(undef, nt, nx)
    for i in 1:nt
        rows = ((i - 1) * nx + 1):(i * nx)
        all(tcol[rows] .== times[i]) || error("time blocks are not contiguous")
        all(xcol[rows] .== x) || error("x grid changes between time samples")
        V[i, :] .= vcol[rows]
    end

    Float64.(times), Float64.(x), V
end

# -----------------------------------------------------------------------------
# Data quality audit required by Task C
# -----------------------------------------------------------------------------

function audit_task_c_data(times, x, V)
    nt, nx = size(V)
    nt == length(times) || throw(DimensionMismatch("V rows must match times"))
    nx == length(x) || throw(DimensionMismatch("V columns must match x"))

    dt = diff(times)
    dx = diff(x)

    finite_count = count(isfinite, V)
    total_count = length(V)

    (
        nt = nt,
        nx = nx,
        nan_count = count(isnan, V),
        inf_count = count(isinf, V),
        all_finite = finite_count == total_count,
        duplicate_times = length(unique(times)) != length(times),
        duplicate_x = length(unique(x)) != length(x),
        time_strictly_increasing = all(dt .> 0),
        x_strictly_increasing = all(dx .> 0),
        dt_min = isempty(dt) ? NaN : minimum(dt),
        dt_max = isempty(dt) ? NaN : maximum(dt),
        dt_mean = isempty(dt) ? NaN : mean(dt),
        dx_min = isempty(dx) ? NaN : minimum(dx),
        dx_max = isempty(dx) ? NaN : maximum(dx),
        dx_mean = isempty(dx) ? NaN : mean(dx),
    )
end

function write_task_c_audit(path, audit)
    open(path, "w") do io
        println(io, "Task C data audit")
        println(io, "=================")
        for name in propertynames(audit)
            println(io, name, " = ", getproperty(audit, name))
        end
    end
    path
end

# -----------------------------------------------------------------------------
# Exact time-weighted averaging on an adaptive time grid
# -----------------------------------------------------------------------------

"""Linear interpolation of all profile columns at one requested time."""
function interpolate_profile_at(times, V, tq::Real)
    tq < times[1] && throw(ArgumentError("requested time is before data start"))
    tq > times[end] && throw(ArgumentError("requested time is after data end"))

    j = searchsortedlast(times, tq)
    if j == length(times) || times[j] == tq
        return copy(V[j, :])
    end

    t0, t1 = times[j], times[j + 1]
    λ = (tq - t0) / (t1 - t0)
    (1 - λ) .* V[j, :] .+ λ .* V[j + 1, :]
end

"""
Return a time slice with exact interpolated endpoints t_start and t_end.
"""
function exact_time_slice(times, V, t_start::Real, t_end::Real)
    t_start < t_end || throw(ArgumentError("t_start must be smaller than t_end"))
    t_start >= times[1] || throw(ArgumentError("t_start precedes available data"))
    t_end <= times[end] || throw(ArgumentError("t_end exceeds available data"))

    middle = findall(t -> t > t_start && t < t_end, times)
    ts = vcat(Float64(t_start), times[middle], Float64(t_end))
    rows = Matrix{Float64}(undef, length(ts), size(V, 2))
    rows[1, :] .= interpolate_profile_at(times, V, t_start)
    for (k, i) in enumerate(middle)
        rows[k + 1, :] .= V[i, :]
    end
    rows[end, :] .= interpolate_profile_at(times, V, t_end)
    ts, rows
end

"""
    time_weighted_profile(times, V, t_start, t_end)

Task-D step 1: true trapezoidal time average at every x location.
"""
function time_weighted_profile(times, V, t_start::Real, t_end::Real)
    ts, Vs = exact_time_slice(times, V, t_start, t_end)
    T = t_end - t_start
    nx = size(Vs, 2)
    out = zeros(Float64, nx)

    for j in 1:nx
        out[j] = trapz1(ts, view(Vs, :, j)) / T
    end
    out
end

"""
    task_d_bulk_velocity(x, vbar_t, xlims)

Task-D step 2: signed bulk velocity using the true non-uniform x coordinates
and no-slip wall endpoints.
"""
function task_d_bulk_velocity(x, vbar_t, xlims)
    xfull, vfull = full_wall_profile(x, vbar_t, xlims)
    Lx = xfull[end] - xfull[1]
    trapz1(xfull, vfull) / Lx
end

"""Useful additional diagnostics; not replacements for the Task-D definition."""
function task_d_velocity_diagnostics(x, vbar_t, xlims)
    xfull, vfull = full_wall_profile(x, vbar_t, xlims)
    Lx = xfull[end] - xfull[1]
    Vbulk = trapz1(xfull, vfull) / Lx
    Vabs = trapz1(xfull, abs.(vfull)) / Lx
    Uwind = maximum(abs.(vfull))
    (; Vbulk, Vabs, Uwind)
end

# -----------------------------------------------------------------------------
# Window diagnostics for Task C
# -----------------------------------------------------------------------------

function weighted_profile_delta(x, a, b, xlims)
    xfa, af = full_wall_profile(x, a, xlims)
    _, bf = full_wall_profile(x, b, xlims)
    num = sqrt(trapz1(xfa, (bf .- af) .^ 2))
    den = sqrt(trapz1(xfa, bf .^ 2)) + eps(Float64)
    num / den
end

"""
Compute independent time-weighted window profiles and their bulk velocities.
"""
function independent_window_statistics(
    times,
    x,
    V,
    xlims;
    window::Real = 20.0,
    first_time::Real = times[1],
    last_time::Real = times[end],
)
    window > 0 || throw(ArgumentError("window must be positive"))

    starts = collect(Float64(first_time):Float64(window):(Float64(last_time) - Float64(window)))
    ends = starts .+ window

    profiles = Vector{Vector{Float64}}()
    Vbulk = Float64[]
    Vrms = Float64[]

    for (a, b) in zip(starts, ends)
        p = time_weighted_profile(times, V, a, b)
        push!(profiles, p)
        push!(Vbulk, task_d_bulk_velocity(x, p, xlims))
        push!(Vrms, profile_rms(x, p, xlims))
    end

    deltas = fill(NaN, length(profiles))
    for k in 2:length(profiles)
        deltas[k] = weighted_profile_delta(x, profiles[k - 1], profiles[k], xlims)
    end

    (; starts, ends, profiles, Vbulk, Vrms, deltas)
end

# -----------------------------------------------------------------------------
# Simple periodicity candidate detector based on the profile RMS signal
# -----------------------------------------------------------------------------

function profile_rms_series(times, x, V, xlims)
    [profile_rms(x, view(V, i, :), xlims) for i in eachindex(times)]
end

"""
Detect well-separated local maxima in an irregularly sampled scalar signal.
This is intentionally conservative: it is a screening tool, not a proof of
periodicity.
"""
function detect_cycle_peaks(times, signal; min_separation::Real = 20.0)
    length(times) == length(signal) || throw(DimensionMismatch("times/signal mismatch"))
    length(signal) < 3 && return Int[]

    # Ignore very small peaks by requiring a value above the 70th percentile.
    threshold = quantile(signal, 0.70)
    candidates = Int[]
    for i in 2:(length(signal) - 1)
        if signal[i] >= signal[i - 1] && signal[i] > signal[i + 1] && signal[i] >= threshold
            push!(candidates, i)
        end
    end

    peaks = Int[]
    for i in candidates
        if isempty(peaks) || times[i] - times[peaks[end]] >= min_separation
            push!(peaks, i)
        elseif signal[i] > signal[peaks[end]]
            peaks[end] = i
        end
    end
    peaks
end

"""
Estimate a repeat period from at least three detected peaks.
Returns `(periodic_candidate, period, cv, peak_times)`.
"""
function estimate_periodicity(times, signal; min_separation::Real = 20.0, cv_limit::Real = 0.10)
    peaks = detect_cycle_peaks(times, signal; min_separation)
    pt = times[peaks]
    if length(pt) < 3
        return (; periodic_candidate = false, period = NaN, cv = NaN, peak_times = pt)
    end

    periods = diff(pt)
    pmean = mean(periods)
    cv = std(periods) / (abs(pmean) + eps(Float64))
    (; periodic_candidate = cv <= cv_limit, period = pmean, cv, peak_times = pt)
end

# -----------------------------------------------------------------------------
# Low-Ra analytic reference (optional validation helper)
# -----------------------------------------------------------------------------

"""
Analytic steady 1-D profile for T(x)=0.5-x/H and

    0 = nu*v'' + g_stream*T + drive_force_y,

with v(0)=v(H)=0.
"""
function analytic_low_ra_v(
    x;
    H::Real = 1.0,
    viscosity::Real,
    g_stream::Real = 1.0,
    drive_force_y::Real = 0.0,
)
    A = Float64(g_stream)
    F = Float64(drive_force_y)
    ν = Float64(viscosity)
    Hf = Float64(H)

    @. (
        A * x^3 / (6Hf) -
        (A / 4 + F / 2) * x^2 +
        (A * Hf / 12 + F * Hf / 2) * x
    ) / ν
end
