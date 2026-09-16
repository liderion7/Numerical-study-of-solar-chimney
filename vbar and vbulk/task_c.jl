# task_c.jl
# -----------------------------------------------------------------------------
# Task C: decide whether vbar(x,t) has entered a reliable statistical regime.
#
# It does NOT calculate the final bulk velocity. It audits the CFD data,
# examines transient/stationary/periodic behaviour, produces V-X diagnostics,
# and, only when justified, exports a clean time interval for Task D.
# -----------------------------------------------------------------------------

using IncompressibleNavierStokes
using GLMakie
using Printf
using Statistics

base = @__DIR__
include(joinpath(base, "..", "shared", "common", "cd_tools.jl"))
include(joinpath(base, "..", "shared", "common", "cd_plots.jl"))

# Default: the most complete production 90-degree case (main_task_cd/data).
const CASE_DIR = get(ENV, "CASE_DIR", joinpath(base, "data", "90degree_Ra1e7_T550"))
const INPUT = get(ENV, "VBAR_FILE", joinpath(CASE_DIR, "vbar_xt_raw.csv"))
const PARAM_FILE = joinpath(CASE_DIR, "case_parameters.txt")
const OUTDIR = joinpath(CASE_DIR, "task_c")
mkpath(OUTDIR)

const WINDOW = parse(Float64, get(ENV, "WINDOW", "20.0"))
const DELTA_LIMIT = parse(Float64, get(ENV, "DELTA_LIMIT", "0.05"))
const RMS_DRIFT_LIMIT = parse(Float64, get(ENV, "RMS_DRIFT_LIMIT", "0.05"))
const N_STABLE_WINDOWS = parse(Int, get(ENV, "N_STABLE_WINDOWS", "3"))
const PERIOD_CV_LIMIT = parse(Float64, get(ENV, "PERIOD_CV_LIMIT", "0.10"))
const PROFILE_CORR_LIMIT = parse(Float64, get(ENV, "PROFILE_CORR_LIMIT", "0.98"))
const CYCLE_VBULK_REL_LIMIT = parse(Float64, get(ENV, "CYCLE_VBULK_REL_LIMIT", "0.10"))

# Optional manual override. If both are supplied, Task C uses them after audit.
const MANUAL_T_START = haskey(ENV, "T_START") ? parse(Float64, ENV["T_START"]) : nothing
const MANUAL_T_END = haskey(ENV, "T_END_STAT") ? parse(Float64, ENV["T_END_STAT"]) : nothing

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

"""Minimum correlation between consecutive same-phase peak profiles."""
function peak_profile_correlation(times, V, peaks)
    length(peaks) < 2 && return NaN
    c = Float64[]
    for i in 1:(length(peaks) - 1)
        p1 = interpolate_profile_at(times, V, peaks[i])
        p2 = interpolate_profile_at(times, V, peaks[i + 1])
        push!(c, cor(p1, p2))
    end
    minimum(c)
end

"""Maximum relative difference between consecutive full-cycle averaged Vbulk."""
function cycle_vbulk_relative_difference(times, x, V, xlims, peaks)
    length(peaks) < 3 && return NaN
    vb = Float64[]
    for i in 1:(length(peaks) - 1)
        p = time_weighted_profile(times, V, peaks[i], peaks[i + 1])
        push!(vb, instantaneous_bulk_velocity(x, p, xlims))
    end
    diffs = [
        abs(vb[i + 1] - vb[i]) / max(abs(vb[i + 1]), eps(Float64))
        for i in 1:(length(vb) - 1)
    ]
    maximum(diffs)
end

params = read_parameter_file(PARAM_FILE)
xmin = parse(Float64, get(params, "x_min", "0.0"))
xmax = parse(Float64, get(params, "x_max", "1.0"))
xlims = (xmin, xmax)

times, x, V = read_vbar_long_csv(INPUT)
audit = audit_task_c_data(times, x, V)
write_task_c_audit(joinpath(OUTDIR, "task_c_data_audit.txt"), audit)

println("Task C audit")
println("------------")
for name in propertynames(audit)
    println(name, " = ", getproperty(audit, name))
end

# Hard data-quality gate.
data_ok = audit.all_finite &&
          !audit.duplicate_times &&
          !audit.duplicate_x &&
          audit.time_strictly_increasing &&
          audit.x_strictly_increasing

if !data_ok
    error("Task C failed data audit. See task_c_data_audit.txt")
end

# Instantaneous one-dimensional diagnostics from vbar(x,t).
Vbulk_inst = [instantaneous_bulk_velocity(x, view(V, i, :), xlims) for i in eachindex(times)]
Vrms_inst = profile_rms_series(times, x, V, xlims)

# Independent time-weighted windows.
wins = independent_window_statistics(
    times,
    x,
    V,
    xlims;
    window = WINDOW,
    first_time = times[1],
    last_time = times[end],
)

# -----------------------------------------------------------------------------
# Stationary-regime test: require several final independent windows to agree.
# -----------------------------------------------------------------------------
stationary = false
stationary_start = NaN
if length(wins.profiles) >= N_STABLE_WINDOWS
    first_k = length(wins.profiles) - N_STABLE_WINDOWS + 1
    relevant_delta = wins.deltas[(first_k + 1):end]

    rms_rel = Float64[]
    for k in (first_k + 1):length(wins.Vrms)
        denom = max(abs(wins.Vrms[k]), eps(Float64))
        push!(rms_rel, abs(wins.Vrms[k] - wins.Vrms[k - 1]) / denom)
    end

    stationary = all(relevant_delta .<= DELTA_LIMIT) && all(rms_rel .<= RMS_DRIFT_LIMIT)
    stationary && (stationary_start = wins.starts[first_k])
end

# -----------------------------------------------------------------------------
# Periodicity screening based on profile RMS. A reliable automatic periodic
# interval needs >= 3 same-phase POST-TRANSIENT peaks -> at least two complete
# mature periods. The first detected peak is treated as the end of the startup
# transient so that startup dynamics are not mixed with mature cycles.
# -----------------------------------------------------------------------------
min_sep = max(3WINDOW, 20.0)
all_peaks = detect_cycle_peaks(times, Vrms_inst; min_separation = min_sep)
transient_end = isempty(all_peaks) ? times[end] : times[all_peaks[1]] + min_sep

post = findall(t -> t >= transient_end, times)
if length(post) >= 3
    period_info = estimate_periodicity(
        times[post],
        Vrms_inst[post];
        min_separation = min_sep,
        cv_limit = PERIOD_CV_LIMIT,
    )
else
    period_info = (; periodic_candidate = false, period = NaN, cv = NaN, peak_times = Float64[])
end

profile_corr = peak_profile_correlation(times, V, period_info.peak_times)
cycle_vbulk_rel = cycle_vbulk_relative_difference(times, x, V, xlims, period_info.peak_times)

regime = "not_ready"
status = "not_ready"
t_start = NaN
t_end = NaN
ready = false

periodic_ok =
    period_info.periodic_candidate &&
    length(period_info.peak_times) >= 3 &&
    isfinite(profile_corr) &&
    profile_corr > PROFILE_CORR_LIMIT &&
    isfinite(cycle_vbulk_rel) &&
    cycle_vbulk_rel < CYCLE_VBULK_REL_LIMIT

if !isnothing(MANUAL_T_START) || !isnothing(MANUAL_T_END)
    (!isnothing(MANUAL_T_START) && !isnothing(MANUAL_T_END)) ||
        error("Manual selection requires both T_START and T_END_STAT")
    MANUAL_T_START < MANUAL_T_END || error("T_START must be smaller than T_END_STAT")
    t_start = MANUAL_T_START
    t_end = MANUAL_T_END
    regime = "manual"
    status = "manual"
    ready = true
elseif stationary
    t_start = stationary_start
    t_end = times[end]
    regime = "stationary"
    status = "converged_stationary"
    ready = true
elseif periodic_ok
    # Last two complete peak-to-peak cycles, same phase at both ends.
    t_start = period_info.peak_times[end - 2]
    t_end = period_info.peak_times[end]
    regime = "periodic"
    status = "converged_periodic"
    ready = true
elseif length(period_info.peak_times) >= 2
    # Fallback for limited compute: a single mature peak-to-peak cycle.
    # This is explicitly provisional, never reported as converged.
    t_start = period_info.peak_times[end - 1]
    t_end = period_info.peak_times[end]
    regime = "provisional_single_cycle"
    status = "provisional_single_cycle"
    ready = true
end

# -----------------------------------------------------------------------------
# Window data output
# -----------------------------------------------------------------------------
open(joinpath(OUTDIR, "window_stats.csv"), "w") do io
    println(io, "t_start,t_end,Vbulk,Vrms,profile_delta")
    for k in eachindex(wins.starts)
        println(
            io,
            wins.starts[k], ',', wins.ends[k], ',', wins.Vbulk[k], ',',
            wins.Vrms[k], ',', wins.deltas[k],
        )
    end
end

# Save window profiles: x plus one column per independent window.
open(joinpath(OUTDIR, "window_profiles.csv"), "w") do io
    print(io, "x")
    for k in eachindex(wins.starts)
        print(io, ",vbar_", wins.starts[k], "_", wins.ends[k])
    end
    println(io)
    for j in eachindex(x)
        print(io, x[j])
        for p in wins.profiles
            print(io, ',', p[j])
        end
        println(io)
    end
end

# -----------------------------------------------------------------------------
# Plots required for Task C report
# -----------------------------------------------------------------------------
# The velocity-vs-time, V-X profile and periodicity figures use the shared
# conventions of common/cd_plots.jl: the profile figure follows the agreed axis
# convention (horizontal = position x running 1 to 0, vertical = velocity), and the
# two time-series figures shade the startup transient / statistically usable
# phases and mark the time at which the flow becomes statistically usable.
phase_split = regime == "stationary" ? t_start : transient_end

plot_task_c_figures(
    OUTDIR;
    times,
    Vrms_inst,
    Vbulk_inst,
    x,
    window_starts = wins.starts,
    window_ends = wins.ends,
    window_profiles = wins.profiles,
    status,
    phase_split,
    t_start,
    t_end,
    transient_end,
    peak_times = period_info.peak_times,
)

# Window-to-window profile change (kept as in the original implementation).
fig3 = Figure(size = (1000, 650))
ax3 = Axis(fig3[1, 1];
    title = "Task C: window-to-window profile change",
    xlabel = "window end time",
    ylabel = "relative profile delta",
)
length(wins.ends) >= 2 && lines!(ax3, wins.ends[2:end], wins.deltas[2:end])
hlines!(ax3, [DELTA_LIMIT]; linestyle = :dash)
save(joinpath(OUTDIR, "task_c_profile_delta.png"), fig3)

# -----------------------------------------------------------------------------
# Clean selected data for Task D
# -----------------------------------------------------------------------------
if ready
    ts_clean, V_clean = exact_time_slice(times, V, t_start, t_end)
    clean_profiles = [vec(copy(V_clean[i, :])) for i in 1:size(V_clean, 1)]
    write_vbar_long_csv(
        joinpath(OUTDIR, "task_c_clean_vbar.csv"),
        ts_clean,
        x,
        clean_profiles,
    )

    open(joinpath(OUTDIR, "task_c_selected_interval.csv"), "w") do io
        println(io, "regime,t_start,t_end,estimated_period,period_cv,status,profile_corr,cycle_vbulk_rel")
        println(
            io,
            regime, ',', t_start, ',', t_end, ',',
            period_info.period, ',', period_info.cv, ',',
            status, ',', profile_corr, ',', cycle_vbulk_rel,
        )
    end
end

# -----------------------------------------------------------------------------
# Human-readable summary
# -----------------------------------------------------------------------------
open(joinpath(OUTDIR, "task_c_summary.txt"), "w") do io
    println(io, "Task C summary")
    println(io, "==============")
    println(io, "input = ", INPUT)
    println(io, "data audit passed = ", data_ok)
    println(io, "stationary candidate = ", stationary)
    println(io, "periodic candidate = ", period_info.periodic_candidate)
    println(io, "transient_end = ", transient_end)
    println(io, "estimated period = ", period_info.period)
    println(io, "period coefficient of variation = ", period_info.cv)
    println(io, "peak times = ", collect(period_info.peak_times))
    println(io, "profile correlation (peak profiles) = ", profile_corr)
    println(io, "cycle Vbulk relative difference = ", cycle_vbulk_rel)
    println(io, "selected regime = ", regime)
    println(io, "status = ", status)
    println(io, "Task D ready = ", ready)
    if ready
        println(io, "t_start = ", t_start)
        println(io, "t_end = ", t_end)
        println(io, "clean file = task_c_clean_vbar.csv")
    else
        println(io, "No statistically justified interval was selected automatically.")
        println(io, "Extend the CFD run or inspect the diagnostics before Task D.")
    end
end

println("------------------------------------------------------------")
println("Task C regime: ", regime)
println("Task C status: ", status)
println("Task D ready: ", ready)
if ready
    println("Selected interval: [", t_start, ", ", t_end, "]")
else
    println("No valid interval yet. Extend the run or use a justified manual override.")
end
println("Outputs: ", OUTDIR)
println("------------------------------------------------------------")
