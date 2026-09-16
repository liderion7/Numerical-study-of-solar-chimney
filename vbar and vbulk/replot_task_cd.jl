# replot_task_cd.jl
# -----------------------------------------------------------------------------
# Read-only redraw of the Task C and Task D figures for ONE case directory.
#
# It never touches the decision outputs. It only READS
#     vbar_xt_raw.csv, case_parameters.txt,
#     task_c/task_c_summary.txt, task_c/task_c_selected_interval.csv,
#     task_c/window_stats.csv, task_c/window_profiles.csv,
#     task_d/vbar_t_profile.csv, task_d_manual/vbar_t_profile.csv
# and it OVERWRITES only these figures (whichever inputs exist):
#     task_c/task_c_velocity_vs_time.png
#     task_c/task_c_periodicity.png
#     task_c/task_c_window_profiles.png
#     task_d/V_X_time_averaged.png
#     task_d_manual/V_X_single_cycle.png
#
# Figure conventions (identical to task_c.jl / task_d.jl, see common/cd_plots.jl):
#   - every profile figure: horizontal = position x (running 1 to 0), vertical = velocity
#   - Task C time-series figures: startup transient / statistically usable phases
#     are shaded differently, and the time at which the flow becomes usable is
#     marked (transient end, t_start, t_end).
#
# Usage:
#   CASE_DIR=data/90degree_Ra1e7_T550 julia replot_task_cd.jl
# -----------------------------------------------------------------------------

using GLMakie
using DelimitedFiles
using Printf

base = @__DIR__
include(joinpath(base, "..", "shared", "common", "cd_tools.jl"))
include(joinpath(base, "..", "shared", "common", "cd_plots.jl"))

const CASE_DIR = get(ENV, "CASE_DIR", joinpath(base, "data", "90degree_Ra1e7_T550"))
const TASKC_DIR = joinpath(CASE_DIR, "task_c")
const RAW_FILE = get(ENV, "VBAR_FILE", joinpath(CASE_DIR, "vbar_xt_raw.csv"))
const PARAM_FILE = joinpath(CASE_DIR, "case_parameters.txt")

const SUMMARY_FILE = joinpath(TASKC_DIR, "task_c_summary.txt")
const SELECTED_FILE = joinpath(TASKC_DIR, "task_c_selected_interval.csv")
const WSTATS_FILE = joinpath(TASKC_DIR, "window_stats.csv")
const WPROF_FILE = joinpath(TASKC_DIR, "window_profiles.csv")

# -----------------------------------------------------------------------------
# Small text parsers (summary / selected-interval files)
# -----------------------------------------------------------------------------
function parse_key_value(path)
    d = Dict{String, String}()
    isfile(path) || return d
    for line in eachline(path)
        occursin('=', line) || continue
        key, value = strip.(split(line, '='; limit = 2))
        haskey(d, key) || (d[key] = value)
    end
    d
end

function kv_float(d, key, default = NaN)
    haskey(d, key) || return default
    v = tryparse(Float64, strip(d[key]))
    isnothing(v) ? default : v
end

function kv_float_list(d, key)
    haskey(d, key) || return Float64[]
    s = replace(strip(d[key]), '[' => "", ']' => "")
    isempty(strip(s)) && return Float64[]
    out = Float64[]
    for tok in split(s, ',')
        v = tryparse(Float64, strip(tok))
        isnothing(v) || push!(out, v)
    end
    out
end

function read_selected_interval(path)
    fallback = (; regime = "unknown", t_start = NaN, t_end = NaN, status = "unknown")
    isfile(path) || return fallback
    raw = readdlm(path, ',', '\n'; skipstart = 1)
    (size(raw, 1) >= 1 && size(raw, 2) >= 6) || return fallback
    (;
        regime = string(raw[1, 1]),
        t_start = Float64(raw[1, 2]),
        t_end = Float64(raw[1, 3]),
        status = string(raw[1, 6]),
    )
end

# -----------------------------------------------------------------------------
# Task C figures (need the raw vbar(x,t) record)
# -----------------------------------------------------------------------------
if isfile(RAW_FILE) && isfile(WPROF_FILE)
    params = parse_key_value(PARAM_FILE)
    xlims = (
        kv_float(params, "x_min", 0.0),
        kv_float(params, "x_max", 1.0),
    )

    println("case dir       = ", CASE_DIR)
    println("raw data       = ", RAW_FILE)
    println("x limits       = ", xlims)

    times, x, V = read_vbar_long_csv(RAW_FILE)
    Vrms_inst = profile_rms_series(times, x, V, xlims)
    Vbulk_inst = [instantaneous_bulk_velocity(x, view(V, i, :), xlims) for i in eachindex(times)]

    summary = parse_key_value(SUMMARY_FILE)
    selected = read_selected_interval(SELECTED_FILE)

    status = haskey(summary, "status") ? strip(summary["status"]) : selected.status
    regime = haskey(summary, "selected regime") ? strip(summary["selected regime"]) : selected.regime
    transient_end = kv_float(summary, "transient_end")
    t_start = isfinite(selected.t_start) ? selected.t_start : kv_float(summary, "t_start")
    t_end = isfinite(selected.t_end) ? selected.t_end : kv_float(summary, "t_end")
    peak_times = kv_float_list(summary, "peak times")

    # Same rule as task_c.jl: a stationary case becomes usable at its first
    # stable window; a periodic case at the end of the startup transient.
    phase_split = regime == "stationary" ? t_start : transient_end

    @printf("status         = %s\n", status)
    @printf("regime         = %s\n", regime)
    @printf("transient end  = %.6f\n", transient_end)
    @printf("t_start        = %.6f\n", t_start)
    @printf("t_end          = %.6f\n", t_end)
    println("phase split    = ", phase_split)
    println("peak times     = ", peak_times)

    wstats = readdlm(WSTATS_FILE, ',', '\n'; skipstart = 1)
    wprof = readdlm(WPROF_FILE, ',', Float64, '\n'; skipstart = 1)
    xp = wprof[:, 1]
    window_profiles = [wprof[:, j] for j in 2:size(wprof, 2)]
    window_starts = size(wstats, 1) == length(window_profiles) ? Float64.(wstats[:, 1]) : fill(NaN, length(window_profiles))
    window_ends = size(wstats, 1) == length(window_profiles) ? Float64.(wstats[:, 2]) : fill(NaN, length(window_profiles))

    plot_task_c_figures(
        TASKC_DIR;
        times,
        Vrms_inst,
        Vbulk_inst,
        x = xp,
        window_starts,
        window_ends,
        window_profiles,
        status,
        phase_split,
        t_start,
        t_end,
        transient_end,
        peak_times,
    )
    println("wrote ", joinpath(TASKC_DIR, "task_c_velocity_vs_time.png"))
    println("wrote ", joinpath(TASKC_DIR, "task_c_periodicity.png"))
    println("wrote ", joinpath(TASKC_DIR, "task_c_window_profiles.png"))
else
    println("skip Task C figures (missing $RAW_FILE or $WPROF_FILE)")
end

# -----------------------------------------------------------------------------
# Task D figures (need the saved time-averaged V-X profile)
# -----------------------------------------------------------------------------
for (subdir, png, title) in (
    ("task_d", "V_X_time_averaged.png", "Task D: time-averaged V-X profile"),
    ("task_d_manual", "V_X_single_cycle.png", "Task D: single-cycle time-averaged V-X profile"),
)
    outdir = joinpath(CASE_DIR, subdir)
    profile_file = joinpath(outdir, "vbar_t_profile.csv")
    if !isfile(profile_file)
        println("skip ", png, " (no ", profile_file, ")")
        continue
    end

    prof = readdlm(profile_file, ',', Float64, '\n'; skipstart = 1)
    xfull = prof[:, 1]
    vfull = prof[:, 2]

    plot_task_d_profile(outdir, png, xfull, vfull; title)
    println("wrote ", joinpath(outdir, png))
end

println("------------------------------------------------------------")
println("Only the figure files listed above were overwritten; every")
println("csv/txt decision output was left untouched.")
println("------------------------------------------------------------")
