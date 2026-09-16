# common/cd_plots.jl
# -----------------------------------------------------------------------------
# Shared Task C / Task D figure generation.
#
# Axis convention for EVERY profile figure (window profiles and the V-X profile).
# The profile is drawn as if the original "velocity on the horizontal axis"
# figure had been rotated 90° counterclockwise, with the text left upright:
#
#     horizontal axis = wall-normal position x, axis title "x axis"
#                       running from x = 1 (left) to x = 0 (right)
#     vertical   axis = velocity vbar_t, axis title "z axis"
#
# Set PROFILE_POSITION_AXIS_REVERSED = false to draw the position axis the
# normal way round (x = 0 on the left) while keeping everything else.
#
# Task C time-series figures additionally mark when the flow becomes
# statistically usable:
#     - "transient end"      : end of the startup transient (periodic screening)
#     - "t_start" / "t_end"  : the averaging interval Task C selected
# and shade the startup transient phase and the usable phase with different
# background colours.
# -----------------------------------------------------------------------------

using GLMakie
using Printf

const AXIS_VELOCITY_TITLE = "z axis"
const AXIS_POSITION_TITLE = "x axis"

const PROFILE_POSITION_AXIS_REVERSED = true

const TRANSIENT_BAND_COLOUR = (:darkorange, 0.10)
const USABLE_BAND_COLOUR = (:seagreen, 0.08)
const SELECTED_BAND_COLOUR = (:seagreen, 0.22)

"""
    use_chinese_font()

Select a Chinese-capable font so the Chinese axis titles render correctly.
Falls back to the Makie default font when no candidate file exists.
"""
function use_chinese_font()
    for f in (
        "C:/Windows/Fonts/simhei.ttf",
        "C:/Windows/Fonts/msyh.ttc",
        "C:/Windows/Fonts/simsun.ttc",
    )
        if isfile(f)
            set_theme!(fonts = (; regular = f, bold = f))
            return f
        end
    end
    nothing
end

"""
    apply_profile_axis_limits!(ax, xvalues)

Fix the horizontal (position) axis of a profile figure. After the 90°
counterclockwise rotation the position runs from 1 on the left to 0 on the
right, which is what `PROFILE_POSITION_AXIS_REVERSED` selects.
"""
function apply_profile_axis_limits!(ax, xvalues)
    xmin = minimum(xvalues)
    xmax = maximum(xvalues)
    if PROFILE_POSITION_AXIS_REVERSED
        xlims!(ax, xmax, xmin)
    else
        xlims!(ax, xmin, xmax)
    end
    nothing
end

"""
    add_phase_bands!(ax, times, phase_split, t_start, t_end)

Shade `[t0, phase_split]` as the startup transient, `[phase_split, t1]` as the
statistically usable phase, and `[t_start, t_end]` as the interval Task C
actually selected for averaging.
"""
function add_phase_bands!(ax, times, phase_split, t_start, t_end)
    t0 = times[1]
    t1 = times[end]

    if isfinite(phase_split) && t0 < phase_split < t1
        vspan!(ax, t0, phase_split; color = TRANSIENT_BAND_COLOUR, label = "startup transient")
        vspan!(ax, phase_split, t1; color = USABLE_BAND_COLOUR, label = "statistically usable phase")
    elseif !isfinite(phase_split)
        # No split known: the whole record is treated as usable.
        vspan!(ax, t0, t1; color = USABLE_BAND_COLOUR, label = "statistically usable phase")
    end

    if isfinite(t_start) && isfinite(t_end) && t_start < t_end
        vspan!(ax, t_start, t_end; color = SELECTED_BAND_COLOUR, label = "selected averaging interval")
    end
    nothing
end

"""Vertical marker line with the numerical time written next to it."""
function mark_event!(ax, t, name, colour, linestyle, ytop)
    isfinite(t) || return nothing
    label = @sprintf("%s = %.1f", name, t)
    vlines!(ax, [t]; color = colour, linestyle = linestyle, linewidth = 2, label = label)
    text!(
        ax,
        t,
        ytop;
        text = string(round(t; digits = 1)),
        align = (:left, :top),
        offset = (4, -4),
        fontsize = 13,
        color = colour,
    )
    nothing
end

"""
    plot_task_c_figures(outdir; ...)

Write the three Task C figures:
    task_c_velocity_vs_time.png, task_c_window_profiles.png, task_c_periodicity.png
"""
function plot_task_c_figures(
    outdir;
    times,
    Vrms_inst,
    Vbulk_inst,
    x,
    window_starts,
    window_ends,
    window_profiles,
    status = "",
    phase_split = NaN,
    t_start = NaN,
    t_end = NaN,
    transient_end = NaN,
    peak_times = Float64[],
)
    use_chinese_font()
    mkpath(outdir)

    # ---- velocity statistics versus time ------------------------------------
    fig1 = Figure(size = (1100, 700))
    ax1 = Axis(
        fig1[1, 1];
        title = "Task C: velocity statistics versus time ($status)",
        xlabel = "t",
        ylabel = "velocity",
    )
    add_phase_bands!(ax1, times, phase_split, t_start, t_end)
    lines!(ax1, times, Vrms_inst; label = "profile RMS")
    lines!(ax1, times, Vbulk_inst; label = "instantaneous signed Vbulk")
    ytop = maximum(Vrms_inst)
    mark_event!(ax1, transient_end, "transient end", :gray30, :dash, ytop)
    mark_event!(ax1, t_start, "t_start", :seagreen, :solid, ytop)
    mark_event!(ax1, t_end, "t_end", :seagreen, :dashdot, ytop)
    axislegend(ax1; position = :lb)
    save(joinpath(outdir, "task_c_velocity_vs_time.png"), fig1)

    # ---- independent window-averaged V-X profiles ---------------------------
    # Data layout is unchanged: position on the horizontal axis (1 → 0),
    # velocity on the vertical axis. Only the axis TITLES are swapped, as
    # requested for the report (horizontal axis named "z axis", vertical "x axis").
    fig2 = Figure(size = (1000, 700))
    ax2 = Axis(
        fig2[1, 1];
        title = "Task C: independent window-averaged V-X profiles",
        xlabel = AXIS_VELOCITY_TITLE,
        ylabel = AXIS_POSITION_TITLE,
    )
    # Plot at most ~12 representative windows to keep the figure readable.
    stride = max(1, cld(length(window_profiles), 12))
    for k in 1:stride:length(window_profiles)
        lines!(
            ax2,
            x,
            window_profiles[k];
            label = @sprintf("%.1f-%.1f", window_starts[k], window_ends[k]),
        )
    end
    apply_profile_axis_limits!(ax2, x)
    isempty(window_profiles) || axislegend(ax2; position = :rt)
    save(joinpath(outdir, "task_c_window_profiles.png"), fig2)

    # ---- periodicity screening signal ---------------------------------------
    fig3 = Figure(size = (1100, 700))
    ax3 = Axis(
        fig3[1, 1];
        title = "Task C: periodicity screening signal ($status)",
        xlabel = "t",
        ylabel = "profile RMS",
    )
    add_phase_bands!(ax3, times, phase_split, t_start, t_end)
    lines!(ax3, times, Vrms_inst; label = "profile RMS")
    if !isempty(peak_times)
        peak_indices = [argmin(abs.(times .- tp)) for tp in peak_times]
        scatter!(ax3, times[peak_indices], Vrms_inst[peak_indices]; label = "detected peaks")
    end
    mark_event!(ax3, transient_end, "transient end", :gray30, :dash, ytop)
    mark_event!(ax3, t_start, "t_start", :seagreen, :solid, ytop)
    mark_event!(ax3, t_end, "t_end", :seagreen, :dashdot, ytop)
    axislegend(ax3; position = :lb)
    save(joinpath(outdir, "task_c_periodicity.png"), fig3)

    nothing
end

"""
    plot_task_d_profile(outdir, filename, xfull, vfull; title)

Write one Task D V-X figure with the shared profile orientation
(horizontal = position x running 1 → 0, vertical = velocity vbar_t).
"""
function plot_task_d_profile(
    outdir,
    filename,
    xfull,
    vfull;
    title = "Task D: time-averaged V-X profile",
)
    use_chinese_font()
    mkpath(outdir)

    fig = Figure(size = (1000, 700))
    ax = Axis(
        fig[1, 1];
        title = title,
        xlabel = AXIS_POSITION_TITLE,
        ylabel = AXIS_VELOCITY_TITLE,
    )
    lines!(ax, xfull, vfull; linewidth = 3)
    hlines!(ax, [0.0]; linestyle = :dash)
    apply_profile_axis_limits!(ax, xfull)
    save(joinpath(outdir, filename), fig)

    nothing
end
