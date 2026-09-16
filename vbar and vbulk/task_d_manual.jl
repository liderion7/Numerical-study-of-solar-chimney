# task_d_manual.jl
#
# Manual Task-D calculation when the CFD run cannot be extended further.
#
# Uses one selected mature cycle from the raw Task-C data.
#
# Current 90-degree recommended interval:
#   t_start = 281.6001522614686
#   t_end   = 443.4000188476378

using IncompressibleNavierStokes
using GLMakie
using Printf

base = @__DIR__

include(
    joinpath(
        base,
        "..",
        "shared",
        "common",
        "cd_tools.jl",
    ),
)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

const CASE_DIR =
    get(
        ENV,
        "CASE_DIR",
        joinpath(
            base,
            "data",
            "90degree_Ra1e7_T550",
        ),
    )

const INPUT =
    get(
        ENV,
        "TASKD_FILE",
        joinpath(
            CASE_DIR,
            "vbar_xt_raw.csv",
        ),
    )

const OUTDIR =
    joinpath(
        CASE_DIR,
        "task_d_manual",
    )

mkpath(OUTDIR)

# ------------------------------------------------------------
# Manual averaging interval
# ------------------------------------------------------------

const T_START =
    parse(
        Float64,
        get(
            ENV,
            "T_START",
            "281.6001522614686",
        ),
    )

const T_END =
    parse(
        Float64,
        get(
            ENV,
            "T_END",
            "443.4000188476378",
        ),
    )

# ------------------------------------------------------------
# Geometry
# ------------------------------------------------------------

const XMIN = 0.0
const XMAX = 1.0
const XLIMS = (XMIN, XMAX)

# ------------------------------------------------------------
# Read raw Task-C data
# ------------------------------------------------------------

println("Reading: ", INPUT)

times, x, V =
    read_vbar_long_csv(INPUT)

audit =
    audit_task_c_data(
        times,
        x,
        V,
    )

audit.all_finite ||
    error("NaN / Inf found")

audit.time_strictly_increasing ||
    error("Time is not strictly increasing")

audit.x_strictly_increasing ||
    error("x is not strictly increasing")

T_START >= times[1] ||
    error("T_START is before available CFD data")

T_END <= times[end] ||
    error("T_END exceeds available CFD data")

# ------------------------------------------------------------
# Task D step 1:
# true time-weighted average
# ------------------------------------------------------------

vbar_t =
    time_weighted_profile(
        times,
        V,
        T_START,
        T_END,
    )

# ------------------------------------------------------------
# Task D step 2:
# spatial integration
# ------------------------------------------------------------

diag =
    task_d_velocity_diagnostics(
        x,
        vbar_t,
        XLIMS,
    )

# Include the no-slip wall points explicitly
xfull, vfull =
    full_wall_profile(
        x,
        vbar_t,
        XLIMS,
    )

flow_area =
    trapz1(
        xfull,
        vfull,
    )

# ------------------------------------------------------------
# Save final V-X profile
# ------------------------------------------------------------

open(
    joinpath(
        OUTDIR,
        "vbar_t_profile.csv",
    ),
    "w",
) do io

    println(io, "x,vbar_t")

    for i in eachindex(xfull)
        println(
            io,
            xfull[i],
            ',',
            vfull[i],
        )
    end
end

# ------------------------------------------------------------
# Save result
# ------------------------------------------------------------

open(
    joinpath(
        OUTDIR,
        "task_d_manual_result.csv",
    ),
    "w",
) do io

    println(
        io,
        "t_start,t_end,Vbulk,flow_area,Vabs,Uwind",
    )

    println(
        io,
        T_START, ',',
        T_END, ',',
        diag.Vbulk, ',',
        flow_area, ',',
        diag.Vabs, ',',
        diag.Uwind,
    )
end

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

open(
    joinpath(
        OUTDIR,
        "task_d_manual_summary.txt",
    ),
    "w",
) do io

    println(io, "Task D manual single-cycle estimate")
    println(io, "===================================")

    println(io)
    println(io, "Averaging interval:")
    println(io, "t_start = ", T_START)
    println(io, "t_end   = ", T_END)
    println(io, "duration = ", T_END - T_START)

    println(io)
    println(io, "Primary result:")
    println(io, "Vbulk = ", diag.Vbulk)

    println(io)
    println(io, "Additional diagnostics:")
    println(io, "flow area = ", flow_area)
    println(io, "Vabs      = ", diag.Vabs)
    println(io, "Uwind     = ", diag.Uwind)

    println(io)
    println(
        io,
        "NOTE: This result is based on one mature-cycle candidate.",
    )
    println(
        io,
        "It should be reported as a provisional single-cycle estimate,",
    )
    println(
        io,
        "not as a fully statistically converged value.",
    )
end

# ------------------------------------------------------------
# Plot
#
# Shared Task C / Task D axis convention:
#   horizontal = position x (running 1 to 0), vertical = velocity.
# ------------------------------------------------------------

include(
    joinpath(
        base,
        "..",
        "shared",
        "common",
        "cd_plots.jl",
    ),
)

plot_task_d_profile(
    OUTDIR,
    "V_X_single_cycle.png",
    xfull,
    vfull;
    title = "Task D: single-cycle time-averaged V-X profile",
)

# ------------------------------------------------------------
# Terminal output
# ------------------------------------------------------------

println()
println("------------------------------------------")
println("Task D manual result")
println("------------------------------------------")

@printf(
    "Interval  = [%.6f, %.6f]\n",
    T_START,
    T_END,
)

@printf(
    "Duration  = %.6f\n",
    T_END - T_START,
)

@printf(
    "Vbulk     = %.8f\n",
    diag.Vbulk,
)

@printf(
    "Vabs      = %.8f\n",
    diag.Vabs,
)

@printf(
    "Uwind     = %.8f\n",
    diag.Uwind,
)

println("------------------------------------------")
println("Output: ", OUTDIR)