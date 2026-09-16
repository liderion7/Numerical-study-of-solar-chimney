using GLMakie
using IncompressibleNavierStokes
using Statistics


# ============================================================
# 1. Tilted Boussinesq force
#    Wang: buoyancy direction = (sinβ, cosβ)
# ============================================================

function boussinesq_tilt!(
    force,
    state,
    t;
    setup,
    cache,
    viscosity,
    conductivity,
    gravity,
    beta,
    dodissipation,
)
    θ = deg2rad(beta)

    gx = gravity * sin(θ)
    gy = gravity * cos(θ)

    # Standard Boussinesq terms + y-direction buoyancy
    IncompressibleNavierStokes.boussinesq!(
        force,
        state,
        t;
        setup,
        cache,
        viscosity,
        conductivity,
        gdir = 2,
        gravity = gy,
        dodissipation,
    )

    # Add x-direction buoyancy
    IncompressibleNavierStokes.applygravity!(
        force.u,
        state.temp,
        setup,
        1,
        gx,
    )

    nothing
end


# Same adaptive timestep rule as built-in Boussinesq model
function IncompressibleNavierStokes.propose_timestep(
    ::typeof(boussinesq_tilt!),
    state,
    setup,
    params,
)
    IncompressibleNavierStokes.propose_timestep(
        IncompressibleNavierStokes.boussinesq!,
        state,
        setup,
        params,
    )
end


# ============================================================
# 2. Helper
#    Find grid index closest to physical coordinate
# ============================================================

function nearest_index(grid, target)
    argmin(abs.(grid .- target))
end


# ============================================================
# 3. Wang Fig.2-style plot
#    Temperature + velocity + mesh + gravity arrow
# ============================================================

function fig2plot(
    state;
    setup,
    beta = 0.0,
    arrow_count = (17, 13),
    mesh_stride = (4, 4),
)
    state isa Observable || (state = Observable(state))

    (; xp, Ip) = setup
    xf = Array.(getindex.(xp, Ip.indices))

    xplot = xf[1]
    yplot = xf[2]

    temp = IncompressibleNavierStokes.observefield(
        state;
        setup,
        fieldname = :temperature,
    )

    velocity = IncompressibleNavierStokes.observefield(
        state;
        setup,
        fieldname = :velocity,
    )


    # --------------------------------------------------------
    # Uniform velocity-arrow positions in physical space
    # --------------------------------------------------------

    nx_arrow, ny_arrow = arrow_count

    x_targets = range(
        first(xplot),
        last(xplot),
        length = nx_arrow,
    )

    y_targets = range(
        first(yplot),
        last(yplot),
        length = ny_arrow,
    )

    ix = [
        nearest_index(xplot, x)
        for x in x_targets
    ]

    iy = [
        nearest_index(yplot, y)
        for y in y_targets
    ]

    x_arrow = xplot[ix]
    y_arrow = yplot[iy]


    # --------------------------------------------------------
    # Velocity vectors
    # --------------------------------------------------------

    scaled_velocity = lift(velocity) do V
        U = Array(V[ix, iy, 1])
        VY = Array(V[ix, iy, 2])

        speed = sqrt.(U .^ 2 .+ VY .^ 2)
        umax = maximum(speed)

        if !isfinite(umax) || umax < 1e-12
            return zeros(size(U)), zeros(size(VY))
        end

        scale = 0.055 / umax

        U .* scale, VY .* scale
    end

    uq = lift(v -> v[1], scaled_velocity)
    vq = lift(v -> v[2], scaled_velocity)


    # --------------------------------------------------------
    # Figure
    # --------------------------------------------------------

    fig = Figure(size = (760, 420))

    ax = Axis(
        fig[1, 1];
        title = "β = $(beta)°",
        xlabel = "x",
        ylabel = "y",
        aspect = DataAspect(),
    )


    # Temperature field
    hm = heatmap!(
        ax,
        xplot,
        yplot,
        temp;
        colormap = :jet,
        colorrange = (0.0, 1.0),
    )


    # --------------------------------------------------------
    # Mesh
    # --------------------------------------------------------

    mesh_stride_x, mesh_stride_y = mesh_stride

    for i in 1:mesh_stride_x:length(xplot)
        x = xplot[i]

        lines!(
            ax,
            [x, x],
            [first(yplot), last(yplot)];
            color = (:black, 0.22),
            linewidth = 0.5,
        )
    end

    for j in 1:mesh_stride_y:length(yplot)
        y = yplot[j]

        lines!(
            ax,
            [first(xplot), last(xplot)],
            [y, y];
            color = (:black, 0.22),
            linewidth = 0.5,
        )
    end


    # Outer boundaries
    lines!(
        ax,
        [first(xplot), first(xplot)],
        [first(yplot), last(yplot)];
        color = :black,
        linewidth = 0.8,
    )

    lines!(
        ax,
        [last(xplot), last(xplot)],
        [first(yplot), last(yplot)];
        color = :black,
        linewidth = 0.8,
    )

    lines!(
        ax,
        [first(xplot), last(xplot)],
        [first(yplot), first(yplot)];
        color = :black,
        linewidth = 0.8,
    )

    lines!(
        ax,
        [first(xplot), last(xplot)],
        [last(yplot), last(yplot)];
        color = :black,
        linewidth = 0.8,
    )


    # --------------------------------------------------------
    # Velocity arrows
    # --------------------------------------------------------

    arrows2d!(
        ax,
        x_arrow,
        y_arrow,
        uq,
        vq;
        color = :black,
        lengthscale = 1.0,
        align = :center,
    )


    # --------------------------------------------------------
    # Gravity arrow
    # --------------------------------------------------------

    θ = deg2rad(beta)
    Lg = 0.23

    arrows2d!(
        ax,
        [0.0],
        [0.58],
        [-Lg * sin(θ)],
        [-Lg * cos(θ)];
        color = :black,
        lengthscale = 1.0,
        align = :tail,
    )


    xlims!(ax, -1.0, 1.0)
    ylims!(ax, 0.0, 1.0)

    Colorbar(
        fig[1, 2],
        hm;
        ticks = [0.0, 0.5, 1.0],
        label = "T",
    )

    fig
end


# ============================================================
# 4. Global Nu history
#
# Store Nu(t) so that after the simulation we can calculate
# the statistically stable mean value.
# ============================================================

nu_times = Float64[]
nu_bottom_values = Float64[]
nu_top_values = Float64[]


# ============================================================
# 5. Nusselt number
#
# Wang:
# Nu is spatially averaged over the hot plate.
#
# Nu(t) = -(1/W) ∫ (dT/dy)_wall dx
#
# Bottom and top are both calculated.
# ============================================================

function nusseltplot(state; setup)

    state isa Observable || (state = Observable(state))

    (; Δ, Δu) = setup

    Δy_bottom = Δu[2][1]
    Δy_top = Δu[2][end - 1]

    Nu_bottom = Observable(Point2f[])
    Nu_top = Observable(Point2f[])

    on(state) do (; temp, t)

        # ----------------------------------------------------
        # Bottom hot wall
        # ----------------------------------------------------

        dTdy_bottom =
            @. (temp[:, 2] - temp[:, 1]) / Δy_bottom

        wx = Δ[1][2:end-1]

        Nu_b =
            sum(
                (-dTdy_bottom[2:end-1]) .* wx
            ) / sum(wx)


        # ----------------------------------------------------
        # Top cold wall
        # ----------------------------------------------------

        dTdy_top =
            @. (temp[:, end] - temp[:, end-1]) / Δy_top

        Nu_t =
            sum(
                (-dTdy_top[2:end-1]) .* wx
            ) / sum(wx)


        # ----------------------------------------------------
        # Save history
        # ----------------------------------------------------

        push!(nu_times, Float64(t))
        push!(nu_bottom_values, Float64(Nu_b))
        push!(nu_top_values, Float64(Nu_t))


        # ----------------------------------------------------
        # Plot history
        # ----------------------------------------------------

        push!(
            Nu_bottom[],
            Point2f(t, Nu_b),
        )

        push!(
            Nu_top[],
            Point2f(t, Nu_t),
        )


        notify(Nu_bottom)
        notify(Nu_top)
    end


    fig = Figure(size = (850, 480))

    ax = Axis(
        fig[1, 1];
        title = "Nusselt number - β = 0°",
        xlabel = "t / t_f",
        ylabel = "Nu",
    )


    lines!(
        ax,
        Nu_bottom;
        label = "Bottom hot plate",
    )

    lines!(
        ax,
        Nu_top;
        label = "Top cold plate",
    )


    # Wang Γ = 2, β = 0°, double-roll result
    hlines!(
        ax,
        [12.13];
        linestyle = :dash,
        label = "Wang DRS: Nu = 12.13",
    )


    axislegend(ax)

    on(_ -> autolimits!(ax), Nu_bottom)

    fig
end


# ============================================================
# 6. Wang parameters
# ============================================================

Ra = 1.0e7
Pr = 0.71

ν = sqrt(Pr / Ra)
κ = 1.0 / sqrt(Ra * Pr)

beta = 0.0


println("Ra = ", Ra)
println("Pr = ", Pr)
println("ν  = ", ν)
println("κ  = ", κ)


# ============================================================
# 7. Non-uniform near-wall mesh
#
# Gamma = W/H = 2
#
# x = [-1, 1]
# y = [0, 1]
# ============================================================

n = 128

xgrid = tanh_grid(
    -1.0,
    1.0,
    2n,
    1.2,
)

ygrid = tanh_grid(
    0.0,
    1.0,
    n,
    1.2,
)


# ============================================================
# 8. Boundary conditions
# ============================================================

setup = Setup(;

    x = (
        xgrid,
        ygrid,
    ),

    boundary_conditions = (;

        # No-slip on all walls
        u = (
            (DirichletBC(), DirichletBC()),
            (DirichletBC(), DirichletBC()),
        ),

        # Left/right insulated
        # Bottom T = 1
        # Top T = 0
        temp = (
            (SymmetricBC(), SymmetricBC()),
            (DirichletBC(1.0), DirichletBC(0.0)),
        ),
    ),
)

# ============================================================
# Pressure solver for non-uniform tanh grid
# ============================================================

psolver = default_psolver(setup)


# ============================================================
# 9. Initial condition
#
# Same perturbation strategy as the successful Fig.3 code:
# - divergence-free roll-shaped initial velocity field
# - small deterministic temperature perturbation
#
# Target: Gamma = 2, beta = 0 deg, double-roll branch
# ============================================================

temperature_amplitude = 0.001
velocity_amplitude = 0.04

xmin_ic = -1.0
xmax_ic = 1.0
ymin_ic = 0.0
ymax_ic = 1.0

W_ic = xmax_ic - xmin_ic
H_ic = ymax_ic - ymin_ic

nroll_ic = 2


function stream_f_ic(x)
    ξ = (x - xmin_ic) / W_ic
    sinpi(ξ)^2 * sin(nroll_ic * pi * ξ)
end


function dstream_f_dx_ic(x)
    ξ = (x - xmin_ic) / W_ic

    term1 =
        pi *
        sin(2pi * ξ) *
        sin(nroll_ic * pi * ξ)

    term2 =
        nroll_ic *
        pi *
        sinpi(ξ)^2 *
        cos(nroll_ic * pi * ξ)

    (term1 + term2) / W_ic
end


function stream_g_ic(y)
    η = (y - ymin_ic) / H_ic
    sinpi(η)^2
end


function dstream_g_dy_ic(y)
    η = (y - ymin_ic) / H_ic
    (pi / H_ic) * sin(2pi * η)
end


start = (;

    # Divergence-free initial velocity field
    u = velocityfield(
        setup,

        (dim, x, y) -> begin

            if dim == 1
                velocity_amplitude *
                stream_f_ic(x) *
                dstream_g_dy_ic(y)
            else
                -velocity_amplitude *
                dstream_f_dx_ic(x) *
                stream_g_ic(y)
            end
        end;

        psolver,
    ),

    # Conduction state + small deterministic perturbation
    temp = temperaturefield(
        setup,

        (x, y) -> begin

            ξ = (x - xmin_ic) / W_ic
            η = (y - ymin_ic) / H_ic

            perturbation =
                (
                    cos(2pi * ξ) +
                    0.37 * cos(5pi * ξ)
                ) *
                sinpi(η)

            1.0 - η +
            temperature_amplitude *
            perturbation
        end,
    ),
)


# ============================================================
# 10. Output paths
# ============================================================

output_dir = @__DIR__

video_path = joinpath(
    output_dir,
    "Wang_beta0deg_finegrid.mp4",
)

nusselt_path = joinpath(
    output_dir,
    "Wang_beta0deg_Nusselt.png",
)

stable_nusselt_path = joinpath(
    output_dir,
    "Wang_beta0deg_Nusselt_stable.png",
)


# ============================================================
# 11. Simulation time
# ============================================================

t_end = 100.0


# ============================================================
# 12. Stable-Nu averaging interval
#
# Average Nu over t = 70 ... 100.
#
# If later you decide that the statistically stationary state
# starts earlier/later, only change this value.
# ============================================================

t_avg_start = 80.0


# ============================================================
# 13. Wang reference value
# ============================================================

Nu_Wang = 12.13


# ============================================================
# 14. Solve
# ============================================================

state, outputs = solve_unsteady(;

    force! = boussinesq_tilt!,

    setup,
    start,

    tlims = (
        0.0,
        t_end,
    ),

    force_cache = nothing,

    params = (;

        viscosity = ν,
        conductivity = κ,

        gravity = 1.0,

        beta = beta,

        # Wang temperature equation has no viscous dissipation
        dodissipation = false,
    ),

    processors = (;

        # ----------------------------------------------------
        # Temperature + velocity + mesh animation
        # ----------------------------------------------------

        anim = animator(;

            setup,

            path = video_path,

            plot = fig2plot,

            beta = beta,

            arrow_count = (16, 10),

            mesh_stride = (4, 4),

            nupdate = 20,

            framerate = 24,
        ),


        # ----------------------------------------------------
        # Nusselt-number history
        # ----------------------------------------------------

        nusselt = realtimeplotter(;

            setup,

            plot = nusseltplot,

            displayfig = false,

            nupdate = 20,
        ),
    ),
)


# ============================================================
# 15. Calculate stable mean Nu
# ============================================================

stable_indices =
    findall(
        t -> t >= t_avg_start,
        nu_times,
    )


if isempty(stable_indices)

    error(
        "No Nu data were recorded after t = $t_avg_start."
    )

end


stable_bottom =
    nu_bottom_values[stable_indices]

stable_top =
    nu_top_values[stable_indices]


Nu_bottom_mean =
    mean(stable_bottom)

Nu_bottom_std =
    std(stable_bottom)


Nu_top_mean =
    mean(stable_top)

Nu_top_std =
    std(stable_top)


# ------------------------------------------------------------
# Compare bottom hot-plate mean Nu with Wang
# ------------------------------------------------------------

error_percent =
    abs(
        Nu_bottom_mean - Nu_Wang
    ) / Nu_Wang * 100


# ============================================================
# 16. Print stable Nu results
# ============================================================

println()
println("==============================================")
println("Stable Nusselt-number results")
println("==============================================")

println(
    "Averaging interval: t = ",
    t_avg_start,
    " to ",
    t_end,
)

println()

println(
    "Mean Nu_bottom = ",
    round(
        Nu_bottom_mean;
        digits = 4,
    ),
)

println(
    "Std  Nu_bottom = ",
    round(
        Nu_bottom_std;
        digits = 4,
    ),
)

println()

println(
    "Mean Nu_top    = ",
    round(
        Nu_top_mean;
        digits = 4,
    ),
)

println(
    "Std  Nu_top    = ",
    round(
        Nu_top_std;
        digits = 4,
    ),
)

println()

println(
    "Wang Nu        = ",
    Nu_Wang,
)

println(
    "Relative error = ",
    round(
        error_percent;
        digits = 2,
    ),
    " %",
)

println("==============================================")


# ============================================================
# 17. Save original Nu(t) plot
# ============================================================

save(
    nusselt_path,
    outputs.nusselt,
)


# ============================================================
# 18. Create Nu plot with stable mean
# ============================================================

fig_stable =
    Figure(
        size = (
            900,
            520,
        ),
    )


ax_stable =
    Axis(
        fig_stable[1, 1];

        title =
            "Nusselt number - β = 0°",

        xlabel =
            "t / t_f",

        ylabel =
            "Nu",
    )


# Full Nu histories
lines!(
    ax_stable,
    nu_times,
    nu_bottom_values;
    label = "Bottom hot plate",
)

lines!(
    ax_stable,
    nu_times,
    nu_top_values;
    label = "Top cold plate",
)


# Wang reference
hlines!(
    ax_stable,
    [Nu_Wang];

    linestyle = :dash,

    label =
        "Wang DRS: Nu = $(Nu_Wang)",
)


# CFD stable mean
hlines!(
    ax_stable,
    [Nu_bottom_mean];

    linestyle = :dot,

    label =
        "CFD mean Nu = $(round(Nu_bottom_mean; digits=3))",
)


# Show start of averaging interval
vlines!(
    ax_stable,
    [t_avg_start];

    linestyle = :dashdot,

    label =
        "Averaging starts at t = $(t_avg_start)",
)


axislegend(
    ax_stable;
    position = :rt,
)


save(
    stable_nusselt_path,
    fig_stable,
)


# ============================================================
# 19. Final output
# ============================================================

println()
println("Simulation finished.")

println()
println("Video:")
println(video_path)

println()
println("Original Nusselt plot:")
println(nusselt_path)

println()
println("Stable-Nu comparison plot:")
println(stable_nusselt_path)

println()
println("Finished.")