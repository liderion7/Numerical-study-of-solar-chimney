using GLMakie
using IncompressibleNavierStokes
using Statistics
using Printf


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
# 2. Physical parameters
# ============================================================

Ra = 1.0e7
Pr = 0.71

ν = sqrt(Pr / Ra)
κ = 1.0 / sqrt(Ra * Pr)

Gamma = 8.0

println("==============================================")
println("Wang Gamma = 8, n = 5 validation")
println("==============================================")
println("Ra = ", Ra)
println("Pr = ", Pr)
println("Gamma = ", Gamma)
println("ν  = ", ν)
println("κ  = ", κ)


# ============================================================
# 3. Gamma = 8 geometry
#
# Gamma = W/H = 8
#
# x = [-4, 4]
# y = [0, 1]
#
# Vertical resolution   = 128
# Horizontal resolution = 8 × 128 = 1024
# ============================================================

n = 128

xmin_ic = -4.0
xmax_ic =  4.0

ymin_ic = 0.0
ymax_ic = 1.0

W_ic = xmax_ic - xmin_ic
H_ic = ymax_ic - ymin_ic


# ============================================================
# 4. Non-uniform near-wall mesh
# ============================================================

xgrid = tanh_grid(
    xmin_ic,
    xmax_ic,
    8n,
    1.2,
)

ygrid = tanh_grid(
    ymin_ic,
    ymax_ic,
    n,
    1.2,
)


# ============================================================
# 5. Boundary conditions
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
        # Bottom T = 0.5
        # Top T = -0.5
        temp = (
            (SymmetricBC(), SymmetricBC()),
            (DirichletBC(0.5), DirichletBC(-0.5)),
        ),
    ),
)


# ============================================================
# 6. Pressure solver
# ============================================================

psolver = default_psolver(setup)


# ============================================================
# 7. Initial-condition settings
#
# All cases use n = 5
# ============================================================

nroll_ic = 5

temperature_amplitude = 0.001
velocity_amplitude = 0.04


# ============================================================
# 8. Initial-condition functions
# ============================================================

function stream_f_ic(x)

    ξ = (x - xmin_ic) / W_ic

    sinpi(ξ)^2 *
    sin(nroll_ic * pi * ξ)
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

    (pi / H_ic) *
    sin(2pi * η)
end


# ============================================================
# 9. Create initial state
#
# Same initial-condition form for every beta
# ============================================================

function make_start()

    start = (;

        # --------------------------------------------------------
        # Divergence-free 5-roll initial velocity field
        # --------------------------------------------------------

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


        # --------------------------------------------------------
        # Conduction state + small deterministic perturbation
        # --------------------------------------------------------

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

                0.5 - η +
                temperature_amplitude *
                perturbation
            end,
        ),
    )

    return start
end


# ============================================================
# 10. Wang Table I
#
# Gamma = 8
# n = 5 for ALL cases
#
# beta =  0° -> Nu = 12.60
# beta =  3° -> Nu = 12.60
# beta =  5° -> Nu = 12.78
# beta = 10° -> Nu = 13.08
# beta = 15° -> Nu = 12.57
# beta = 20° -> Nu = 12.87
# ============================================================

betas = [
    0.0,
    3.0,
    5.0,
    10.0,
    15.0,
    20.0,
]

Nu_Wang = [
    12.60,
    12.60,
    12.78,
    13.08,
    12.57,
    12.87,
]


# ============================================================
# 11. Simulation time
# ============================================================

t_end = 100.0

t_avg_start = 80.0


# ============================================================
# 12. Output arrays
# ============================================================

Ncases = length(betas)

Nu_bottom_means = zeros(Float64, Ncases)
Nu_bottom_stds = zeros(Float64, Ncases)

Nu_top_means = zeros(Float64, Ncases)
Nu_top_stds = zeros(Float64, Ncases)

error_percents = zeros(Float64, Ncases)


# ============================================================
# 13. Run all beta cases
# ============================================================

for i in eachindex(betas)

    beta = betas[i]
    Nu_ref = Nu_Wang[i]

    println()
    println("============================================================")
    println(
        "Running case $(i)/$(Ncases): ",
        "Gamma = 8, beta = $(beta)°, n = 5"
    )
    println("============================================================")


    # --------------------------------------------------------
    # Nu histories for this case
    # --------------------------------------------------------

    nu_times = Float64[]
    nu_bottom_values = Float64[]
    nu_top_values = Float64[]


    # --------------------------------------------------------
    # Nusselt recorder
    # --------------------------------------------------------

    function nusseltplot_case(state; setup)

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
            # Save Nu history
            # ----------------------------------------------------

            push!(nu_times, Float64(t))
            push!(nu_bottom_values, Float64(Nu_b))
            push!(nu_top_values, Float64(Nu_t))


            # Dummy observables required by realtimeplotter
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


        # Figure is not displayed or saved
        fig = Figure(size = (400, 250))

        ax = Axis(
            fig[1, 1];
            xlabel = "t",
            ylabel = "Nu",
        )

        lines!(ax, Nu_bottom)
        lines!(ax, Nu_top)

        fig
    end


    # --------------------------------------------------------
    # Fresh initial state for every beta
    # --------------------------------------------------------

    start = make_start()


    # --------------------------------------------------------
    # Solve
    # --------------------------------------------------------

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

            nusselt = realtimeplotter(;

                setup,

                plot = nusseltplot_case,

                displayfig = false,

                nupdate = 20,
            ),
        ),
    )


    # --------------------------------------------------------
    # Stable interval t = 80 ... 100
    # --------------------------------------------------------

    stable_indices =
        findall(
            t -> t >= t_avg_start,
            nu_times,
        )


    if isempty(stable_indices)

        error(
            "No Nu data recorded after t = $(t_avg_start) " *
            "for beta = $(beta)°."
        )
    end


    stable_bottom =
        nu_bottom_values[stable_indices]

    stable_top =
        nu_top_values[stable_indices]


    # --------------------------------------------------------
    # Mean and standard deviation
    # --------------------------------------------------------

    Nu_bottom_mean =
        mean(stable_bottom)

    Nu_bottom_std =
        std(stable_bottom)


    Nu_top_mean =
        mean(stable_top)

    Nu_top_std =
        std(stable_top)


    # --------------------------------------------------------
    # Relative error
    # --------------------------------------------------------

    error_percent =
        abs(
            Nu_bottom_mean - Nu_ref
        ) / Nu_ref * 100


    # --------------------------------------------------------
    # Save results
    # --------------------------------------------------------

    Nu_bottom_means[i] = Nu_bottom_mean
    Nu_bottom_stds[i] = Nu_bottom_std

    Nu_top_means[i] = Nu_top_mean
    Nu_top_stds[i] = Nu_top_std

    error_percents[i] = error_percent


    # --------------------------------------------------------
    # Text output for this case
    # --------------------------------------------------------

    println()
    println("==============================================")
    println("Stable Nusselt-number results")
    println("==============================================")

    println(
        "Target state: Gamma = 8, beta = ",
        beta,
        " deg, n = 5",
    )

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
        Nu_ref,
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
end


# ============================================================
# 14. Final summary
# ============================================================

println()
println()
println("============================================================================")
println("Gamma = 8, n = 5 validation summary")
println("============================================================================")

@printf(
    "%-8s %-6s %-12s %-12s %-12s %-12s %-10s %-10s\n",
    "beta",
    "n",
    "Nu_bottom",
    "Std_bottom",
    "Nu_top",
    "Std_top",
    "Wang Nu",
    "Error (%)",
)

println("----------------------------------------------------------------------------")

for i in eachindex(betas)

    @printf(
        "%-8.1f %-6d %-12.4f %-12.4f %-12.4f %-12.4f %-10.2f %-10.2f\n",
        betas[i],
        5,
        Nu_bottom_means[i],
        Nu_bottom_stds[i],
        Nu_top_means[i],
        Nu_top_stds[i],
        Nu_Wang[i],
        error_percents[i],
    )
end

println("============================================================================")

println()

println(
    "Mean relative error = ",
    round(
        mean(error_percents);
        digits = 2,
    ),
    " %",
)


# ============================================================
# 15. Final Nu-beta comparison plot
#
# ONLY ONE plot is saved.
# ============================================================

output_dir = @__DIR__

stable_nusselt_path =
    joinpath(
        output_dir,
        "Wang_Gamma8_n5_Nusselt_stable.png",
    )


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
            "Nusselt number comparison - Γ = 8, n = 5",

        xlabel =
            "Inclination angle β (°)",

        ylabel =
            "Nu",

        xticks =
            betas,

        yticks =
            0:5:20,
    )


# ============================================================
# Our CFD
# ============================================================

lines!(
    ax_stable,
    betas,
    Nu_bottom_means;
    linewidth = 2,
    label = "Our CFD",
)

scatter!(
    ax_stable,
    betas,
    Nu_bottom_means;
    markersize = 12,
)


# ============================================================
# Wang et al.
# ============================================================

lines!(
    ax_stable,
    betas,
    Nu_Wang;
    linewidth = 2,
    linestyle = :dash,
    label = "Wang et al.",
)

scatter!(
    ax_stable,
    betas,
    Nu_Wang;
    markersize = 12,
    marker = :rect,
)


# ============================================================
# Error labels
# ============================================================

for i in eachindex(betas)

    text!(
        ax_stable,
        betas[i],
        3.0;
        text =
            "$(round(error_percents[i]; digits = 2))%",
        align = (:center, :center),
        fontsize = 13,
    )
end


# ============================================================
# Axis limits
# ============================================================

xlims!(
    ax_stable,
    -1,
    21,
)

ylims!(
    ax_stable,
    0,
    20,
)


axislegend(
    ax_stable;
    position = :rb,
)


# ============================================================
# Save ONLY final comparison figure
# ============================================================

save(
    stable_nusselt_path,
    fig_stable,
)


# ============================================================
# 16. Final message
# ============================================================

println()
println("All Gamma = 8, n = 5 simulations finished.")

println()

println(
    "Stable Nu comparison figure:"
)

println(
    stable_nusselt_path
)

println()
println("Finished.")
