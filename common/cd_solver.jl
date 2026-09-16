using IncompressibleNavierStokes

"""
    boussinesq_cd!(force, state, t; ...)

Task C/D Boussinesq forcing under the project coordinate convention:

- x is wall-normal (hot -> cold),
- y is vertical / streamwise,
- theta is measured from +x toward +y,
- theta = 90° therefore makes the chimney axis parallel to +y.

The buoyancy direction is therefore

    e(theta) = (cos(theta), sin(theta)).

`drive_force_y` is an optional uniform streamwise forcing, representing the
additional pressure-gradient / external driving force discussed in the
meetings. Leave it at 0.0 until a physical value is specified.

IMPORTANT: There is deliberately NO zero-mass-flux constraint in this Task C/D
solver. The purpose of Task D is to measure the resulting bulk velocity, so a
zero-flux constraint would prescribe the answer to be zero.
"""
function boussinesq_cd!(
    force,
    state,
    t;
    setup,
    cache,
    viscosity,
    conductivity,
    gravity,
    theta,
    drive_force_y,
    dodissipation,
)
    (; u, temp) = state

    fill!(force.u, 0)
    fill!(force.temp, 0)

    # Momentum convection + viscous diffusion.
    IncompressibleNavierStokes.convectiondiffusion!(
        force.u,
        u,
        setup,
        viscosity,
    )

    # Coordinate convention: theta = 90 deg -> buoyancy parallel to +y.
    angle = deg2rad(theta)
    gx = gravity * cos(angle)
    gy = gravity * sin(angle)

    IncompressibleNavierStokes.applygravity!(
        force.u,
        temp,
        setup,
        1,
        gx,
    )
    IncompressibleNavierStokes.applygravity!(
        force.u,
        temp,
        setup,
        2,
        gy,
    )

    # Optional uniform vertical / streamwise driving force.
    # This is equivalent to a constant pressure-gradient body force.
    if !iszero(drive_force_y)
        @views force.u[:, :, 2] .+= drive_force_y
    end

    # Temperature convection + diffusion.
    IncompressibleNavierStokes.convection_diffusion_temp!(
        force.temp,
        u,
        temp,
        setup,
        conductivity,
    )

    if dodissipation
        IncompressibleNavierStokes.dissipation!(
            force.temp,
            u,
            setup,
            viscosity,
        )
    end

    nothing
end

IncompressibleNavierStokes.get_cache(::typeof(boussinesq_cd!), setup) = nothing

IncompressibleNavierStokes.propose_timestep(
    ::typeof(boussinesq_cd!),
    state,
    setup,
    params,
) = IncompressibleNavierStokes.propose_timestep(
    IncompressibleNavierStokes.boussinesq!,
    state,
    setup,
    params,
)

"""
    cd_dimensionless_coefficients(Ra, Pr)

Free-fall nondimensional coefficients used by the existing project code.
"""
function cd_dimensionless_coefficients(Ra::Real, Pr::Real)
    Ra > 0 || throw(ArgumentError("Ra must be positive"))
    Pr > 0 || throw(ArgumentError("Pr must be positive"))

    viscosity = sqrt(Pr / Ra)
    conductivity = 1 / sqrt(Ra * Pr)

    (; viscosity, conductivity)
end
