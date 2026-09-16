"""
    build_cd_setup(n = 128; gap = 1.0, periodic_length = 2.0,
                   clustering = 1.2, backend = nothing)

2-D Task C/D geometry using the coordinate convention agreed for the project:

- `x`: wall-normal coordinate, from hot wall to cold wall.
- `y`: vertical / streamwise coordinate.
- `theta = 90°`: chimney axis is parallel to +y.

The x walls are no-slip and isothermal:
    x = 0       : T = +0.5 (hot)
    x = gap     : T = -0.5 (cold)
The y direction is periodic.

The x mesh is tanh-clustered near the hot/cold walls. The periodic y mesh is
uniform. With this convention, the conduction temperature is a one-dimensional
function of x,

    T_0(x) = 0.5 - x/gap,

which is exactly the form required for the low-Rayleigh-number validation.
"""
function build_cd_setup(
    n::Integer = 128;
    gap::Real = 1.0,
    periodic_length::Real = 2.0,
    clustering::Real = 1.2,
    ny::Union{Nothing, Integer} = nothing,
    backend = nothing,
)
    n > 0 || throw(ArgumentError("n must be positive"))
    gap > 0 || throw(ArgumentError("gap must be positive"))
    periodic_length > 0 || throw(ArgumentError("periodic_length must be positive"))
    clustering > 0 || throw(ArgumentError("clustering must be positive"))

    # Preserve roughly the same physical resolution in x and y.
    ny_cells = isnothing(ny) ? max(8, round(Int, n * periodic_length / gap)) : Int(ny)
    ny_cells > 0 || throw(ArgumentError("ny must be positive"))

    # Wall-normal x: clustered near both physical walls.
    xgrid = tanh_grid(0.0, Float64(gap), n, clustering)

    # Vertical / streamwise y: periodic and uniform.
    ygrid = range(0.0, Float64(periodic_length); length = ny_cells + 1)

    setup_kwargs = isnothing(backend) ? NamedTuple() : (; backend)

    Setup(;
        x = (xgrid, ygrid),
        boundary_conditions = (;
            u = (
                # x = 0, gap: physical no-slip walls
                (DirichletBC(), DirichletBC()),
                # y: periodic streamwise direction
                (PeriodicBC(), PeriodicBC()),
            ),
            temp = (
                # Hot wall -> cold wall, so T decreases with x.
                (DirichletBC(0.5), DirichletBC(-0.5)),
                # Periodic in y.
                (PeriodicBC(), PeriodicBC()),
            ),
        ),
        setup_kwargs...,
    )
end

"""
    build_cd_initial_state(setup, psolver; gap = 1.0,
                           periodic_length = 2.0, ϵ = 1e-3)

Stationary initial velocity plus a conduction temperature field with a small
periodic perturbation. The perturbation is zero at both x walls, so it does not
alter the prescribed wall temperatures.
"""
function build_cd_initial_state(
    setup,
    psolver;
    gap::Real = 1.0,
    periodic_length::Real = 2.0,
    ϵ::Real = 1e-3,
)
    (
        u = velocityfield(
            setup,
            (dim, x, y) -> zero(x);
            psolver,
        ),
        temp = temperaturefield(
            setup,
            (x, y) ->
                0.5 - x / gap +
                ϵ * sinpi(x / gap) * sin(2π * y / periodic_length),
        ),
    )
end
