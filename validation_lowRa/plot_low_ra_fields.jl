# plot_low_ra_fields.jl
# -----------------------------------------------------------------------------
# Low-Ra visualisation for EVERY angle that has field snapshots.
#
# Panel convention for these figures (as requested for the report):
#     horizontal axis = x   (wall-normal: hot wall at x = 0, cold wall at x = 1)
#     vertical   axis = z   (chimney axis / streamwise direction, periodic)
#
# So code `x` is drawn horizontally and code `y` (the periodic direction) is
# drawn vertically and labelled `z`.
#
# For every <case>/field_t<TIME>.csv under ROOT (produced by run_cd_case.jl with
# SNAPSHOT_TIMES set) this writes three-panel figures to
#
#     <case>/validation_2d/<case>_<TAG>.png
#
#     (1) temperature field, blue-red  (-0.5 cold ... +0.5 hot)
#     (2) velocity magnitude |u|       (common colour scale across all angles)
#     (3) streamlines (streamfunction contours) over the temperature field
#
# Zero-flow cases (theta = 0 deg, where the buoyancy acts along the wall normal
# so no convection starts) are labelled explicitly instead of showing the
# amplification of machine noise.
#
# An overview figure with the last snapshot of every angle is written to
#
#     <ROOT>/lowRa_fields_overview.png
#
# Usage:
#   julia plot_low_ra_fields.jl
#   ROOT=data/lowRa_fields julia plot_low_ra_fields.jl
#
# Optional environment variables:
#   ROOT        : results root to scan          (default data)
#   PATTERN     : snapshot glob                 (default field_t*.csv)
#   OUTDIRNAME  : sub-directory for the figures (default validation_2d)
#   FLOW_TOL    : below this max|u| a case is reported as "no flow" (default 1e-6)
# -----------------------------------------------------------------------------

using GLMakie
using Printf

base = @__DIR__
const ROOT = get(ENV, "ROOT", joinpath(base, "data"))
const PATTERN = get(ENV, "PATTERN", "field_t*.csv")
const OUTDIRNAME = get(ENV, "OUTDIRNAME", "validation_2d")
const FLOW_TOL = parse(Float64, get(ENV, "FLOW_TOL", "1e-6"))

isdir(ROOT) || error("results root not found: $ROOT")

# -----------------------------------------------------------------------------
# Data reading
# -----------------------------------------------------------------------------
"""
    read_field_snapshot(path)

Return `(t, x, y, T, ux, uy)` from one snapshot CSV. `T`, `ux` and `uy` are
matrices indexed `[i, j]` with `i` along code `x` (wall-normal) and `j` along
code `y` (the periodic / vertical direction).
"""
function read_field_snapshot(path)
    t = NaN
    rows = Vector{NTuple{7, Float64}}()

    open(path) do io
        for line in eachline(io)
            if startswith(line, "#")
                m = match(r"t\s*=\s*([0-9eE\+\-\.]+)", line)
                m !== nothing && (t = parse(Float64, m.captures[1]))
                continue
            end
            startswith(line, "ix") && continue
            isempty(strip(line)) && continue
            parts = split(line, ',')
            length(parts) == 7 || continue
            push!(rows, Tuple(parse(Float64, p) for p in parts))
        end
    end

    isempty(rows) && error("no data rows found in $path")

    nx = maximum(r -> Int(r[1]), rows)
    ny = maximum(r -> Int(r[2]), rows)

    X = zeros(nx, ny); Y = zeros(nx, ny); T = zeros(nx, ny)
    U = zeros(nx, ny)   # code u_x : wall-normal velocity
    V = zeros(nx, ny)   # code u_y : streamwise velocity (the "z" component)

    for r in rows
        i = Int(r[1]); j = Int(r[2])
        X[i, j] = r[3]; Y[i, j] = r[4]; T[i, j] = r[5]
        U[i, j] = r[6]; V[i, j] = r[7]
    end

    t, X[:, 1], Y[1, :], T, U, V
end

"""
    streamfunction(x, uy)

Integrate `u_y` along the wall-normal direction to obtain psi with
`u_y = d psi / d x` and `u_x = -d psi / d y` (incompressibility). The wall-normal
velocity vanishes on both walls, so psi is constant there.
"""
function streamfunction(x::AbstractVector, uy::AbstractMatrix)
    psi = zeros(size(uy))
    for j in axes(uy, 2)
        acc = 0.0
        for i in 2:size(uy, 1)
            acc += 0.5 * (uy[i, j] + uy[i - 1, j]) * (x[i] - x[i - 1])
            psi[i, j] = acc
        end
    end
    psi
end

function glob_to_regex(glob::AbstractString)
    out = IOBuffer()
    print(out, '^')
    for c in glob
        if c == '*'
            print(out, ".*")
        elseif c in ('.', '+', '(', ')', '[', ']', '{', '}', '^', '$', '|', '\\')
            print(out, '\\', c)
        else
            print(out, c)
        end
    end
    print(out, '$')
    Regex(String(take!(out)))
end

# -----------------------------------------------------------------------------
# Scan and load everything once (so the colour scale can be shared)
# -----------------------------------------------------------------------------
rx = glob_to_regex(PATTERN)

snaps = NamedTuple[]
for d in sort(readdir(ROOT; join = true))
    isdir(d) || continue
    case = basename(d)
    for f in sort([f for f in readdir(d) if occursin(rx, f)])
        t, x, y, T, U, V = read_field_snapshot(joinpath(d, f))
        push!(snaps, (
            case = case,
            file = f,
            path = joinpath(d, f),
            t = t,
            x = x,
            y = y,
            T = T,
            U = U,
            V = V,
            speed = sqrt.(U .^ 2 .+ V .^ 2),
            psi = streamfunction(x, V),
        ))
    end
end

isempty(snaps) && error("no snapshots matching $PATTERN found under $ROOT")

global_speed = maximum(s -> maximum(s.speed), snaps)
flowing_any = global_speed > FLOW_TOL
speed_max = flowing_any ? global_speed : FLOW_TOL

println("root             = ", ROOT)
println("snapshots loaded = ", length(snaps), " in ",
        length(unique(s -> s.case, snaps)), " cases")
@printf("global max |u|   = %.6e  (flow threshold %.1e)\n", global_speed, FLOW_TOL)

# -----------------------------------------------------------------------------
# Per-snapshot figure
# -----------------------------------------------------------------------------
function plot_snapshot(s, speed_max)
    outdir = joinpath(dirname(s.path), OUTDIRNAME)
    mkpath(outdir)

    local_max = maximum(s.speed)
    flowing = local_max > FLOW_TOL
    tag = @sprintf("%s  (t = %.1f)", s.case, s.t)

    fig = Figure(size = (1560, 640))

    # (1) temperature, blue-red
    ax1 = Axis(fig[1, 1]; title = "temperature T — $tag",
               xlabel = "x", ylabel = "z", aspect = DataAspect())
    hm1 = heatmap!(ax1, s.x, s.y, s.T; colormap = :balance, colorrange = (-0.5, 0.5))
    Colorbar(fig[1, 2], hm1; label = "T")

    # (2) velocity magnitude, common colour scale
    ax2 = Axis(fig[1, 3]; title = "|u| — $tag",
               xlabel = "x", ylabel = "z", aspect = DataAspect())
    hm2 = heatmap!(ax2, s.x, s.y, s.speed; colormap = :viridis,
                   colorrange = (0.0, speed_max))
    Colorbar(fig[1, 4], hm2; label = "|u|")

    # (3) temperature with streamlines
    ax3 = Axis(fig[1, 5]; title = "streamlines over T — $tag",
               xlabel = "x", ylabel = "z", aspect = DataAspect())
    heatmap!(ax3, s.x, s.y, s.T; colormap = :balance, colorrange = (-0.5, 0.5))

    xmid = 0.5 * (s.x[1] + s.x[end])
    ymid = 0.5 * (s.y[1] + s.y[end])
    if flowing
        contour!(ax3, s.x, s.y, s.psi; levels = 13, color = (:black, 0.55), linewidth = 1.0)
    else
        text!(ax3, xmid, ymid;
              text = @sprintf("no flow\nmax |u| = %.2e\n(pure conduction)", local_max),
              align = (:center, :center), fontsize = 16)
    end

    if !flowing
        text!(ax2, xmid, ymid;
              text = @sprintf("no flow\nmax |u| = %.2e", local_max),
              align = (:center, :center), fontsize = 14, color = :white)
    end

    png = joinpath(outdir, s.case * "_" * replace(basename(s.file), ".csv" => "") * ".png")
    save(png, fig)
    png
end

# -----------------------------------------------------------------------------
# Overview: last snapshot of every angle
# -----------------------------------------------------------------------------
function plot_overview(snaps, speed_max)
    # last snapshot per case, ordered by case name
    by_case = Dict{String, NamedTuple}()
    for s in snaps
        if !haskey(by_case, s.case) || s.t > by_case[s.case].t
            by_case[s.case] = s
        end
    end
    names = sort(collect(keys(by_case)))
    isempty(names) && return nothing

    ncol = min(length(names), 5)
    nrow = cld(length(names), ncol)

    fig = Figure(size = (320 * ncol + 140, 460 * nrow))
    for (k, name) in enumerate(names)
        s = by_case[name]
        r = (k - 1) ÷ ncol + 1
        c = (k - 1) % ncol + 1

        ax = Axis(fig[r, c];
                  title = @sprintf("%s  (t = %.0f)", name, s.t),
                  xlabel = "x", ylabel = "z", aspect = DataAspect())
        heatmap!(ax, s.x, s.y, s.speed; colormap = :viridis, colorrange = (0.0, speed_max))
        if maximum(s.speed) > FLOW_TOL
            contour!(ax, s.x, s.y, s.psi; levels = 11,
                     color = (:white, 0.6), linewidth = 0.8)
        else
            text!(ax, 0.5 * (s.x[1] + s.x[end]), 0.5 * (s.y[1] + s.y[end]);
                  text = "no flow", align = (:center, :center),
                  fontsize = 13, color = :white)
        end
    end

    png = joinpath(ROOT, "lowRa_fields_overview.png")
    save(png, fig)
    png
end

# -----------------------------------------------------------------------------
# Write everything
# -----------------------------------------------------------------------------
for s in snaps
    println("wrote ", plot_snapshot(s, speed_max))
end

ov = plot_overview(snaps, speed_max)
ov === nothing || println("overview         = ", ov)
println("done.")
