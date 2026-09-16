# plot_huang_fields.jl
# -----------------------------------------------------------------------------
# Plot the full-field snapshots written by make_huang_field_recorder
# (see common/cd_validation.jl).
#
# For every `field_t*.csv` in the case directory this writes a three-panel
# figure so that the convection-roll pattern can be read off directly:
#
#   (1) temperature T(x, y) with contours
#   (2) velocity magnitude |u|(x, y)
#   (3) streamlines, drawn as iso-contours of the streamfunction psi
#
# psi is obtained by integrating u_x along the wall-normal direction, so that
# u_x = d psi / d y and u_y = -d psi / d x. Makie's own `streamplot` cannot be
# used here: the wall-normal grid is tanh-clustered, and streamplot rejects the
# resulting non-uniform coordinate vectors.
#
# Usage:
#   CASE_DIR=validation_huang/Huang_NS_PD_Ra1e06_Pr4p3_N64_perturbed_W1 julia plot_huang_fields.jl
#
# The horizontal axis is the periodic direction (code y), the vertical axis is
# the wall-normal direction (code x), so the 2:1 domain aspect ratio shows up
# horizontally, as in the paper.
#
# Optional environment variables:
#   OUTDIR     : where the PNGs go (default <case>/figures)
#   PATTERN    : glob for snapshots (default "field_t*.csv")
#   PSI_LEVELS : number of streamline levels (default 13)
# -----------------------------------------------------------------------------

using GLMakie
using Printf

base = @__DIR__

const CASE_DIR = get(ENV, "CASE_DIR", joinpath(base, "data"))
isdir(CASE_DIR) || error("case directory not found: $CASE_DIR (set CASE_DIR)")

const OUTDIR = get(ENV, "OUTDIR", joinpath(CASE_DIR, "figures"))
const PATTERN = get(ENV, "PATTERN", "field_t*.csv")
const PSI_LEVELS = parse(Int, get(ENV, "PSI_LEVELS", "13"))

"""
    read_field_snapshot(path)

Return `(t, x, y, T, ux, uy)` from one snapshot CSV, where `T`, `ux` and `uy`
are matrices indexed `[i, j]`: `i` runs along the wall-normal direction (code x)
and `j` along the periodic direction (code y).
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

    X = zeros(nx, ny)
    Y = zeros(nx, ny)
    T = zeros(nx, ny)
    U = zeros(nx, ny)   # code u_x = paper vertical velocity
    V = zeros(nx, ny)   # code u_y = paper horizontal velocity

    for r in rows
        i = Int(r[1])
        j = Int(r[2])
        X[i, j] = r[3]
        Y[i, j] = r[4]
        T[i, j] = r[5]
        U[i, j] = r[6]
        V[i, j] = r[7]
    end

    x = X[:, 1]
    y = Y[1, :]
    t, x, y, T, U, V
end

"""
    streamfunction(x, ux)

Integrate `u_x` along the wall-normal direction to obtain psi(x, y) with
`u_x = d psi / d y` and `u_y = -d psi / d x`. The wall-normal velocity vanishes
at both plates, so psi is constant along them (taken as zero here).
"""
function streamfunction(x::AbstractVector, ux::AbstractMatrix)
    psi = zeros(size(ux))
    for j in axes(ux, 2)
        acc = 0.0
        for i in 2:size(ux, 1)
            acc += 0.5 * (ux[i, j] + ux[i - 1, j]) * (x[i] - x[i - 1])
            psi[i, j] = acc
        end
    end
    psi
end

"""Turn a simple glob such as `field_t*.csv` into an anchored Regex."""
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

"""
    plot_field_snapshot(path)

Write the three-panel figure for one snapshot and return the PNG path.
"""
function plot_field_snapshot(path)
    t, x, y, T, U, V = read_field_snapshot(path)

    # Matrices handed to Makie must be (length(y_horizontal), length(x_vertical)).
    Tp = permutedims(T)
    speed = permutedims(sqrt.(U .^ 2 .+ V .^ 2))
    psip = permutedims(streamfunction(x, U))

    tag = replace(basename(path), ".csv" => "")
    fig = Figure(size = (1750, 540))

    ax1 = Axis(
        fig[1, 1];
        title = @sprintf("T  (t = %.1f)", t),
        xlabel = "periodic direction (code y)",
        ylabel = "wall-normal (code x)",
        aspect = DataAspect(),
    )
    hm1 = heatmap!(ax1, y, x, Tp; colormap = :balance, colorrange = (-0.5, 0.5))
    contour!(ax1, y, x, Tp; levels = 11, color = (:black, 0.45), linewidth = 0.8)
    Colorbar(fig[1, 2], hm1; label = "T")

    ax2 = Axis(
        fig[1, 3];
        title = "|u|",
        xlabel = "periodic direction (code y)",
        ylabel = "wall-normal (code x)",
        aspect = DataAspect(),
    )
    hm2 = heatmap!(ax2, y, x, speed; colormap = :viridis)
    Colorbar(fig[1, 4], hm2; label = "|u|")

    ax3 = Axis(
        fig[1, 5];
        title = "streamlines (psi contours)",
        xlabel = "periodic direction (code y)",
        ylabel = "wall-normal (code x)",
        aspect = DataAspect(),
    )
    contour!(ax3, y, x, psip; levels = PSI_LEVELS, linewidth = 1.2)

    png = joinpath(OUTDIR, tag * ".png")
    mkpath(OUTDIR)
    save(png, fig)
    png
end

rx = glob_to_regex(PATTERN)
files = sort([f for f in readdir(CASE_DIR) if occursin(rx, f)])
isempty(files) && error("no snapshots matching $PATTERN in $CASE_DIR")

println("case directory  = ", CASE_DIR)
println("output directory= ", OUTDIR)
println("snapshots found = ", length(files))

for f in files
    png = plot_field_snapshot(joinpath(CASE_DIR, f))
    println("wrote ", png)
end

println("done.")
