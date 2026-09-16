# make_low_ra_video.jl
# -----------------------------------------------------------------------------
# Build an animation from a dense set of field snapshots.
#
# Reads every `field_t*.csv` in CASE_DIR (written by run_cd_case.jl with
# SNAPSHOT_TIMES covering many times), keeps them in memory and records a
# two-panel movie:
#
#     left  panel : temperature T, blue-red  (-0.5 cold ... +0.5 hot)
#     right panel : velocity magnitude |u| with streamfunction contours
#
# Axes follow the report convention used elsewhere in this project:
#     horizontal axis = x  (wall-normal, hot wall x = 0 -> cold wall x = 1)
#     vertical   axis = z  (chimney axis / streamwise direction, periodic)
#
# The output format is taken from the file extension: use `.mp4` (needs ffmpeg,
# shipped with Makie) or `.gif`. The file is written next to the snapshots.
#
# Usage:
#   CASE_DIR=data/lowRa_movie/90degree julia make_low_ra_video.jl
#
# Optional environment variables:
#   CASE_DIR  : directory holding the snapshots (required unless the default exists)
#   OUT       : output file name (default "lowRa_<case>_T_u.mp4" inside CASE_DIR)
#   FPS       : frames per second (default 12)
#   STRIDE    : use every STRIDE-th snapshot (default 1)
#   RES       : figure resolution, e.g. "1400x700" (default 1400x700)
# -----------------------------------------------------------------------------

using GLMakie
using Printf

base = @__DIR__
const CASE_DIR = get(ENV, "CASE_DIR", joinpath(base, "data", "lowRa_movie", "90degree"))
const FPS = parse(Int, get(ENV, "FPS", "12"))
const STRIDE = parse(Int, get(ENV, "STRIDE", "1"))
const RES = let s = split(get(ENV, "RES", "1400x700"), 'x')
    (parse(Int, s[1]), parse(Int, s[2]))
end

isdir(CASE_DIR) || error("case directory not found: $CASE_DIR (set CASE_DIR)")

const OUT = get(ENV, "OUT", joinpath(CASE_DIR, "lowRa_$(basename(CASE_DIR))_T_u.mp4"))

# -----------------------------------------------------------------------------
# Read one snapshot (same format as the other plotting scripts)
# -----------------------------------------------------------------------------
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
    U = zeros(nx, ny); V = zeros(nx, ny)

    for r in rows
        i = Int(r[1]); j = Int(r[2])
        X[i, j] = r[3]; Y[i, j] = r[4]; T[i, j] = r[5]
        U[i, j] = r[6]; V[i, j] = r[7]
    end

    t, X[:, 1], Y[1, :], T, U, V
end

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

# -----------------------------------------------------------------------------
# Collect the frames in time order
# -----------------------------------------------------------------------------
files = [f for f in readdir(CASE_DIR) if occursin(r"^field_t.*\.csv$", f)]
isempty(files) && error("no field_t*.csv snapshots in $CASE_DIR")

times = [parse(Float64, match(r"field_t([0-9eE\+\-\.]+)\.csv", f).captures[1]) for f in files]
order = sortperm(times)
files = files[order]
times = times[order]

if STRIDE > 1
    keep = 1:STRIDE:length(files)
    files = files[keep]
    times = times[keep]
end

println("case directory = ", CASE_DIR)
println("frames         = ", length(files), "  (t = ", times[1], " ... ", times[end], ")")
println("output         = ", OUT)
println("fps            = ", FPS, "   resolution = ", RES)

frames = [read_field_snapshot(joinpath(CASE_DIR, f)) for f in files]

_, x, y, T1, U1, V1 = frames[1]
speeds = [sqrt.(f[5] .^ 2 .+ f[6] .^ 2) for f in frames]
psis = [streamfunction(x, f[6]) for f in frames]
speed_max = maximum(maximum(s) for s in speeds)
flowing = speed_max > 1e-6

@printf("global max |u| = %.6e%s\n", speed_max, flowing ? "" : "   (no flow)")

# -----------------------------------------------------------------------------
# Record
# -----------------------------------------------------------------------------
T_obs = Observable(T1)
speed_obs = Observable(speeds[1])
psi_obs = Observable(psis[1])

fig = Figure(size = RES)
ax1 = Axis(fig[1, 1]; title = "temperature T", xlabel = "x", ylabel = "z",
           aspect = DataAspect())
hm1 = heatmap!(ax1, x, y, T_obs; colormap = :balance, colorrange = (-0.5, 0.5))
Colorbar(fig[1, 2], hm1; label = "T")

ax2 = Axis(fig[1, 3]; title = "|u| with streamlines", xlabel = "x", ylabel = "z",
           aspect = DataAspect())
hm2 = heatmap!(ax2, x, y, speed_obs; colormap = :viridis, colorrange = (0.0, speed_max))
Colorbar(fig[1, 4], hm2; label = "|u|")

# Streamfunction contours. Updated in place when the observable changes; on
# Makie versions that cannot update a contour recipe in place the call simply
# rebuilds the contour plot, which is fine for the frame counts used here.
cont = Observable{Any}(nothing)
if flowing
    contour!(ax2, x, y, psi_obs; levels = 13, color = (:white, 0.6), linewidth = 0.9)
end

title_obs = Observable(@sprintf("t = %.1f", times[1]))
Label(fig[0, 1:4], title_obs; tellwidth = false, fontsize = 18)

n = length(frames)
record(fig, OUT, 1:n; framerate = FPS) do i
    T_obs[] = frames[i][4]
    speed_obs[] = speeds[i]
    flowing && (psi_obs[] = psis[i])
    title_obs[] = @sprintf("Ra = 1000, theta = 90 deg,   t = %.1f", times[i])
end

println("wrote ", OUT)
