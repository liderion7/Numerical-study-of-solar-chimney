# plot_huang_curves.jl
# -----------------------------------------------------------------------------
# Validation curves for the Huang et al. (2022) benchmark.
#
# Scans every case under validation_huang/, reads its
# huang_validation_comparison.csv and plots
#
#     Nu, Re, Rex, Rey   versus   Ra
#
# for the CFD runs and for the published reference values on the same axes.
#
# For each Ra the record with the longest averaging window (tavg) is used as the
# representative curve point; every other record is still listed in the CSV with
# a `selected` flag so nothing is silently dropped.
#
# Usage:
#   julia plot_huang_curves.jl
#   ROOT=validation_huang julia plot_huang_curves.jl
#
# Outputs (next to the cases, i.e. inside ROOT):
#   huang_validation_curves.png
#   huang_validation_curves.csv
# -----------------------------------------------------------------------------

using GLMakie
using DelimitedFiles
using Printf

base = @__DIR__
const ROOT = get(ENV, "ROOT", joinpath(base, "data"))
isdir(ROOT) || error("case root not found: $ROOT")

const QUANTITIES = ("Nu", "Re", "Rex", "Rey")

# -----------------------------------------------------------------------------
# 1. collect every result row
# -----------------------------------------------------------------------------
records = NamedTuple[]

for d in sort(readdir(ROOT; join = true))
    isdir(d) || continue
    f = joinpath(d, "huang_validation_comparison.csv")
    isfile(f) || continue

    raw = readdlm(f, ',', '\n'; skipstart = 1)
    (size(raw, 1) >= 1 && size(raw, 2) >= 17) || continue

    for i in 1:size(raw, 1)
        r = raw[i, :]
        push!(records, (
            dir = basename(d),
            Ra = Float64(r[1]),
            Pr = Float64(r[2]),
            N = Int(Float64(r[3])),
            Ny = Int(Float64(r[4])),
            tavg = Float64(r[5]),
            Nu = Float64(r[6]),
            Nu_ref = Float64(r[7]),
            Re = Float64(r[9]),
            Re_ref = Float64(r[10]),
            Rex = Float64(r[12]),
            Rex_ref = Float64(r[13]),
            Rey = Float64(r[15]),
            Rey_ref = Float64(r[16]),
        ))
    end
end

isempty(records) && error("no huang_validation_comparison.csv found under $ROOT")

# -----------------------------------------------------------------------------
# 2. pick one representative record per Ra (longest averaging window wins)
# -----------------------------------------------------------------------------
groups = Dict{Float64, Vector{NamedTuple}}()
for r in records
    push!(get!(groups, r.Ra, NamedTuple[]), r)
end

selected = NamedTuple[]
excluded = NamedTuple[]
for ra in sort(collect(keys(groups)))
    g = groups[ra]
    # Degenerate records are dropped first: a run whose flow never started
    # (Re = 0, e.g. the u = 0 / theta = 0 initial condition) carries no
    # information and must not become the representative point of its Ra.
    keep = findall(r -> r.Re > 1e-8, g)
    isempty(keep) && (keep = collect(eachindex(g)))
    # Longest averaging window wins; ties broken by directory name.
    sel = keep[argmax([(g[i].tavg, g[i].dir) for i in keep])]
    push!(selected, g[sel])
    for i in eachindex(g)
        i == sel || push!(excluded, g[i])
    end
end

println("records excluded (degenerate or shorter window) = ", length(excluded))
for r in excluded
    @printf("  excluded: %-46s Ra=%-8.3g tavg=%.1f  Re=%.4g\n", r.dir, r.Ra, r.tavg, r.Re)
end

println("records found    = ", length(records))
println("distinct Ra      = ", length(selected))
for r in selected
    @printf("  Ra = %-8.3g  representative dir = %-46s tavg = %.1f  N = %d\n",
            r.Ra, r.dir, r.tavg, r.N)
end

# -----------------------------------------------------------------------------
# 3. summary CSV (all records, with the selection flag)
# -----------------------------------------------------------------------------
selected_keys = Set((r.dir, r.Ra) for r in selected)

csv = joinpath(ROOT, "huang_validation_curves.csv")
open(csv, "w") do io
    println(io, "Ra,N_wall,Ny,tavg,",
                "Nu_cfd,Nu_ref,Re_cfd,Re_ref,Rex_cfd,Rex_ref,Rey_cfd,Rey_ref,selected,source_dir")
    for r in sort(records; by = r -> (r.Ra, -r.tavg, r.dir))
        println(
            io, r.Ra, ',', r.N, ',', r.Ny, ',', r.tavg, ',',
            r.Nu, ',', r.Nu_ref, ',',
            r.Re, ',', r.Re_ref, ',',
            r.Rex, ',', r.Rex_ref, ',',
            r.Rey, ',', r.Rey_ref, ',',
            (r.dir, r.Ra) in selected_keys,
            ',', r.dir,
        )
    end
end
println("wrote ", csv)

# -----------------------------------------------------------------------------
# 4. figure: one panel per quantity
# -----------------------------------------------------------------------------
ras = [r.Ra for r in selected]
nda = length(unique(r.N for r in selected)) == 1 ? selected[1].N : 0

fig = Figure(size = (1400, 950))

for (k, q) in enumerate(QUANTITIES)
    row = (k - 1) ÷ 2 + 1
    col = (k - 1) % 2 + 1

    ax = Axis(
        fig[row, col];
        title = q,
        xlabel = "Ra",
        ylabel = q,
        xscale = log10,
    )

    cfd = [getproperty(r, Symbol(q)) for r in selected]
    ref = [getproperty(r, Symbol(q * "_ref")) for r in selected]

    scatterlines!(
        ax, ras, cfd;
        marker = :circle,
        label = nda > 0 ? "CFD (N = $nda)" : "CFD",
    )
    lines!(ax, ras, ref; linestyle = :dash, linewidth = 2, label = "Huang (2022)")
    scatter!(ax, ras, ref; marker = :rect)

    if k == 1
        axislegend(ax; position = :lt)
    end
end

save(joinpath(ROOT, "huang_validation_curves.png"), fig)
println("wrote ", joinpath(ROOT, "huang_validation_curves.png"))

# -----------------------------------------------------------------------------
# 5. compact error table on stdout
# -----------------------------------------------------------------------------
println()
@printf("%-8s %10s %10s %8s %10s %10s %8s\n", "Ra", "Nu", "Nu_ref", "err%", "Re", "Re_ref", "err%")
for r in selected
    @printf("%-8.3g %10.4f %10.4f %8.2f %10.4f %10.4f %8.2f\n",
            r.Ra, r.Nu, r.Nu_ref, 100 * (r.Nu - r.Nu_ref) / r.Nu_ref,
            r.Re, r.Re_ref, 100 * (r.Re - r.Re_ref) / r.Re_ref)
end
println()
@printf("%-8s %10s %10s %8s %10s %10s %8s\n", "Ra", "Rex", "Rex_ref", "err%", "Rey", "Rey_ref", "err%")
for r in selected
    @printf("%-8.3g %10.4f %10.4f %8.2f %10.4f %10.4f %8.2f\n",
            r.Ra, r.Rex, r.Rex_ref, 100 * (r.Rex - r.Rex_ref) / r.Rex_ref,
            r.Rey, r.Rey_ref, 100 * (r.Rey - r.Rey_ref) / r.Rey_ref)
end

println("done.")
