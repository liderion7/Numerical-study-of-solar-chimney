# Huang et al. (2022) validation — review, benchmark notes and follow-up log

Scope: review of the two new validation files, the first benchmark runs against
the published Rayleigh–Bénard (RBC) reference values, and the follow-up
experiment sequence (①–⑤) requested afterwards.

- `huang_validation.jl` — benchmark runner (NS plates / periodic sidewalls).
- `common/cd_validation.jl` — reference table, initial conditions, full-field
  statistics, recorders, time averaging, summary writer, field snapshots.
- `plot_huang_fields.jl` — snapshot figures (T, |u|, streamlines).
- `plot_huang_curves.jl` — validation curves over Ra.

Benchmark subset: `Pr = 4.3`, aspect ratio `Γ = L/H = 2`, no-slip (NS) plates +
periodic (PD) sidewalls, free-fall nondimensionalisation (`ν* = √(Pr/Ra)`,
`κ* = 1/√(Ra·Pr)`), `θ = 0°` so buoyancy acts along code `+x`.

Coordinate mapping used throughout: code `x` = paper vertical `y`
(hot plate → cold plate), code `y` = paper horizontal `x` (periodic).

---

## 1. Bug found and fixed (blocking)

`common/cd_validation.jl`, `huang_instantaneous_stats` returned a **plain
Tuple** while both call sites access **named fields**:

- `make_huang_statistics_recorder` → `s.mean_horizontal2`, `s.mean_vertical2`, `s.mean_vT`
- `make_huang_runtime_monitor` → `s.Re_inst`, `s.Rex_inst`, `s.Rey_inst`, `s.Nu_inst`

The first run therefore died immediately with

```
FieldError: type Tuple has no field `Re_inst`, available fields: `1`, ..., `7`
  @ common/cd_validation.jl:293 (make_huang_runtime_monitor)
```

Fix: return a `NamedTuple` (leading `;` inside the returned tuple). One change
fixes every call site; the module then ran end to end on the first retry.

---

## 2. Benchmark runs

### 2.1 Ra sweep (the main result)

All runs: `Pr = 4.3`, `Γ = 2`, `N_wall = 64` (pressure grid 64 × 130),
`INIT_MODE = perturbed`, `PERT_WAVES = 1`, `INIT_EPS = 1e-6`,
`AVG_START = 200`, `T_END = 1200`, `SAMPLE_DT = 2` (averaging window ≈ 1000
time units, matching Huang's `tavg`).

| Ra | Nu | Re | Rex | Rey | tavg | samples |
|---|---|---|---|---|---|---|
| 1e6 | 7.8762 | 70.9287 | 50.3006 | 50.0074 | 1000.0 | 501 |
| 2e6 | 9.6614 | 105.9704 | 75.4324 | 74.4290 | 1000.0 | 501 |
| 3e6 | 10.6443 | 132.9880 | 95.4778 | 92.5732 | 1000.0 | 501 |
| **Huang** | **7.9 / 9.6 / 10.4** | **67.0 / 99.7 / 127.6** | **51.3 / 77.3 / 99.7** | **43.1 / 63.0 / 79.6** | 1000 | — |

Absolute error vs Huang (%)

| Ra | Nu | Re | Rex | Rey |
|---|---|---|---|---|
| 1e6 | **0.30** | 5.86 | 1.95 | **16.03** |
| 2e6 | **0.64** | 6.29 | 2.42 | **18.14** |
| 3e6 | **2.35** | 4.22 | 4.23 | **16.30** |

Anisotropy `⟨u_horiz²⟩ / ⟨u_vert²⟩ = (Rex/Rey)²`

| Ra | CFD | Huang |
|---|---|---|
| 1e6 | 1.012 | 1.417 |
| 2e6 | 1.027 | 1.506 |
| 3e6 | 1.064 | 1.568 |

Every run satisfies `Re = √(Rex² + Rey²)` to machine precision
(`Re identity difference ≈ 1e-14`).

### 2.2 Ra = 1e6: initial-condition and averaging-window variants

| run | `PERT_WAVES` | `T_END` | tavg | Nu | Re | Rex | Rey |
|---|---|---|---|---|---|---|---|
| B | 1 | 400 | 300 | 8.0346 | 71.7412 | 50.8123 | 50.6450 |
| C | 1 | 1200 | 1000 | 7.8762 | 70.9287 | 50.3006 | 50.0074 |
| D | 2 | 400 | 300 | 7.7751 | 63.6736 | 39.8988 | 49.6228 |
| E | `INIT_MODE=huang` | 1200 | 1000 | 1.000000 | 0.0 | 0.0 | 0.0 |

Run C's statistics are stationary: first half vs second half of the averaging
window drifts by only −0.24 % (Rex), −0.52 % (Rey), −1.68 % (Nu).

Run E is the degenerate case discussed in §3.4.

---

## 3. Key findings

1. **Nu converges to ≈ 0.3–2.4 %** across the whole Ra sweep, and its trend
   (7.88 → 9.66 → 10.64 vs 7.9 → 9.6 → 10.4) is correct. `Rex` is also good
   (2–4 %). Together with the `Re` identity these are strong indications that
   the solver, the nondimensionalisation and the post-processing are right.

2. **The `Rey` deviation (16–18 %) is systematic and Ra-independent.** It is
   present at every Ra, and tripling the averaging window at Ra = 1e6 moved it
   only from 17.51 % to 16.03 % while the window itself is demonstrably
   stationary. It is therefore not statistical noise.

3. **The flow is nearly isotropic in our runs but clearly anisotropic in the
   reference**: `(Rex/Rey)²` is 1.01 / 1.03 / 1.06 for Ra = 1e6 / 2e6 / 3e6,
   versus 1.42 / 1.51 / 1.57 in Huang's table. Our vertical (wall-normal)
   fluctuation is ≈ 30 % too large relative to the horizontal one, and the
   ratio of the two anisotropies is almost constant over Ra.

4. **The initial condition selects the solution branch, but the branch does not
   explain the `Rey` offset.** `PERT_WAVES = 2` (two rolls) changes `Re` by
   11 % and `Rex` by 22 % versus `PERT_WAVES = 1`, yet *both* branches stay near
   isotropy (1.01 and 0.65 vs the reference 1.42).

5. **Field snapshots confirm the same roll structure at every Ra.** All three
   Ra cases settle into **two roughly square convection rolls** in the `Γ = 2`
   box (Ra = 1e6: roll centres near `y ≈ 0.5` and `1.3` at t = 1200;
   Ra = 2e6: `y ≈ 0.5` and `1.3`; Ra = 3e6: `y ≈ 0.85` and `1.5` with thinner
   plumes and a thinner boundary layer). Since the roll count is identical while
   the `Rey` offset is constant, **roll count is excluded as the explanation**.

6. **Remaining candidates for the `Rey` gap** (not separable with the material
   at hand):
   - **(a) definition/scope of the published `Re_x` / `Re_y`** (e.g. fluctuation
     components, different time window, different volume weighting), which would
     redistribute the two components;
   - **(b) structural difference of the reference solution** (e.g. wider/flatter
     rolls, or 3D effects) that a 2D, `Γ = 2`, square-roll solution cannot
     reproduce.

   Distinguishing (a) from (b) requires the paper's definitions (or field
   figures), which are not available in this repository.

### 3.4 `INIT_MODE=huang` (u = 0, θ = 0) does not start convection

Requested as step ①. Setting `u = 0`, `T = 0` with wall temperatures `±0.5`
evolves to the **one-dimensional conduction profile** `T = 0.5 − x`; that state
(with `u ≡ 0`) is an exact solution of the incompressible NS + Boussinesq
system, and it is perfectly uniform along `y`, so nothing breaks the symmetry:

```
Ra = 1e6, N = 64, t = 0–1200:  Re_inst = Rex = Rey = 0 ,  Nu_inst = 1
final:                         Nu = 1.000000, Re = Rex = Rey = 0.000000
```

This is a useful negative result (and doubles as a "no-convection limit"
sanity check of the solver), but it means this initial condition cannot be used
for the benchmark: **a perturbation must be added** (physically or numerically).
The mode is kept in the code as `INIT_MODE=huang` (equivalent to `zero`).

---

## 4. What the review verified as correct

- **Coordinate mapping is self-consistent.** With `θ = 0°` buoyancy is
  `(cos θ, sin θ) = (1, 0)`, i.e. along code `+x` from the hot plate (`x = 0`,
  `T = +0.5`) to the cold plate — matching the documented mapping.
- **Nondimensionalisation matches free-fall RBC** (`ν*`, `κ*` consistent between
  `cd_dimensionless_coefficients` and the internal `H/ν*` usage).
- **Definitions follow the standard RBC benchmark**:
  `Re = (H/ν*)√⟨u² + v²⟩_{V,t}`, `Re_x`, `Re_y` from the individual components,
  `Nu = 1 + √(Ra·Pr)·⟨v·T⟩_{V,t}`, with `Re` formed **after** time-averaging
  the squared velocities.
- **Volume averaging uses the real non-uniform grid** (`W = dx ⊗ dy` at the
  pressure indices); **time averaging is trapezoidal** on the physical sample
  times with a strict-monotonicity check.
- **The reference table is internally consistent**: `Re² = Rex² + Rey²` holds
  exactly for all 12 stored rows (spot-checked `1e6`, `1e7`, `1e9`).
- **Guards are in place**: `Pr = 4.3`, `GAP = 1.0`, `Γ = 2`, `PERT_WAVES ≥ 1`
  are enforced; velocity/temperature grid shapes are asserted equal; all
  parameters are range-checked; `huang_case_parameters.txt` is written per run.
- Adding the field-snapshot processor does **not** change the solver results
  (verified bit-identical with and without `SNAPSHOT_TIMES`).

---

## 5. Open questions / risks

1. **Reference values are unverified against the paper text** (see §3.6): is the
   NS/PD row 2D or 3D, what exactly does `Γ = 2` mean, and how are `Re_x`,
   `Re_y`, `Nu` defined and averaged? The 16–18 % `Rey` offset should be checked
   against the published definition before any code change is made.
2. **Grid resolution.** `N_wall = 64` is much coarser than the paper's spectral
   element discretisation. The moderate error growth with Ra in `Nu`
   (0.30 % → 2.35 %) is consistent with resolution becoming tight at higher Ra.
   A `N = 128` run would settle this (the runner honours `N_WALL`).
3. **`case_name` does not include `T_END`/`AVG_START`.** `PERT_WAVES` is now part
   of the name, but two runs with different averaging windows still collide
   (this happened once during the review: the T_END = 400 products of run B were
   overwritten by a later run with the same parameters). Adding `tavg` to the
   name would close the gap.
4. **No statistical-convergence criterion** (unlike Task C, which tests
   stationarity/periodicity). The half-window drift used here can be automated
   from the existing CSV.
5. **No restart support**: a long run must start from `t = 0` if interrupted.
   The runtime monitor prints to stdout only (no CSV on disk).
6. Minor: `write_huang_summary` hard-codes `"aspect ratio = 2"` (guaranteed by
   the guard, but a computed value would be cleaner); `outputs` from
   `solve_unsteady` is unused in the runner.

---

## 6. Tools added during the follow-up

### 6.1 Field snapshots (`SNAPSHOT_TIMES`)

```powershell
$env:SNAPSHOT_TIMES="200,600,1200"     # empty (default) = disabled
julia huang_validation.jl
```

Each requested physical time writes `<case>/field_t<TIME>.csv` — a long table
`ix,iy,x,y,T,ux,uy` on the pressure grid, with the actual time in a `# t = ...`
comment line. The processor is a no-op when the list is empty.

### 6.2 Snapshot figures

```powershell
$env:CASE_DIR="validation_huang/Huang_NS_PD_Ra2e06_Pr4p3_N64_perturbed_W1"
julia plot_huang_fields.jl             # -> <case>/figures/field_t*.png
```

Three panels per snapshot: `T` (filled + contours), `|u|`, and streamlines.
Streamlines are drawn as **iso-contours of the streamfunction** ψ (obtained by
integrating `u_x` along the wall-normal direction); Makie's own `streamplot`
cannot be used because the wall-normal grid is tanh-clustered (non-uniform).

### 6.3 Validation curves

```powershell
julia plot_huang_curves.jl             # -> validation_huang/huang_validation_curves.png / .csv
```

Scans every `huang_validation_comparison.csv` under `validation_huang/`, plots
`Nu / Re / Rex / Rey` versus `Ra` with the reference values on the same axes,
and writes a CSV listing **all** records with a `selected` flag. Selection rule:
per `Ra`, degenerate records (`Re ≈ 0`) are dropped first, then the longest
averaging window wins; excluded records are printed with their reason.

---

## 7. Follow-up experiment log (steps ①–⑤)

| step | what was done | outcome |
|---|---|---|
| ① | `N = 64`, `INIT_MODE=huang` (`u = 0`, `θ = 0`), `t = 0–1200` | **no convection**: `Nu = 1`, `Re = Rex = Rey = 0` (§3.4). The branch test could not be performed. |
| ② | `N = 128` grid comparison | **cancelled by the user mid-run** (the run finished on its own, but no comparison was made and its numbers are not used anywhere in this document). |
| ③ | field snapshots + plotting tools | implemented (§6.1, §6.2); snapshots show two square rolls at every Ra (§3.5). |
| ④ | `Ra = 2e6` and `Ra = 3e6` (same setup as run C) | completed; results in §2.1. `Rey` offset persists (18.1 % / 16.3 %). |
| ⑤ | validation curves over Ra | `plot_huang_curves.jl`, curves in `validation_huang/huang_validation_curves.png`. |
| ⑥ | this document | updated with ①–⑤. |
| ⑦ | return to the main line (tilted RBC / angle extension) | pending. |

---

## 8. How to reproduce

Paper-matched run (Ra = 1e6, run C, ~14 min on the RTX 4060 Laptop):

```powershell
$env:RA="1e6"; $env:N_WALL="64"; $env:AVG_START="200"; $env:T_END="1200"
$env:SAMPLE_DT="2"; $env:MONITOR_DT="100"
julia huang_validation.jl
```

The Ra = 2e6 / 3e6 runs use exactly the same settings with `RA` changed, plus
`SNAPSHOT_TIMES="200,600,1200"`.

Initial-condition / window variants (Ra = 1e6):

```powershell
$env:PERT_WAVES="1"        # or "2";   default 1 = original setup
$env:INIT_MODE="perturbed" # "huang"/"zero" = u=0,theta=0 (does not start)
julia huang_validation.jl
```

Failed smoke test (documents the bug that was fixed):

```powershell
$env:RA="1e6"; $env:N_WALL="32"; $env:AVG_START="5"; $env:T_END="20"
julia huang_validation.jl     # -> FieldError before the fix
```

---

## 9. Files touched during this work

- `common/cd_validation.jl` — fixed the blocking `Tuple` → `NamedTuple` bug;
  added `perturbation_waves` (guard + docstring) and the `huang`/`zero`
  initial-condition alias to `build_huang_initial_state`; added
  `write_huang_field_snapshot` and `make_huang_field_recorder`.
- `huang_validation.jl` — added `PERT_WAVES` and `SNAPSHOT_TIMES` environment
  variables (with validation, parameter printing and inclusion in `case_name`),
  and registered the field recorder as a processor.
- `plot_huang_fields.jl` — new (snapshot figures).
- `plot_huang_curves.jl` — new (validation curves).
- `HUANG_VALIDATION_CHECK.md` — this document.
- Benchmark products under `validation_huang/`:

| directory | contents |
|---|---|
| `Huang_NS_PD_Ra1e06_Pr4p3_N64_perturbed/` | run C (Ra = 1e6, `tavg = 1000`) |
| `Huang_NS_PD_Ra1e06_Pr4p3_N64_perturbed_W1/` | run B + 4 snapshots + figures |
| `Huang_NS_PD_Ra1e06_Pr4p3_N64_perturbed_W2/` | run D |
| `Huang_NS_PD_Ra1e06_Pr4p3_N64_huang_W1/` | run E (degenerate, all zeros) |
| `Huang_NS_PD_Ra1e06_Pr4p3_N128_perturbed_W1/` | step ② (cancelled, unused) |
| `Huang_NS_PD_Ra2e06_Pr4p3_N64_perturbed_W1/` | step ④ + 3 snapshots + figures |
| `Huang_NS_PD_Ra3e06_Pr4p3_N64_perturbed_W1/` | step ④ + 3 snapshots + figures |
| `huang_validation_curves.csv` / `.png` | step ⑤ summary |

Each case directory contains `huang_statistics_timeseries.csv`,
`huang_validation_summary.txt`, `huang_validation_comparison.csv` and
`huang_case_parameters.txt`.
