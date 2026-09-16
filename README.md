# task_cd_code — 2-D chimney model (Task C / Task D) and three validation tracks

This directory is organised by **task goal**: one main-line CFD production task (Task C/D) plus three validation tasks (low-Ra analytical validation, the Huang et al. (2022) literature benchmark, and the Wang et al. (2018) tilted-convection literature benchmark).

## Task goals

| Directory | Task goal | Key parameters | Main products |
|---|---|---|---|
| `main_task_cd/` | **Main line**: 2-D chimney natural convection, swept over the inclination angle θ. Task C decides which sampling interval is statistically reliable; Task D reports `Vbulk(θ)` and finds the optimal angle. | `Ra=1e7`, `Pr=0.71`, `drive_force_y=0`, no zero-mass-flux constraint | `data/<case>/task_c`, `data/<case>/task_d`, `data/bulk_velocity_vs_angle.csv`, `data/optimal_angle.txt` |
| `validation_lowRa/` | **Low-Ra validation**: verify that the solver agrees with the analytical solution (at `Ra=1000`, `T(x)=0.5-x` and `v` is the analytical cubic), plus flow-field visualisation | `Ra=1000`, `Pr=0.71`, `T_END=100`, `drive_force_y=0` | `data/<θ>degree/validation/`, `data/lowRa_fields*`, `data/lowRa_movie/`, `data/Ra1000_summary/` |
| `validation_huang/` | **Literature benchmark**: against the NS-plate + periodic-sidewall RBC case of Huang et al. (2022) Table 1 | `Pr=4.3`, `Γ=L/H=2`, `θ=0°` | `data/Huang_*/`, `data/huang_validation_curves.{csv,png}` |
| `validation_wang/` | **Literature benchmark (unfinished)**: against the Nu of Wang et al. (2018) tilted convection Table I (Γ=4, β=0°/3°, n=4) | Hard-coded inside the scripts: `Ra=1e7`, `Pr=0.71`, `n=128`, `T_END=100` | `Wang_beta{0,3}deg_*` in the script's own directory (not produced yet) |
| `shared/common/` | Code shared by all task lines (mesh, solver, Task C/D numerics, plotting conventions, Huang validation library) | — | — |
| `docs/` | Original project documents (kept verbatim from before the tidy-up; their example paths are mapped in the migration table below; still written in Chinese) | — | — |

## Directory structure

```
task_cd_code/
├── README.md                    This file: task goals, structure, how to run, known issues
├── README.zh-CN.md              Simplified-Chinese original of this file
├── reasonix.toml                Tool configuration
├── docs/                        Original documents (not rewritten)
│   ├── README_TASK_CD.md            Main-line Task C/D conventions and workflow
│   ├── FINAL_TASK_CD_SUMMARY.md     Main-line angle-sweep summary (note: source of the numbers, see "Known issues")
│   └── HUANG_VALIDATION_CHECK.md    Review record and conclusions for the Huang benchmark
├── shared/common/               Shared library
│   ├── cd_grid.jl                   2-D geometry and initial fields (tanh wall clustering, periodic y)
│   ├── cd_solver.jl                 Boussinesq force e(θ)=(cosθ, sinθ), non-dimensional coefficients
│   ├── cd_tools.jl                  Task C/D numerics (integration, auditing, periodicity detection, recorders)
│   ├── cd_plots.jl                  Plotting conventions (x axis reversed 1→0, velocity on the vertical axis)
│   └── cd_validation.jl             Huang benchmark reference table plus statistics/snapshot tools
├── main_task_cd/                Main-line code + data
│   ├── run_cd_case.jl               [1] CFD production (the only CFD generator, shared by both task lines)
│   ├── task_c.jl                    [2] Data audit + stationarity/periodicity decision, exports the clean interval
│   ├── task_d.jl                    [3] Time-averaged profiles and Vbulk / Vabs / Uwind
│   ├── task_d_manual.jl             Task D for a single manually chosen mature cycle
│   ├── collect_cd_results.jl        [4] Collect all Ra≈1e7 cases → Vbulk(θ)
│   ├── replot_task_cd.jl            Redraw figures only, without touching decision files
│   ├── run_n128_sweep.ps1           Batch: θ=0,15,30,45,60,75 running 1→2→3 in sequence
│   └── data/                        Main-line data (see below)
├── validation_lowRa/            Low-Ra validation code + data
│   ├── validate_low_ra.jl           Part A analytical self-checks (A1–A7) + Part B CFD vs analytical
│   ├── plot_low_ra_fields.jl        Per-angle snapshots → validation_2d/ triples + overview figure
│   ├── make_low_ra_video.jl         Snapshot sequence → mp4/gif animation
│   ├── collect_ra1000_results.jl    Collect all Ra=1000 cases → Ra1000_summary/
│   └── data/
├── validation_huang/            Huang benchmark code + data
│   ├── huang_validation.jl          Benchmark runner (Ra=1e6/2e6/3e6)
│   ├── plot_huang_fields.jl         Snapshot triples (T, |u|, stream function)
│   ├── plot_huang_curves.jl         All cases → Nu/Re/Rex/Rey vs Ra curves
│   └── data/
└── validation_wang/             Wang (2018) tilted-convection benchmark code (self-contained, not run yet)
    ├── Nu0degree_myBC_myPerturbation.jl   β=0°: own geometry/BC/initial condition + wall Nu
    └── Nu3degree_myBC_myPerturbation.jl   β=3°: against Wang Table I Nu=13.06
```

Data directory contents:

```
main_task_cd/data/
├── coarse_N64/          Ra=1e7, N=64, T_END=550 coarse-grid angle sweep (0,15,30,45,60,60_T300,75 degree)
├── 15degree/  75degree/ Ra=1e7, N=128, T_END=550 single-point supplements
├── 90degree_Ra1e7/      Ra=1e7, N=128 (Task C products only)
├── 90degree_Ra1e7_T550/ Ra=1e7, N=128, T_END=550 (Task C + Task D + manual)
├── bulk_velocity_vs_angle.csv  summary produced by collect_cd_results.jl
├── optimal_angle.txt / Vbulk_vs_angle.png
validation_lowRa/data/
├── 0,10,20,30,40,45,50,60,70,80,90degree/  Ra=1000, N=128, T_END=100 (some contain field_t*.csv)
├── 90degree_Ra1e3/      Ra=1000, N=32
├── lowRa_fields/        Ra=1000, N=64 snapshots + validation_2d figures + lowRa_fields_overview.png
├── lowRa_movie/90degree/ Dense snapshots for Ra=1000, N=64 (used for the animation)
├── Ra1000_summary/      Summary produced by collect_ra1000_results.jl (includes copies of the 12 cases)
└── lowRa_Ra1e3_summary.csv, lowRa_Ra1e3_sweep*.log
validation_huang/data/
├── Huang_NS_PD_Ra{1e06,2e06,3e06}_Pr4p3_N64_perturbed[_W1,_W2]/
├── Huang_NS_PD_Ra1e06_Pr4p3_N64_huang_W1/ (degenerate case: u=0 never starts convection)
├── Huang_NS_PD_Ra1e06_Pr4p3_N128_perturbed_W1/ (cancelled by the user, not used for any conclusion)
└── huang_validation_curves.{csv,png}
```

## How to run

All paths are relative to the script itself (`@__DIR__` / `$PSScriptRoot`), so the whole directory can be moved as a unit. Parameters are always passed through environment variables.

Main line (Task C/D):

```powershell
cd main_task_cd
$env:THETA="90"; $env:T_END="550"; $env:RA="1e7"; $env:N_WALL="128"
julia run_cd_case.jl                 # -> data/90degree/{vbar_xt_raw.csv,...}
$env:CASE_DIR="data/90degree"; julia task_c.jl    # decides regime, writes task_c/task_c_clean_vbar.csv
julia task_d.jl                                    # -> data/90degree/task_d/... + Vbulk
julia collect_cd_results.jl                        # -> data/bulk_velocity_vs_angle.csv, optimal_angle.txt
```

Low-Ra validation (note: it reuses `main_task_cd/run_cd_case.jl`, but the results root must be pointed at this task):

```powershell
cd main_task_cd
$env:RA="1000"; $env:THETA="90"; $env:T_END="100"; $env:N_WALL="128"
$env:RESULTS_ROOT="..\validation_lowRa\data"
julia run_cd_case.jl
$env:CASE_DIR="..\validation_lowRa\data\90degree"; julia task_c.jl; julia task_d.jl
cd ..\validation_lowRa
julia validate_low_ra.jl             # Part A + Part B -> data/90degree/validation/
julia collect_ra1000_results.jl      # -> data/Ra1000_summary/
```

Huang benchmark:

```powershell
cd validation_huang
$env:RA="1e6"; $env:N_WALL="64"; $env:AVG_START="200"; $env:T_END="1200"
julia huang_validation.jl            # -> data/Huang_NS_PD_Ra1e06_.../
julia plot_huang_curves.jl           # -> data/huang_validation_curves.png
```

Wang benchmark (`validation_wang/`, self-contained scripts, no output yet):

```powershell
cd validation_wang
julia Nu0degree_myBC_myPerturbation.jl   # β=0° -> Wang_beta0deg_*.mp4/png in this directory
julia Nu3degree_myBC_myPerturbation.jl   # β=3° -> against Wang Table I Nu=13.06
```

## Current results (as of the tidy-up)

- **Main line, Ra=1e7**: every completed angle is `provisional_single_cycle` (single-cycle estimate), **not** a statistically converged value.
  - `coarse_N64` (N=64): 15° `Vbulk=-0.0329`, 45° `-3.2998`, 60° `3.3913`, 75° `0.4866`; 0°/30° are still not mature at `T_END=550` (`not_ready`).
  - 90°, N=128, `T_END=550`: `Vbulk=3.7992`, `Vabs=5.1038`, `Uwind=10.3004`.
  - Therefore "optimal angle = 90°" can only be read as a **coarse trend**, and N=64 vs N=128 even disagree in sign at 45°/60°/75° (see "Known issues").
- **Low Ra, Ra=1000**: all cases are `converged_stationary`; `Vbulk` is a floating-point zero of order 1e-17 (an inevitable consequence of the symmetric counter-flow, not a net flow). The meaningful quantities are `Vabs`/`Uwind`, which grow monotonically with θ and reach `0.1955`/`0.3009` at 90°. Analytical validation: all Part A self-checks pass to machine precision (1e-14–1e-16).
- **Huang benchmark**: at Ra=1e6/2e6/3e6 the Nu errors are 0.30%/0.64%/2.35%, `Rex` is 2–4% off, and `Re=√(Rex²+Rey²)` matches to machine precision; however **`Rey` is systematically 16–18% too large and independent of Ra**, with an anisotropy ratio of ≈1.01–1.06 (paper: 1.42–1.57). The N=128 run was cancelled and is not used for any conclusion.
- **Wang benchmark**: not run yet, no results.

## Known issues (the tidy-up changed no logic, it only records them)

1. **`bulk_velocity_vs_angle.csv` and `docs/FINAL_TASK_CD_SUMMARY.md` disagree**: the former was collected by `collect_cd_results.jl` from the `results_cd/` top-level cases of that time (the N=128 batch: 15° -0.0030, 30° +1.5286, 45° +2.6393, 60° -1.5050, 75° -2.6515), whereas the latter corresponds to the N=64 `coarse_N64` sweep plus the N=128 90° case. The two sets disagree in sign at 45°/60°/75°, which shows the grid is far from converged; moreover the raw case directories behind the former were later overwritten by the Ra=1000 run, so it can no longer be reproduced from the current data.
2. **Case directory names do not include `Ra`/`N_WALL`/`T_END`** (`OUTDIR = data/<θ>degree` in `run_cd_case.jl`), so cases with the same name overwrite each other — this is exactly why the top-level low-Ra directories once overwrote the Ra=1e7 angle-sweep results. Recommendation: put the grid and Ra into the directory name.
3. **No restart/checkpoint**: extending a case means re-running from `t=0`; `T_END=550` is not enough for 0°/30°.
4. In `shared/common/cd_plots.jl` the x/y axis titles of the Task C window-profile figure are explicitly swapped (`xlabel="z axis"`), which is the opposite of the convention in `plot_task_d_profile`; the two figures do not use a consistent axis-title convention.
5. `main_task_cd/data/90degree/task_d/CFD_analytic_comparison.{csv,png}` has no generating code in any script — an orphan product that this tidy-up did not delete.
6. The table in `docs/FINAL_TASK_CD_SUMMARY.md` is now annotated as coming from `coarse_N64`, but the file itself was not rewritten.
7. **The two scripts in `validation_wang/` overlap heavily with `shared/common`**: the tilted buoyancy force, the non-dimensional coefficients (`ν=√(Pr/Ra)`, `κ=1/√(Ra·Pr)`), the `Setup`/`BC`/initial condition, and the Nu statistics and field-plot framework all duplicate the shared library; moreover the angle convention is the **opposite** of `cd_solver.jl` (the scripts use `(sinβ, cosβ)`, the shared library uses `(cosθ, sinθ)` — the same number means something 90° different). The two scripts are also about 95% identical line by line. Not run yet, no output.

## Follow-up: archiving the Wang benchmark scripts

- The two self-contained scripts `Nu0degree_myBC_myPerturbation.jl` and `Nu3degree_myBC_myPerturbation.jl` were moved from the top level of the directory into the newly created `validation_wang/`, alongside `validation_huang/`.
- Neither is referenced by any file and neither has any output (they never produced results), so the move does not affect any existing workflow; the scripts contain no `include` and do not depend on their original location.
- The scripts write their products into their own directory (`output_dir = @__DIR__`), so after the move the video/figures land in `validation_wang/`.

## Migration table (old → new)

| Old path | New path |
|---|---|
| `common/*.jl` | `shared/common/*.jl` |
| `README_TASK_CD.md`, `FINAL_TASK_CD_SUMMARY.md`, `HUANG_VALIDATION_CHECK.md` | `docs/…` (same names) |
| `run_cd_case.jl`, `task_c.jl`, `task_d.jl`, `task_d_manual.jl`, `collect_cd_results.jl`, `replot_task_cd.jl`, `run_n128_sweep.ps1` | `main_task_cd/…` (same names) |
| `validate_low_ra.jl`, `plot_low_ra_fields.jl`, `make_low_ra_video.jl`, `collect_ra1000_results.jl` | `validation_lowRa/…` (same names) |
| `huang_validation.jl`, `plot_huang_fields.jl`, `plot_huang_curves.jl` | `validation_huang/…` (same names) |
| `results_cd/coarse_N64/`, `results_cd/15degree/`, `results_cd/75degree/`, `results_cd/90degree_Ra1e7/`, `results_cd/90degree_Ra1e7_T550/` | `main_task_cd/data/…` (same names) |
| `results_cd/bulk_velocity_vs_angle.csv`, `optimal_angle.txt`, `Vbulk_vs_angle.png` | `main_task_cd/data/…` (same names) |
| `results_cd/{0,10,20,30,40,45,50,60,70,80,90}degree/` (Ra=1000), `results_cd/90degree_Ra1e3/`, `results_cd/lowRa_fields/`, `results_cd/lowRa_movie/`, `results_cd/Ra1000_summary/`, `results_cd/lowRa_Ra1e3_*` | `validation_lowRa/data/…` (same names) |
| `validation_huang/Huang_*/`, `validation_huang/huang_validation_curves.*` | `validation_huang/data/…` (same names) |
| `.reasonix/`, `results_cd/Ra1000_summary.rar` | deleted |
