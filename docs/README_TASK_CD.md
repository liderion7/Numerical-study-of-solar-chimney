# Task C / Task D implementation for the 2-D chimney model

## Agreed coordinate convention

- `x`: horizontal / wall-normal direction, from hot wall to cold wall.
- `y`: vertical / streamwise direction.
- `theta = 90°`: chimney axis parallel to `y`.
- Temperature is one-dimensional in the conduction reference state:
  `T(x) = 0.5 - x/Lx`.
- The velocity needed by Task C/D is the **y velocity component** `v = u_y`.
- Therefore the profile to analyse is a function of x:

  `vbar(x,t) = (1/Ly) ∫ v(x,y,t) dy`.

Task D then computes

`vbar_t(x) = (1/T) ∫ vbar(x,t) dt`

and

`Vbulk = (1/Lx) ∫ vbar_t(x) dx`.

Both integrations are trapezoidal and use the true coordinates.

## Files

- `common/cd_grid.jl`: Task-C/D geometry and initial state.
- `common/cd_solver.jl`: angle convention and optional streamwise pressure/body force.
- `common/cd_tools.jl`: profile extraction, raw-data recorder, audit, time/spatial integration, window and periodicity tools.
- `run_cd_case.jl`: production CFD data generator.
- `task_c.jl`: data audit + transient/stationary/periodic screening; exports clean data only when suitable.
- `task_d.jl`: final time-averaged V-X profile and `Vbulk`.
- `collect_cd_results.jl`: combines completed angles and finds the maximum `Vbulk`.

## Important model distinction

The previous 90° Ng-style validation used a zero-mass-flux constraint. Do **not** use that constraint for the final Task C/D chimney calculation, because Task D is supposed to measure the net bulk velocity. `boussinesq_cd!` therefore has no zero-flux correction.

The meeting notes mention an additional upward force / pressure-gradient term, but do not give its magnitude. The new code exposes this as `DRIVE_FORCE_Y`; it defaults to zero and should only be assigned a nonzero value when the project specifies one.

## Recommended workflow

### 1. Run one CFD case

Example 90° case with no extra driving force:

```bash
THETA=90 T_END=400 DRIVE_FORCE_Y=0.0 julia run_cd_case.jl
```

If the final chimney model specifies a pressure-gradient/body force, e.g. `F0`:

```bash
THETA=90 T_END=400 DRIVE_FORCE_Y=F0 julia run_cd_case.jl
```

Output:

`results_cd/90degree/vbar_xt_raw.csv`

with columns `time,x,vbar`.

### 2. Run Task C

```bash
CASE_DIR=results_cd/90degree julia task_c.jl
```

Task C checks data quality, independent windows, stationarity and periodicity. If it finds a reliable regime it writes:

`results_cd/90degree/task_c/task_c_clean_vbar.csv`

and sets `Task D ready = true` in `task_c_summary.txt`.

If the flow is clearly periodic but fewer than two complete cycles have been recorded, extend `T_END` and rerun rather than forcing Task D.

A justified manual interval can be supplied as:

```bash
CASE_DIR=results_cd/90degree T_START=200 T_END_STAT=400 julia task_c.jl
```

### 3. Run Task D

Only after Task C passes:

```bash
CASE_DIR=results_cd/90degree julia task_d.jl
```

Outputs include:

- `vbar_t_profile.csv`
- `V_X_time_averaged.png`
- `task_d_result.csv`
- `task_d_summary.txt`

The primary Task-D result is signed `Vbulk`.

### 4. Repeat for angles

Run the CFD + C + D chain for all required angles, then:

```bash
julia collect_cd_results.jl
```

This creates:

- `results_cd/bulk_velocity_vs_angle.csv`
- `results_cd/Vbulk_vs_angle.png`
- `results_cd/optimal_angle.txt`

## Expected 90° low-Ra validation

At 90°:

- buoyancy has no x component and is entirely along y;
- `T=T(x)` is linear at low Ra;
- `v=u_y` is a function of x;
- the steady low-Ra profile is cubic because `nu d²v/dx²` is forced by linear `T(x)`.

`analytic_low_ra_v` in `cd_tools.jl` is included for this optional validation.
