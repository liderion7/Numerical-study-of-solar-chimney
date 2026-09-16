# FINAL_TASK_CD_SUMMARY

## Simulation setup

- 2-D chimney model (no 3D extension).
- Coordinates: `x` = wall-normal (hot wall `x=0`, `T=+0.5`; cold wall `x=1`, `T=-0.5`), `y` = streamwise/vertical, periodic.
- `theta = 90 deg` => chimney axis parallel to `y`.
- Production parameters: `Ra = 1e7`, `Pr = 0.71`, `drive_force_y = 0.0`, `zero_mass_flux = false`.
- Coarse sweep grid: `Nx = 64`, tanh wall clustering `1.2`, periodic `Ly = 2`.
- 90-degree baseline: `Nx = 128`, `T_END = 550`.
- Solver: IncompressibleNavierStokes.jl v5 on GPU (`NVIDIA GeForce RTX 4060 Laptop GPU`).
- Sampling: `vbar(x,t)` every `0.2` time units; runtime monitor every `5.0` time units.

## Task C methodology

- Data audit: NaN/Inf count, duplicate time/x, strict monotonicity, `dt`/`dx` statistics.
- Stationary candidate: independent 20-unit windows; requires >= 3 late windows with profile delta `< 0.05` and RMS drift `< 0.05`.
- Periodic candidate: exclude startup transient (first RMS peak + min separation), then require >= 3 post-transient peaks, period CV `< 0.10`, peak-profile correlation `> 0.98`, and consecutive cycle-averaged Vbulk relative difference `< 10%`.
- Fallback: if a single mature peak-to-peak cycle exists, mark `provisional_single_cycle` and allow Task D on that one cycle. This is never reported as converged.

## Task D methodology

- `vbar_t(x) = (1/T) * trapezoidal_time_integral(vbar(x,t))` on the true adaptive time grid.
- Add no-slip wall values `vbar_t(0)=vbar_t(1)=0`.
- `Vbulk = (1/Lx) * trapezoidal_x_integral(vbar_t(x))` on the true tanh x grid.
- Also reported: `Vabs = (1/Lx)*integral(|vbar_t|)`, `Uwind = max(|vbar_t|)`.

## Angle sweep table

| theta | Vbulk | Vabs | Uwind | t_start | t_end | regime | status | grid_n | runtime (s) |
|---|---|---|---|---|---|---|---|---|---|
| 0 | — | — | — | — | — | not_ready | not_ready | 64 | 111.0 |
| 15 | -0.03291 | 0.07006 | 0.14714 | 466.218 | 545.013 | provisional_single_cycle | provisional_single_cycle | 64 | 218.0 |
| 30 | — | — | — | — | — | not_ready | not_ready | 64 | 738.4 |
| 45 | -3.29979 | 4.15022 | 8.29770 | 294.001 | 453.600 | provisional_single_cycle | provisional_single_cycle | 64 | 1227.0 |
| 60 | 3.39127 | 4.38756 | 8.77930 | 258.200 | 436.200 | provisional_single_cycle | provisional_single_cycle | 64 | 1939.7 |
| 75 | 0.48664 | 2.80133 | 4.83283 | 261.601 | 522.000 | provisional_single_cycle | provisional_single_cycle | 64 | 1368.8 |
| 90 | 3.79920 | 5.10381 | 10.30043 | 281.600 | 443.400 | provisional_single_cycle | provisional_single_cycle | 128 | (earlier run) |

## Vbulk vs angle

- Coarse optimum: `theta_opt = 90 deg`
- `max Vbulk = 3.7992`
- Values available: 15, 45, 60, 75, 90 deg (0 and 30 deg did not form a mature cycle by `T_END=550`).

## Optimal angle

- `theta_opt = 90 deg` (from the completed sweep).
- 60 deg is the second largest completed value: `Vbulk = 3.3913`.

## Computational limitations

- No restart/checkpoint support in the current runner: every extension was a full restart from `t=0`.
- `T_END=550` was insufficient for 0 and 30 deg to produce a mature cycle; they are marked `not_ready`, not converged.
- All completed sweep results are `provisional_single_cycle` estimates from one mature cycle, not fully statistically converged values.
- 90 deg uses `Nx=128` while other angles use `Nx=64`; a direct grid-difference check at the optimum is still needed.

## Recommended final result

- Coarse trend: `Vbulk(θ)` peaks near `90 deg` with `Vbulk ≈ 3.799` (provisional).
- 60 deg is close (`3.391`), so the true optimum may lie between 60 and 90 deg.
- Recommended next step: fine sweep around `75–90 deg` (e.g. `80, 85, 90`) and a high-resolution `Nx=128` confirmation of the optimum, then a grid-difference check.
- Do not use the 0/30 deg cases (or any provisional value) as converged final results.
