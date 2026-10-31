# DiffuseCode: correct and cheaper structure refinement (REFINE, RMC)

**Branch:** [`perf/refine`](https://github.com/jmsweng/DiffuseCode/tree/perf/refine) on `jmsweng/DiffuseCode`, a fork of `tproffen/DiffuseCode`
**Base:** `master` @ `6fd5f33` (upstream, 28 Oct 2026)
**Status (31 Oct 2026):** every commit was checked against the one before it in the cloud (gfortran 13.3) and on a laptop (Ubuntu 22.04 under WSL2, gfortran 11.4), with identical results. MPI was checked in both environments with OpenMPI.

## Summary

A REFINE cycle runs the user macro once at the current parameters and again for every derivative (2 per parameter for 3-point, 4 for 5-point). In a real refinement each run is a full DISCUS simulation, so the number of macro runs is the cost. REFINE also depends on the random numbers being identical in all runs of a cycle, and on every MPI worker returning its whole result.

This branch:

1. fixes five bugs that gave wrong derivatives, wrong convergence figures or no result at all: MPI step sizes (1), the seed restore (2), RMC with lots (3), the `dP/sigma` measure (4), and two MPI problems found while testing (8, 9);
2. removes macro runs that do not change the result: the second run at the centre point (5) and derivatives at rejected steps (6);
3. adds an opt-in one-sided derivative that needs half the runs (7).

Items 5 and 6 give bit-identical results. Items 1, 2, 3, 4, 8 and 9 change results only where the old code was wrong.

## Changes

| # | Change | Files | When active | Effect |
|---|---|---|---|---|
| 1 | **MPI derivative step sizes** (bug fix). In `refine_calc_deriv_mpi` three assignments wrote `delta = ...` (the whole array) instead of `delta(k) = ...`. The 5-point derivative divides by `delta(k)` after all children return and so used the last parameter's step. | `refine/prog/refine_run.f90` | `points:5` under MPI. 3-point was not affected | First cycle of a 2-parameter `points:5` fit under MPI: chi²/(N-P) 0.50458, equal to serial; 244 before |
| 2 | **Seed restore.** `put_seeds` restored the `RANDOM_NUMBER` state but not the cached second value of the Box-Muller generator `GASDEV` (`iset`/`gset`). REFINE restores the seeds so the centre and all derivative points see the same random numbers; with a stale cache the first Gaussian number depended on what ran before, and the derivatives were noisy | `lib_f90/random_state_mod.f90`, `refine/prog/refine_random_mod.f90` | Any model that draws Gaussian random numbers after a seed restore (`gran()`, thermal displacements, MC) | Fit with `gran(0.02)` noise: 21 cycles to a wrong minimum (chi² 7992) before, 4 cycles to the true one (chi² 76) after. Deterministic fits unchanged |
| 3 | **RMC with lots.** `rmc_inlot` tested the Fourier lot type `ilots` instead of `rmc_ilots`, so no move was made unless `fourier lots` was set; its periodic wrap could never trigger; `rmc_strucf` added every move to lot 1; `rmc_inten` counted lot 1 twice; `rmc_calc_scat` subtracted the average over `num(1)×num(1)` points instead of `num(1)×num(2)` | `discus/prog/rmc_sup_mod.f90` | `set lots` in RMC. Without lots nothing changes | Cu 6×6×6, one hk0 plane, 300 accepted moves: incremental and freshly calculated chi² agree to 8e-8 (they differed by 0.5 % with `fourier lots`) |
| 4 | **`dP/sigma` measure.** `refine_mrq` divided the parameter shift by `SQRT(cl(k,k))`, but after `mrqmin` `cl` holds the curvature matrix of the new point, not its inverse. The printed `MAX(dP/sig)` had the wrong units, and convergence criterion 1 (`pshift`, default 0.005) could be met while parameters still moved by many sigma | `refine/prog/refine_run.f90` | Always | Test fit, cycle 0: printed 5e-3, actual shift 215 sigma. A fit that used to stop on criterion 1 may now run more cycles (correctly) |
| 5 | **Skip the repeated centre run.** After the derivatives, `refine_theory` ran the macro a second time at the centre point to leave the program state at the current parameters. Now done only with a per-cycle plot macro (`run ..., plot:`) or under MPI; otherwise the random generator is restored to the state the repeated run would have left. Also fixes a read past the array bound for a single parameter (`gl_is_der(1) .and. gl_is_der(2)`) | `refine/prog/refine_run.f90`, `refine/prog/refine_random_mod.f90` | Always, except with plot macro or MPI | −1 run per cycle, bit-identical |
| 6 | **Derivatives only for accepted steps.** `mrqmin` computed all derivatives at every trial point and discarded them when chi² did not improve. The trial point is now evaluated with one run; only an accepted step computes the derivatives, reusing the trial run as their centre | `refine/prog/refine_run.f90` | Always | A rejected step costs 1 run instead of 2N+1 (3-point) or 4N+1 (5-point), bit-identical |
| 7 | **`points:2`** (new, opt-in). `(f(p+h) − f(p))/h` needs one extra run per parameter instead of two. The default stays `points:3`. Help text and the error message for invalid `points:` values are updated; `refine_calc_deriv_mpi` handles it | `refine/prog/refine_run.f90`, `refine_add_param.f90`, `refine_err_appl.f90`, `appl_ref.hlp` | Only with `newparam ..., points:2` | N instead of 2N runs per Jacobian. Suited to smooth parameters such as scale factors (error ∝ h instead of h²) |
| 8 | **KUPLOT result for the MPI master** (bug fix). Each worker returns its last KUPLOT data set as `DISCUS_SUITE_DERIVATIVES/data.nnnn`. `kuplot_exit` called `check_form` with no limits, which for the format H5 leaves the limits at zero, and `do_save` then wrote only the range 0..0, the single point at x = 0, for any data set not read from an HDF5 file. The master loaded a 1×1×1 data set and read every derivative run as the original data | `kuplot/prog/kuplot_exit.f90` | MPI refinement of a KUPLOT curve or NIPL map | Before, no such MPI refinement finished (`All derivatives are zero`) |
| 9 | **HDF5 write** (bug fix). `hdferr` was tested uninitialised when the HDF5 library was already open, so a worker's second write failed with a false `Error opening file`. `hdferr = 0` is now set before the tests (two places) | `lib_f90/hdf_write.f90` | Any HDF5 write after the first one | Multi-cycle MPI fits work |

Also changed: the help text for `newparam` (`points:2`).

### Commits

| Commit | Description |
|---|---|
| `2e95b93` | REFINE: fix MPI derivative step sizes being overwritten |
| `df2b814` | Make seed restore reproduce Gaussian random numbers |
| `a299cfa` | RMC: fix refinement with lots |
| `368094c` | REFINE: compute the dP/sigma convergence measure from the parameter sigmas |
| `5d65bcb` | REFINE: skip the repeated centre run after the derivatives |
| `863bb4f` | REFINE: compute derivatives only for accepted steps |
| `f40b2b5` | REFINE: add two-point (one-sided) derivatives, `points:2` |
| `d1edc81` | KUPLOT: save the complete data set when it is returned to the REFINE master |
| `0d73677` | HDF5 write: initialise hdferr before testing it |

9 source files changed, +232 / −70 lines. Each commit builds on its own; 5 → 6 → 7 depend on each other. Any prefix of the series can be merged.

## Cost

Number of macro runs for a 4-parameter Gaussian fit (fixed seeds, KUPLOT-only model, so it runs in seconds). "Base" is `master`; the other columns are the cumulative state after that commit.

| case | base | after 4 | after 5 | after 6 | after 7 |
|---|---|---|---|---|---|
| easy start | 52 | 52 | 46 | 46 | 46 |
| far start (6 of 15 steps rejected) | 162 | 162 | 145 | 97 | 97 |
| 25 fixed cycles (13 rejected) | 262 | 262 | 235 | 131 | 131 |
| with `gran()` noise | 222 | 52 | 46 | 38 | 38 |
| noise, 25 fixed cycles | 262 | 262 | 235 | 91 | 91 |
| 5-point derivatives | 92 | 92 | 86 | 86 | 86 |
| `points:2` (easy / far / noise) | error | error | error | error | 26 / 63 / 22 |

Parameters, esd's and chi² of the deterministic cases are the same in every column; only the number of runs changes. The "noise" rows change with commit 2 because the fit now converges. The saving grows with the share of rejected steps and with the number of parameters, and is largest when the macro is slow.

## Validation

### Serial checks (commits 1–7)
For each commit, the binary built before it and the one built after it run the same fits. Compared: cycle tables, parameters, esd's, chi² and macro-run counts. Fits: a 4-parameter Gaussian (easy start, far start, 25 fixed cycles, `gran()` noise, 5-point, `points:2`) and an RMC test on Cu 6×6×6 with box lots and with `fourier lots`. Result: all checks pass, in both environments, with identical run counts.

### MPI checks (commits 1, 8, 9, and the speed-ups under MPI)
An MPI build of the tip runs four fits serially and with `mpirun -n 5` (4 workers, so each worker handles several derivative runs per cycle):

| case | what it tests |
|---|---|
| `det` | easy start, 3-point derivatives |
| `far` | far start with rejected steps (derivatives only for accepted steps) |
| `p5` | 5-point derivatives (16 derivative runs per cycle) |
| `p2far` | `points:2`, far start |

Pass criteria: same number of cycles, parameters equal to 1e-4 (relative), no error lines, worker log files present. A second MPI build with item 1 put back is run on `p5`: its chi²/(N-P) in cycle 0 is 365.29 against 109.18 serial, so the check detects the old bug. Result: all four cases pass and the old bug is detected, in both environments, with the same numbers.

Two things to know when running MPI refinements by hand:
- Without a terminal the program cannot find its parent process (it runs `ps j`), silently decides MPI is off and runs the whole fit serially. Start `mpirun` under `script -qec "..." /dev/null` and check that `DISCUS_SUITE_DERIVATIVES/LOGFILE.*` exist.
- MPI workers start with empty KUPLOT, so the work macro must create its own data (`func ..., -5, 5, 0.1`) and not use the loaded data set.

### Not tested
- MPI with DISCUS (3D) models (`lh5` data); the checks used a KUPLOT model.
- Multi-data-set REFINE (`mrqcof` loops over data set 1's size for all sets); left alone.
- RMC with ellipsoidal lots; box lots and box lots with `fourier lots` were tested, on one hk0 plane.
- `points:2` and `points:5` with more than one data set.

## Behaviour changes to be aware of

Which results can change, in short:

| Item | Serial REFINE | MPI REFINE | RMC | Numbers differ from `master`? |
|---|---|---|---|---|
| 1 | — | `points:5` | — | yes, where `master` was wrong |
| 2 | models with Gaussian random numbers | same | — | yes, where `master` was wrong |
| 3 | — | — | `set lots` | yes, where `master` was wrong |
| 4 | printed column, stop criterion 1 | same | — | printed values: yes; fit path: only if it stopped on criterion 1 |
| 5, 6 | number of macro runs | same | — | no |
| 7 | only if `points:2` is used | same | — | no |
| 8, 9 | — | every MPI refinement of a KUPLOT curve | — | `master` gave no result |

### Where the numbers change (bug fixes)

- **Item 1 (MPI step sizes).** Only `points:5` under MPI is affected. With `master`, every 5-point derivative was divided by the step of the last parameter, so with parameters of different size the Jacobian, the steps and the esd's were wrong. An MPI fit that converged anyway will take a different path and may end at slightly different parameters and a different number of cycles. Serial fits and `points:3` under MPI are unchanged.
- **Item 2 (seed restore).** Any refinement, or any macro calculation using the seed restore, that draws Gaussian random numbers changes. This includes `gran()` in the macro, thermal displacements and Monte Carlo steps that call the Gaussian generator. Before, the first Gaussian number after a restore depended on what had run before, so the centre run and the derivative runs saw different noise and the derivatives were noisy. Now they see the same numbers. A fit that used to stop early on a wrong minimum can now run more cycles to a different, better result. Fits with uniform random numbers or none are bit-identical. Fixed seeds still give reproducible results, but not the same numbers as before this fix.
- **Item 3 (RMC lots).** Only RMC runs with `set lots` are affected; without lots the code path is unchanged. Before, no move was made unless `fourier lots` was also set, so those runs did not refine anything; with `fourier lots` the lot sums were wrong by about 0.5 %. RMC runs with lots now move atoms, accept and reject steps, and reach different chi² values. Results of earlier RMC runs with lots should not be reused.
- **Items 8 and 9 (MPI results).** Only MPI refinements are affected. Before, a refinement of a KUPLOT curve under MPI stopped with `All derivatives are zero`, and from the second cycle on a worker could fail with `Error opening file`. Now they run. No result that was correct before changes.

### Where the printed output or the stopping point changes

- **Item 4 (`MAX(dP/sig)`).** The column in the cycle table is now the largest parameter shift in units of the parameter sigma (from the covariance matrix), as the name says. Before, it was divided by the square root of the curvature matrix diagonal, which has the wrong units, and the values were far too small (5e-3 where the real shift was 215 sigma; the test fit shows 215, 98, 17, 0.43 instead of 5e-3 … 8e-5). Convergence criterion 1 (`pshift`, default 0.005) uses this number. A fit that stopped on criterion 1 while parameters still moved by many sigma will now keep going, and a script that reads the column must expect larger values. Fits that stop on the chi² criterion stop on the same cycle with the same parameters, esd's and chi². If a parameter has no sigma (`dp(k)` = 0), it is skipped in the maximum; before, a zero curvature diagonal was skipped.
- **Items 5 and 6 (fewer macro runs).** Results are bit-identical, but the user macro runs less often. This matters only if the macro has side effects: a counter or a global variable that is increased per run, files it writes, a log it appends to, or a plot macro that expects to run at a rejected trial point. After these items:
  - Each cycle runs the macro once less at the centre point.
  - A rejected Levenberg–Marquardt step costs one run instead of 2N+1 (3-point) or 4N+1 (5-point), so the number of runs per cycle varies.
  - The program state after a cycle is the one a repeated centre run would have left. The one exception is a per-cycle plot macro (`run ..., plot:`) and MPI, which keep the repeated centre run, because the plot reads the data in KUPLOT and the MPI master loads the derivative results into KUPLOT.
  - Derivatives that the user supplies analytically are not affected.
  - Anything that counted macro runs to find out how many cycles a fit took should count the lines of the cycle table instead.
- **Item 7 (`points:2`).** Nothing changes unless a parameter is declared with `points:2`. The default stays `points:3`. The error message for an invalid `points:` value now also names the allowed value 2. Use it for smooth parameters; the one-sided difference has an error proportional to the step, not its square, so a parameter with strong curvature or a small step on a noisy model can give worse derivatives than `points:3`.

## Not included (found during the review)

- MPI: send the centre run together with the derivative runs.
- Broyden-type Jacobian updates.
- RMC: avoid copying all amplitudes on every move.
- `seed <n>` with a single value is not reproducible (the other seed words come from the clock); the tests therefore set all 8.

## Test kit

The scripts, test inputs and per-commit handoff documents are not part of this branch. They are kept at the tag `perf-refine-handoffs-final` (`git checkout perf-refine-handoffs-final`, folder `handoffs/perf-refine/`; start with its `README.md` and `00_SETUP.md`).
