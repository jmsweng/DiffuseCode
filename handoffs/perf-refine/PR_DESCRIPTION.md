## REFINE / RMC / KUPLOT / HDF5: six bug fixes and three speed-ups for structure refinement

Nine small, independent commits on top of `master` (6fd5f33). Each commit message holds the details and its own test numbers; this page summarises them. The fixes come first, so any prefix of the series can be merged. Commits 5 → 6 → 7 build on each other.

### Bug fixes

| # | Commit | Problem | Effect of the fix |
|---|---|---|---|
| 1 | `REFINE: fix MPI derivative step sizes being overwritten` | In `refine_calc_deriv_mpi` three assignments wrote `delta = ...` (the whole array) instead of `delta(k) = ...`. The 5-point derivative divides by `delta(k)` after all children return and therefore used the last parameter's step. 3-point derivatives were not affected. | 5-point derivatives under MPI use each parameter's own step. `points:5` fit under MPI, cycle 0: chi²/(N-P) 109.18 as in the serial fit, 365.29 without this fix (commits 8 and 9 applied in both builds). |
| 2 | `Make seed restore reproduce Gaussian random numbers` | `put_seeds` restored the `RANDOM_NUMBER` state but not the cached second value of the Box-Muller generator `GASDEV` (`iset`/`gset`). After a restore the first Gaussian number depended on whatever ran before. REFINE restores the seeds so that the centre and all derivative points see the same random numbers, so any model with Gaussian randomness (thermal displacements, MC, `gran()`) got noisy derivatives. | A fit with `gran(0.02)` noise: 21 cycles and a wrong minimum (chi² 7992) before, 4 cycles and the true parameters (chi² 76) after. Deterministic fits are unchanged. |
| 3 | `RMC: fix refinement with lots` | `rmc_inlot` tested the Fourier lot type `ilots` instead of `rmc_ilots` (no move was ever made unless `fourier lots` was set), its periodic wrap could never trigger, `rmc_strucf` added every move to lot 1, `rmc_inten` counted lot 1 twice, and `rmc_calc_scat` subtracted the average over `num(1) x num(1)` points instead of `num(1) x num(2)`. | Cu 6x6x6, one hk0 plane, 300 accepted moves: incremental and freshly calculated chi² now agree to 8e-8 (they differed by 0.5 % with `fourier lots`; without `fourier lots` no move was possible). Without lots nothing changes. |
| 4 | `REFINE: compute the dP/sigma convergence measure from the parameter sigmas` | `refine_mrq` divided the parameter shift by `SQRT(cl(k,k))`, but after `mrqmin` `cl` holds the curvature matrix of the new point, not its inverse. The printed `MAX(dP/sig)` therefore had the wrong units (e.g. 5e-3 where the real shift was 215 sigma), and convergence criterion 1 (`pshift`, default 0.005) could be met while parameters still moved by many sigma. | The column now divides by the parameter sigma. See "Behaviour changes". |
| 8 | `KUPLOT: save the complete data set when it is returned to the REFINE master` | With MPI each worker returns its last KUPLOT data set as `DISCUS_SUITE_DERIVATIVES/data.nnnn`. `kuplot_exit` called `check_form` with no limits, which for the format H5 leaves the limits at zero; `do_save` then writes only the range 0..0, i.e. the single point at x = 0, for any data set not read from an HDF5 file (a curve from `func`/`load xy`, a NIPL map). The master loaded a 1x1x1 data set, KUPLOT did not register it, and every derivative run was read as the original data: `All derivatives are zero`. | `kuplot_exit` passes the full x (and y) range for such data sets. Before, no MPI refinement of a KUPLOT curve finished. |
| 9 | `HDF5 write: initialise hdferr before testing it` | The two write routines in `hdf_write.f90` open the HDF5 library only the first time (`hdf5_get_state`) and then test `hdferr`. When the library is already open, `H5open_f` is skipped and `hdferr` is an uninitialised local variable, so the test could report `Error opening file` (`ier_num = -2`) although nothing went wrong. With MPI each worker writes its result file at the end of every derivative run, so from the second cycle on the runs failed. | `hdferr = 0` before the test. Multi-cycle MPI refinements work. |

### Speed-ups

All cost figures are the number of times the user macro runs, which is what dominates a real refinement (each run is a full DISCUS simulation).

| # | Commit | Change | Effect |
|---|---|---|---|
| 5 | `REFINE: skip the repeated centre run after the derivatives` | After the derivatives, `refine_theory` ran the macro a second time at the centre point to leave the program state at the current parameters. That only matters if something reads the state before the next run. It is now done only with a per-cycle plot macro (`run ..., plot:`) or under MPI; otherwise the random generator is restored to the state the repeated run would have left. Also fixes a read past the array bound for a single parameter (`gl_is_der(1) .and. gl_is_der(2)`). | −1 run per cycle, bit-identical results. |
| 6 | `REFINE: compute derivatives only for accepted steps` | `mrqmin` computed all derivatives at every trial point and discarded them when chi² did not improve. The trial point is now evaluated with one run; only an accepted step computes the derivatives, reusing the trial run as their centre. | A rejected step costs 1 run instead of 2N+1 (3-point) or 4N+1 (5-point). Bit-identical results. |
| 7 | `REFINE: add two-point (one-sided) derivatives, 'points:2'` | New opt-in `points:2` for `newparam`: `(f(p+h) - f(p))/h` needs one extra run per parameter instead of two. The default stays `points:3`; help text and the error message for invalid `points:` values are updated. | N instead of 2N runs per Jacobian. Suited to smooth parameters such as scale factors (error ∝ h instead of h²). |

### Measured effect

Number of user-macro runs for a 4-parameter Gaussian fit (fixed seeds, KUPLOT-only model, so it runs in seconds). Base is `master`; the other columns are the cumulative state after that commit.

| case | base | after 4 | after 5 | after 6 | after 7 |
|---|---|---|---|---|---|
| easy start | 52 | 52 | 46 | 46 | 46 |
| far start (6 of 15 steps rejected) | 162 | 162 | 145 | 97 | 97 |
| 25 fixed cycles (13 rejected) | 262 | 262 | 235 | 131 | 131 |
| with `gran()` noise | 222 | 52 | 46 | 38 | 38 |
| noise, 25 fixed cycles | 262 | 262 | 235 | 91 | 91 |
| 5-point derivatives | 92 | 92 | 86 | 86 | 86 |
| `points:2` (easy / far / noise) | error | error | error | error | 26 / 63 / 22 |

Parameters, esd's and chi² of the deterministic cases are the same in every column; only the number of runs changes. The "noise" rows change with commit 2 because the fit now converges.

### Behaviour changes to be aware of

- **Commit 4** changes the printed `MAX(dP/sig)` column (215, 98, 17, 0.43 sigma in the test fit instead of 5e-3 … 8e-5). Stopping cycle, parameters and esd's of the test fits are unchanged because they stop on the chi² criterion, but a fit that used to stop on criterion 1 may now run more cycles (correctly).
- **Commit 2** changes the results of any refinement or MC calculation that draws Gaussian random numbers after a seed restore. This is the point of the fix, but numbers will differ from earlier runs.
- **Commit 5** is bit-identical in the test fits; the restore of the random state after the last derivative run replaces the repeated centre run.
- Everything else is bit-identical (cycle tables, parameters, esd's) to the previous commit; commits 5 and 6 were also checked with a noisy model and fixed seeds.

### How it was tested

A set of test fits (a Gaussian peak with 4 parameters, in KUPLOT only, plus an RMC test on Cu 6x6x6) compares the binary built before each commit with the one built after it: cycle tables, parameters, esd's, chi² and macro-run counts. The kit is on the author's fork, `jmsweng/DiffuseCode`, branch `perf/refine-handoffs`, under `handoffs/perf-refine/` (not part of this PR). It was run twice, with the same results:

- cloud, Ubuntu, gfortran 13.3
- the author's laptop, Ubuntu 22.04 under WSL2, gfortran 11.4.0

Both runs: all checks PASS for commits 1–7 (the numbers in the table above come from the first run; the second gave identical run counts). The MPI check (commits 1, 8 and 9) has so far been run in the cloud only.

### What is not tested

- **MPI was tested only with a KUPLOT model and with OpenMPI 4.1.6 / gfortran 13.3, in one environment so far.** Four 4-parameter fits (easy start, far start with rejected steps, `points:5`, `points:2` far start) were run serially and with `mpirun -n 5` (4 workers): same number of cycles, same parameters and esd's (the last digit of an esd can differ), no errors. A build with the old step-size bug gives chi²/(N-P) 365 against 109 in cycle 0 of the `points:5` fit. Parallel REFINE with DISCUS models (3D data) was not tested.
- Only `points:2/3/5` with one data set were tested. Multi-data-set REFINE (`mrqcof` loops over data set 1's size for all sets) was left alone.
- The RMC fix was tested with box lots and with box lots plus `fourier lots`, on one hk0 plane; ellipsoidal lots were not tested.

### Not included (found during the review)

- MPI: send the centre run together with the derivative runs (needs an MPI build to test).
- Broyden-type Jacobian updates.
- RMC: avoid copying all amplitudes on every move.
- `seed <n>` with a single value is not reproducible (the other seed words come from the clock); the tests therefore set all 8.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

https://claude.ai/code/session_018swRnCiTTBsHJGR5Q1FsFt
