# perf/refine — faster and correct REFINE/RMC refinements

Branch `perf/refine` of jmsweng/DiffuseCode starts at upstream tproffen/DiffuseCode `master` (6fd5f33, 2026-10-28) and holds one commit per item, fixes first. It contains only source changes, so it can be offered upstream as it is, or item by item (each commit is self-contained; 05→06→07 build on each other).

This folder lives only on `perf/refine-handoffs` (= `perf/refine` + this folder) and is not meant for upstream.

| Handoff | Commit | Type | Change | Effect in the test fits |
|---|---|---|---|---|
| [00](00_SETUP.md) | — | setup | Build one binary per commit (c0 = upstream base … c7) | — |
| [01](01_mpi_deriv_delta.md) | c1 | fix | MPI: 5-point derivatives used another parameter's step size | compile/source check; serial fits unchanged |
| [02](02_gaussian_seed_restore.md) | c2 | fix | Seed restore did not reset the cached Gaussian random number → noisy derivatives for any model with Gaussian randomness | noisy fit: 21 cycles to a wrong minimum (chi² 7992) → 4 cycles to the true one (chi² 76) |
| [03](03_rmc_lots.md) | c3 | fix | RMC with `set lots`: no move possible, wrong lot updated, lot 1 counted twice, wrap never worked | incremental vs fresh chi² 0.45 % → 8e-8 |
| [04](04_dp_sigma.md) | c4 | fix | Printed/used `MAX(dP/sig)` was computed from the curvature matrix, not the sigmas | 0.005 → 215 sigma in cycle 0; may stop fits later (correctly) |
| [05](05_skip_centre_rerun.md) | c5 | perf | Skip the second macro run at the centre point each cycle | −1 run per cycle, bit-identical |
| [06](06_lazy_jacobian.md) | c6 | perf | Derivatives only for accepted LM steps | rejected step 2N+1 → 1 run, bit-identical |
| [07](07_points2.md) | c7 | perf | New `points:2` one-sided derivatives | N instead of 2N runs per Jacobian (opt-in) |
| [08](08_SUMMARY.md) | — | report | Macro-run table c0…c7 | — |

Macro runs (user-macro executions) per fit, measured in the cloud (gfortran 13.3), 4-parameter Gaussian fits, fixed seeds:

| case | c0 base | c4 (fixes) | c5 | c6 | c7 |
|---|---|---|---|---|---|
| det (easy start) | 52 | 52 | 46 | 46 | 46 |
| far (6 of 15 steps rejected) | 162 | 162 | 145 | 97 | 97 |
| farfix (25 cycles, 13 rejected) | 262 | 262 | 235 | 131 | 131 |
| ran (gran() noise) | 222 | 52 | 46 | 38 | 38 |
| ranfix (25 cycles, noise) | 262 | 262 | 235 | 91 | 91 |
| p5 (5-point) | 92 | 92 | 86 | 86 | 86 |
| p2 / p2far / p2ran (points:2) | error | error | error | error | 26 / 63 / 22 |

In a real refinement each macro run is a full DISCUS simulation, so wall time scales with these counts.

## For the local agent
Do the handoffs in order 00 → 08. Each one is a single PowerShell command (a `wsl` call) that ends with `RESULT: PASS` or `RESULT: FAIL`. Never edit code; on FAIL, report and stop.

## Kit
`kit/setup.sh` builds the binaries listed in `kit/COMMITS`; `kit/check.py N` runs the checks for commit cN against c(N-1); `kit/check.py all` prints the run table; `kit/dump.sh N` condenses raw outputs for a FAIL report. Test inputs are in `kit/refine/` (REFINE, KUPLOT-only model, runs in seconds) and `kit/rmc/` (RMC, Cu 6×6×6).

## Not in this branch (from the review)
- MPI: send the centre run together with the derivative runs (needs an MPI build to test).
- Broyden-type Jacobian updates.
- Multi-data-set REFINE (`mrqcof` loops over data set 1's size for all sets) — upstream work in progress.
- RMC: avoid copying all amplitudes on every move.
- `seed <n>` with a single value is not reproducible (other seed words come from the clock); the kit therefore sets all 8 values.
