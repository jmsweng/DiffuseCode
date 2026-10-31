# Handoff 09 — check the MPI derivative fix with a real MPI build (experimental)

Commit c1 (`REFINE: fix MPI derivative step sizes being overwritten`). Requires handoff 00 (the serial binaries and the kit). Run it after 01–08; it is independent of their results.

Handoff 01 could only check the source, because the serial build has no MPI. This one builds two MPI-enabled binaries and runs a real `mpirun` refinement with each:
- `mpi_c0`: upstream code (bug present)
- `mpi_c1`: with the fix

The bug only changes 5-point derivatives under MPI, so the test fit uses `points:5` with two parameters (P_amp and P_bg, whose step sizes differ by a factor of 30). The pass criterion is: the fit under MPI with the fix reaches the same parameters as the serial fit, and the unfixed MPI build does not.

## Status of this check (read this first)
The owner's cloud session could build MPI binaries and start MPI runs, but **no MPI fit has ever completed there**: with a self-contained work macro all worker runs succeed and write their results, but the master then stops with `All derivatives are zero; check user macro`. That happens identically with `mpi_c0` and `mpi_c1`, so it is not caused by the fix. The cause is unknown (a mismatch between the test macro and what upstream's MPI refine expects, a problem of the cloud environment, or an upstream MPI problem). Because of that, **INCONCLUSIVE is a likely and acceptable outcome**. Do not try to make it pass. Report the raw output; the owner decides what to do next.

Things the cloud session did learn, which the script already handles:
1. Without a terminal the program cannot find its parent process (it runs `ps j`), silently decides MPI is off, and runs the whole fit serially. A `wsl -e bash -lc ...` call from PowerShell has no terminal. The script therefore starts `mpirun` through `script -qec`. The script checks the worker log files to make sure MPI was really on.
2. The MPI workers start with empty KUPLOT, so the work macro must create its own data (`func ..., -5, 5, 0.1`) instead of using `F_DATA`.
3. With fewer workers than derivative runs, a worker's second HDF5 write failed in the cloud (`Error opening file`). The script uses 9 ranks for the 8 derivative runs of this fit, so each worker runs once.

## Requirements
OpenMPI and `script` (package `bsdutils`):
```powershell
wsl -e bash -lc "command -v mpif90 mpirun script || echo MISSING"
```
If it prints `MISSING`, stop and tell the owner to run `sudo apt install -y libopenmpi-dev openmpi-bin` in WSL (do not install anything yourself). The build needs about 15–30 minutes and 2–4 GB of RAM; 9 MPI ranks run at the same time, so close other programs.

## Run (PowerShell, one command)
The first run builds the two binaries; if the tool call times out, run the same command again (finished binaries are kept in `~/refine_validation/mpi`).
```powershell
wsl -e bash -lc "git -C ~/refine_validation/src fetch -q origin && git -C ~/refine_validation/src checkout -q perf/refine-handoffs && git -C ~/refine_validation/src reset -q --hard origin/perf/refine-handoffs && cp ~/refine_validation/src/handoffs/perf-refine/kit/mpi_check.sh ~/refine_validation/kit/ && cd ~/refine_validation && { echo; echo '## 09 mpi'; echo; echo '~~~'; bash kit/mpi_check.sh; echo '~~~'; } 2>&1 | tee -a REPORT.md | tail -n 40"
```

## Possible results
| Last line | Meaning |
|---|---|
| `RESULT: PASS` | Fixed MPI fit equals the serial fit (and the unfixed one does not). The MPI fix is confirmed. |
| `RESULT: FAIL (...)` | Stop and report. `MPI was not active` means the run was serial, not a failure of the fix. |
| `RESULT: INCONCLUSIVE (...)` | The MPI fit did not finish (likely, see above) or the case cannot tell fixed from unfixed. Not a failure of the fix. |

## Report back
1. The full output of the command. It prints parameters, cycle counts, the first error lines and the number of worker logs for `serial`, `mpi_c0` and `mpi_c1`; on `FAIL` or `INCONCLUSIVE` it also prints the end of the three raw outputs.
2. If you had to restart the build because of a tool timeout, say so.
