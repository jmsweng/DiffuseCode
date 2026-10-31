# Handoff 09 — check the MPI derivative fix with a real MPI build

Commit c1 (`REFINE: fix MPI derivative step sizes being overwritten`), plus the KUPLOT fix `d1edc81` (`KUPLOT: save the complete data set when it is returned to the REFINE master`) that makes MPI refinements of a KUPLOT curve possible at all. Requires handoff 00 (serial binaries and the kit). Independent of 01–08.

Handoff 01 could only check the source, because the serial build has no MPI. This one builds two MPI-enabled binaries and runs a real `mpirun` refinement with each:
- `mpi_c0k`: upstream `refine_run.f90` (step-size bug present)
- `mpi_c1k`: with the c1 fix

Both contain the KUPLOT patch `kit/mpi_kuplot_exit.patch` (the same change as commit `d1edc81`). Without it, every MPI worker returns a single point and the master stops with `All derivatives are zero`. That was why the first version of this handoff (run on 2026-10-30) was inconclusive.

The test fit has 2 parameters (P_amp, P_bg), both with `points:5`, so 8 derivative runs per cycle. The pass criterion is: in cycle 0 the chi² of the fit under MPI **with** the fix equals the serial chi², and the one **without** the fix does not. (Cloud result: serial 0.50458, with fix 0.50458, without fix 244.26.)

## What to know
1. Without a terminal the program cannot find its parent process (it runs `ps j`), silently decides MPI is off, and runs the whole fit serially. A `wsl -e bash -lc ...` call from PowerShell has no terminal. The script starts `mpirun` through `script -qec` and checks that worker log files exist.
2. The MPI workers start with empty KUPLOT, so the work macro creates its own data (`func ..., -5, 5, 0.1`) instead of using `F_DATA`.
3. 9 ranks are used for the 8 derivative runs, so each worker runs once per cycle.
4. **Known, not checked here:** from cycle 1 on, `Error opening file` appears when a worker writes its result file a second time. The cause is not known. That is why the check compares cycle 0 only.

## Requirements
OpenMPI and `script` (package `bsdutils`):
```powershell
wsl -e bash -lc "command -v mpif90 mpirun script || echo MISSING"
```
If it prints `MISSING`, stop and tell the owner to run `sudo apt install -y libopenmpi-dev openmpi-bin` in WSL (do not install anything yourself). The build needs about 15–30 minutes and 2–4 GB of RAM; 9 MPI ranks run at the same time, so close other programs. If an earlier run of this handoff built `mpi_c0`/`mpi_c1` in `~/refine_validation/mpi`, they are not reused (the new binaries have the suffix `k`); the existing build folder is reused, so the build is faster.

## Run (PowerShell, one command)
If the tool call times out, run the same command again (finished binaries are kept).
```powershell
wsl -e bash -lc "git -C ~/refine_validation/src fetch -q origin && git -C ~/refine_validation/src checkout -q perf/refine-handoffs && git -C ~/refine_validation/src reset -q --hard origin/perf/refine-handoffs && cp ~/refine_validation/src/handoffs/perf-refine/kit/mpi_check.sh ~/refine_validation/src/handoffs/perf-refine/kit/mpi_kuplot_exit.patch ~/refine_validation/kit/ && cd ~/refine_validation && { echo; echo '## 09 mpi'; echo; echo '~~~'; bash kit/mpi_check.sh; echo '~~~'; } 2>&1 | tee -a REPORT.md | tail -n 40"
```

## Possible results
| Last line | Meaning |
|---|---|
| `RESULT: PASS` | Cycle 0 with the fix equals the serial fit and without the fix it does not. The MPI fix is confirmed. |
| `RESULT: FAIL (...)` | Stop and report. `MPI was not active` means the run was serial, not a failure of the fix. |
| `RESULT: INCONCLUSIVE (...)` | The MPI fit has no cycle 0, or the unfixed build also matches the serial fit. Report the output. |

## Report back
1. The full output of the command (parameters, chi² of cycle 0, first error lines and the number of worker logs for `serial`, `mpi_c0` and `mpi_c1`; on `FAIL` or `INCONCLUSIVE` also the end of the three raw outputs).
2. If you had to restart the build because of a tool timeout, say so.
