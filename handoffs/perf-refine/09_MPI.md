# Handoff 09 — check REFINE under MPI (serial against `mpirun`)

Checks three commits of `perf/refine` with a real MPI build:
- `2e95b93` MPI step sizes (the 5-point derivative used the wrong step)
- `d1edc81` KUPLOT: a worker returned only one point of its curve (`All derivatives are zero`)
- `0d73677` HDF5 write: a worker's second result write failed with a false `Error opening file`

and, because the build is from the tip of `perf/refine`, also that the speed-ups (skipped centre run, derivatives only for accepted steps, `points:2`) give the same result under MPI as in a serial run.

Requires handoff 00 (it uses `bins/c7`, the serial binary, as the reference). Independent of 01–08.

## What the check does
It builds two MPI-enabled `discus_suite` binaries in `~/refine_validation/mpi`: one from the tip of `origin/perf/refine`, and one that has the old step-size bug put back (`kit/mpi_delta_unfix.patch`). Then it runs four fits of a 4-parameter Gaussian, serially and with `mpirun -n 5` (4 workers, so each worker runs several derivative jobs per cycle):

| case | what it tests |
|---|---|
| `det` | easy start, 3-point derivatives |
| `far` | far start with rejected steps (derivatives only for accepted steps) |
| `p5` | 5-point derivatives (16 derivative runs per cycle) |
| `p2far` | `points:2`, far start |

PASS for a case: same number of cycles, same parameters (relative 1e-4), no error lines, and worker log files present (so MPI really was on). In addition `p5` is run with the old-bug binary; its chi² in cycle 0 must differ from the serial one (cloud: 365 against 109), otherwise the case cannot detect the bug and the result is `INCONCLUSIVE`.

Cloud result (OpenMPI 4.1.6, gfortran 13.3): all four cases PASS, bug detected.

## What to know
1. Without a terminal the program cannot find its parent process (it runs `ps j`), silently decides MPI is off, and runs the whole fit serially. A `wsl -e bash -lc ...` call from PowerShell has no terminal. The script starts `mpirun` through `script -qec` and checks the worker log files.
2. MPI workers start with empty KUPLOT, so the work macro creates its own data (`func ..., -5, 5, 0.1`) instead of using `F_DATA`.

## Requirements
OpenMPI and `script` (package `bsdutils`):
```powershell
wsl -e bash -lc "command -v mpif90 mpirun script || echo MISSING"
```
If it prints `MISSING`, stop and tell the owner to run `sudo apt install -y libopenmpi-dev openmpi-bin` in WSL (do not install anything yourself). The build needs about 15–30 minutes and 2–4 GB of RAM. The fits start 5 MPI processes at once; close other programs.

## Run (PowerShell, one command)
If the tool call times out, run the same command again (finished binaries are kept).
```powershell
wsl -e bash -lc "git -C ~/refine_validation/src fetch -q origin && git -C ~/refine_validation/src checkout -q perf/refine-handoffs && git -C ~/refine_validation/src reset -q --hard origin/perf/refine-handoffs && cp ~/refine_validation/src/handoffs/perf-refine/kit/mpi_check.sh ~/refine_validation/src/handoffs/perf-refine/kit/mpi_delta_unfix.patch ~/refine_validation/kit/ && cd ~/refine_validation && { echo; echo '## 09 mpi'; echo; echo '~~~'; bash kit/mpi_check.sh; echo '~~~'; } 2>&1 | tee -a REPORT.md | tail -n 40"
```
An earlier build in `~/refine_validation/mpi` is reused as the build folder, but its binaries are not (the new ones are `mpi_tip` and `mpi_tip_unfixed`).

## Possible results
| Last line | Meaning |
|---|---|
| `RESULT: PASS` | MPI and serial fits agree in all four cases and the old bug is detected. |
| `RESULT: FAIL` | A case differs, has errors, or MPI was not active (`worker logs 0`). Stop and report; the script prints the end of the MPI outputs. |
| `RESULT: INCONCLUSIVE` | All cases agree, but the old bug is not detected. Report the output. |

## Report back
1. The full output of the command.
2. If you had to restart the build because of a tool timeout, say so.
