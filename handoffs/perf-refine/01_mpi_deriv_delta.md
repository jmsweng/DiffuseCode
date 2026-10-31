# Handoff 01 — REFINE: fix MPI derivative step sizes being overwritten

Branch `perf/refine`, commit c1 (see `kit/COMMITS`). Type: bug fix.
Compares binary `c0` (before) with `c1` (after). Requires handoff 00 to be done.

## What changed
In `refine_calc_deriv_mpi` (`refine/prog/refine_run.f90`) three lines wrote `delta = ...` to the whole step-size array instead of `delta(k) = ...`. Each parameter overwrote the steps of all others, so with MPI the 5-point derivative of every parameter was divided by the wrong step. 3-point derivatives and runs without MPI are not affected.

## What the check does
The owner's build has no MPI, so this check cannot run an MPI refinement. It confirms the fix is in the source, that everything compiles, and that serial refinements are unchanged.

## Run (PowerShell, one command)
```powershell
wsl -e bash -lc "cd ~/refine_validation/kit && { echo; echo '## 01 mpi_deriv_delta'; echo; echo '~~~'; python3 check.py 1; echo '~~~'; } 2>&1 | tee -a ../REPORT.md | tail -n 40"
```
It takes about 1–3 minutes. Raw outputs are kept in `~/refine_validation/runs/c0/` and `~/refine_validation/runs/c1/`.

## Expected (numbers may differ slightly with gfortran 11; PASS/FAIL is what counts)
~~~
PASS  refine_calc_deriv_mpi has no whole-array "delta =" assignment (0 found)
PASS  det     serial refinement identical to base (MPI-only change), runs 52/52
PASS  p5 ..., PASS  ran ...
RESULT: PASS
~~~

## Report back
The full output of the command. If the last line is `RESULT: PASS`, continue with the next handoff.
If it is `RESULT: FAIL`, stop. Add the output of this command to your report and end the handoff (do not try to fix anything):
```powershell
wsl -e bash -lc "bash ~/refine_validation/kit/dump.sh 1"
```
