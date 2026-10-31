# Handoff 05 — REFINE: skip the repeated centre run after the derivatives

Branch `perf/refine`, commit c5 (see `kit/COMMITS`). Type: performance.
Compares binary `c4` (before) with `c5` (after). Requires handoff 00 to be done.

## What changed
After the finite-difference runs, `refine_theory` ran the user macro a second time at the centre point. The centre values were already stored, so that run is now skipped unless a plot macro is given (`run ..., plot:`) or MPI is active. The random generator is put back exactly where the repeated run would have left it, so results are bit-identical. Saves one macro run per cycle.

## What the check does
Every fit must give identical cycle tables, parameters and esd's with fewer macro runs. With a plot macro only the final (unneeded) repeat is saved.

## Run (PowerShell, one command)
```powershell
wsl -e bash -lc "cd ~/refine_validation/kit && { echo; echo '## 05 skip_centre_rerun'; echo; echo '~~~'; python3 check.py 5; echo '~~~'; } 2>&1 | tee -a ../REPORT.md | tail -n 40"
```
It takes about 1–3 minutes. Raw outputs are kept in `~/refine_validation/runs/c4/` and `~/refine_validation/runs/c5/`.

## Expected (numbers may differ slightly with gfortran 11; PASS/FAIL is what counts)
~~~
PASS  det     identical results, macro runs 52 -> 46 (similar lines for every case)
PASS  plot    with plot: only the final repeat is saved (52 -> 51), plots 4 -> 4
RESULT: PASS
~~~

## Report back
The full output of the command. If the last line is `RESULT: PASS`, continue with the next handoff.
If it is `RESULT: FAIL`, stop. Add the output of this command to your report and end the handoff (do not try to fix anything):
```powershell
wsl -e bash -lc "bash ~/refine_validation/kit/dump.sh 5"
```
