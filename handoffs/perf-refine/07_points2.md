# Handoff 07 — REFINE: add two-point (one-sided) derivatives, 'points:2'

Branch `perf/refine`, commit c7 (see `kit/COMMITS`). Type: performance.
Compares binary `c6` (before) with `c7` (after). Requires handoff 00 to be done.

## What changed
New option `newparam <name>, points:2`: the derivative is (f(p+h) - f(p))/h, one macro run per parameter instead of two (points:3) or four (points:5). Near the upper end of a range p-h is used. The default stays points:3; the help text explains the accuracy trade-off.

## What the check does
All existing fits must be unchanged. The `p2*` fits must be refused by the previous binary and, with the new one, reach the same minimum as their points:3 counterparts with fewer runs.

## Run (PowerShell, one command)
```powershell
wsl -e bash -lc "cd ~/refine_validation/kit && { echo; echo '## 07 points2'; echo; echo '~~~'; python3 check.py 7; echo '~~~'; } 2>&1 | tee -a ../REPORT.md | tail -n 40"
```
It takes about 1–3 minutes. Raw outputs are kept in `~/refine_validation/runs/c6/` and `~/refine_validation/runs/c7/`.

## Expected (numbers may differ slightly with gfortran 11; PASS/FAIL is what counts)
~~~
PASS  <each old case> points:3/5 fit unchanged
PASS  p2      before: points:2 is rejected with an error (and p2far, p2range, p2ran)
PASS  p2      converges to the points:3 result within ~0.001 esd, macro runs 26 (points:3: 46)
PASS  p2far ... 63 (97), p2range ... 21 (37), p2ran ... 22 (38)
RESULT: PASS
~~~

## Report back
The full output of the command. If the last line is `RESULT: PASS`, continue with the next handoff.
If it is `RESULT: FAIL`, stop. Add the output of this command to your report and end the handoff (do not try to fix anything):
```powershell
wsl -e bash -lc "bash ~/refine_validation/kit/dump.sh 7"
```
