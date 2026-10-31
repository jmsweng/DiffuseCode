# Handoff 06 — REFINE: compute derivatives only for accepted steps

Branch `perf/refine`, commit c6 (see `kit/COMMITS`). Type: performance.
Compares binary `c5` (before) with `c6` (after). Requires handoff 00 to be done.

## What changed
The Levenberg-Marquardt step computed all derivatives (2N or 4N macro runs) at every trial point, and threw them away when chi^2 did not improve. The trial point is now evaluated with one macro run; derivatives are computed only if the step is accepted, reusing that run as their centre. A rejected step costs 1 run instead of 2N+1.

## What the check does
Results must be identical. Fits with rejected steps (`far`, `farfix`, `ranfix`) need far fewer runs; fits without rejections stay the same.

## Run (PowerShell, one command)
```powershell
wsl -e bash -lc "cd ~/refine_validation/kit && { echo; echo '## 06 lazy_jacobian'; echo; echo '~~~'; python3 check.py 6; echo '~~~'; } 2>&1 | tee -a ../REPORT.md | tail -n 40"
```
It takes about 1–3 minutes. Raw outputs are kept in `~/refine_validation/runs/c5/` and `~/refine_validation/runs/c6/`.

## Expected (numbers may differ slightly with gfortran 11; PASS/FAIL is what counts)
~~~
PASS  far     identical results, macro runs 145 -> 97
PASS  farfix  identical results, macro runs 235 -> 131
PASS  ranfix  identical results, macro runs 235 -> 91
(det, p5, one, range: equal runs)
RESULT: PASS
~~~

## Report back
The full output of the command. If the last line is `RESULT: PASS`, continue with the next handoff.
If it is `RESULT: FAIL`, stop. Add the output of this command to your report and end the handoff (do not try to fix anything):
```powershell
wsl -e bash -lc "bash ~/refine_validation/kit/dump.sh 6"
```
