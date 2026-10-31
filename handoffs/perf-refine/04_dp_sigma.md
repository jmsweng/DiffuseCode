# Handoff 04 — REFINE: compute the dP/sigma convergence measure from the parameter sigmas

Branch `perf/refine`, commit c4 (see `kit/COMMITS`). Type: bug fix.
Compares binary `c3` (before) with `c4` (after). Requires handoff 00 to be done.

## What changed
`refine_mrq` divided each parameter shift by `SQRT(cl(k,k))`, but at that point `cl` holds the curvature matrix, not the covariance. The printed `MAX(dP/sig)` was therefore meaningless (about 5e-3 when parameters moved by 200 sigma), and convergence criterion 1 (`pshift`) could stop a fit too early. It now divides by the parameter sigma `dp(k)`.

## What the check does
In the test fits only the dP/sig column changes; cycles, chi^2, parameters and esd's must be identical. Real fits that used to stop on criterion 1 may now run a few more cycles: that is the intended correction.

## Run (PowerShell, one command)
```powershell
wsl -e bash -lc "cd ~/refine_validation/kit && { echo; echo '## 04 dp_sigma'; echo; echo '~~~'; python3 check.py 4; echo '~~~'; } 2>&1 | tee -a ../REPORT.md | tail -n 40"
```
It takes about 1–3 minutes. Raw outputs are kept in `~/refine_validation/runs/c3/` and `~/refine_validation/runs/c4/`.

## Expected (numbers may differ slightly with gfortran 11; PASS/FAIL is what counts)
~~~
PASS  <each case> cycles, chi^2, parameters, esd unchanged apart from the dP/sig column
INFO  det dP/sig before: [0.0050712, 0.018693, ...]
INFO  det dP/sig after : [214.98, 97.938, ...]
PASS  det     first-cycle dP/sig is now a multiple of sigma (> 10)
RESULT: PASS
~~~

## Report back
The full output of the command. If the last line is `RESULT: PASS`, continue with the next handoff.
If it is `RESULT: FAIL`, stop. Add the output of this command to your report and end the handoff (do not try to fix anything):
```powershell
wsl -e bash -lc "bash ~/refine_validation/kit/dump.sh 4"
```
