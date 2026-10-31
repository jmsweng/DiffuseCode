# Handoff 02 — Make seed restore reproduce Gaussian random numbers

Branch `perf/refine`, commit c2 (see `kit/COMMITS`). Type: bug fix.
Compares binary `c1` (before) with `c2` (after). Requires handoff 00 to be done.

## What changed
REFINE restores the random seeds before each finite-difference run so all runs see the same random numbers. The Gaussian generator `gasdev` keeps a cached second value that `put_seeds` did not reset, so after a restore the Gaussian sequence was shifted. Any model with Gaussian random numbers (thermal displacements, MC, KUPLOT `gran()`) got noisy derivatives. `put_seeds` now resets the cache, and the centre run starts from the same normalised state.

## What the check does
The test fits a Gaussian peak plus `gran(0.02)` noise. Before the fix the fit wanders off (wrong background, very large chi^2) or needs many cycles; after it, it reaches the true parameters in a few cycles. Deterministic fits must be unchanged.

## Run (PowerShell, one command)
```powershell
wsl -e bash -lc "cd ~/refine_validation/kit && { echo; echo '## 02 gaussian_seed_restore'; echo; echo '~~~'; python3 check.py 2; echo '~~~'; } 2>&1 | tee -a ../REPORT.md | tail -n 40"
```
It takes about 1–3 minutes. Raw outputs are kept in `~/refine_validation/runs/c1/` and `~/refine_validation/runs/c2/`.

## Expected (numbers may differ slightly with gfortran 11; PASS/FAIL is what counts)
~~~
PASS  det / far / p5 / range  deterministic fit unchanged
INFO  ran before: ... (expect large chi2 and/or wrong P_bg)
INFO  ran after : cycles 4 chi2 ~60-80 pars close to 8, 0.7, 1.3, 0.5
PASS  ran     noisy fit reaches the true parameters
PASS  ran     chi^2 < 150
RESULT: PASS
~~~

## Report back
The full output of the command. If the last line is `RESULT: PASS`, continue with the next handoff.
If it is `RESULT: FAIL`, stop. Add the output of this command to your report and end the handoff (do not try to fix anything):
```powershell
wsl -e bash -lc "bash ~/refine_validation/kit/dump.sh 2"
```
