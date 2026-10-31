# Handoff 03 — RMC: fix refinement with lots

Branch `perf/refine`, commit c3 (see `kit/COMMITS`). Type: bug fix.
Compares binary `c2` (before) with `c3` (after). Requires handoff 00 to be done.

## What changed
RMC with `set lots` was broken in four places in `discus/prog/rmc_sup_mod.f90`: the lot test used the Fourier setting `ilots` instead of `rmc_ilots` (so no move was ever accepted), its periodic wrap could never succeed, every move was added to lot 1, lot 1 was counted twice, and the average amplitude was subtracted over the wrong number of rows.

## What the check does
The test runs 300 accepted RMC moves (amplitudes updated incrementally), saves the structure, recalculates everything from scratch, and compares the two chi^2 values. They must agree to about 1e-7 (the structure file rounds positions). `fourlots` also sets `fourier lots`, the old workaround that let moves happen.

## Run (PowerShell, one command)
```powershell
wsl -e bash -lc "cd ~/refine_validation/kit && { echo; echo '## 03 rmc_lots'; echo; echo '~~~'; python3 check.py 3; echo '~~~'; } 2>&1 | tee -a ../REPORT.md | tail -n 40"
```
It takes about 1–3 minutes. Raw outputs are kept in `~/refine_validation/runs/c2/` and `~/refine_validation/runs/c3/`.

## Expected (numbers may differ slightly with gfortran 11; PASS/FAIL is what counts)
~~~
INFO  nolots/lots/fourlots before and after values
PASS  nolots   after: incremental and fresh chi^2 agree
PASS  nolots   before: also consistent (control)
PASS  lots     after: ... agree
PASS  fourlots after: ... agree
(before: 'lots' shows 'no valid move' / None, 'fourlots' differs by ~5e-3)
RESULT: PASS
~~~

## Report back
The full output of the command. If the last line is `RESULT: PASS`, continue with the next handoff.
If it is `RESULT: FAIL`, stop. Add the output of this command to your report and end the handoff (do not try to fix anything):
```powershell
wsl -e bash -lc "bash ~/refine_validation/kit/dump.sh 3"
```
