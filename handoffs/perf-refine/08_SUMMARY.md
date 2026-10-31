# Handoff 08 — summary table and final report

Requires handoffs 00–07 to be done (all `RESULT: PASS`). Runs every test fit with all 8 binaries and prints how many times the user macro ran.

## Run (PowerShell, one command, about 5 minutes)
```powershell
wsl -e bash -lc "cd ~/refine_validation/kit && { echo; echo '## 08 summary'; echo; echo '~~~'; python3 check.py all; echo '~~~'; } 2>&1 | tee -a ../REPORT.md | tail -n 20"
```

## Expected
A table like the one in `README.md`: every column from c5 on has fewer runs than c0, `ran` drops from about 222 to 38, and the `p2*` rows are `err` up to c6 and numbers in c7.

## Report back
1. The table.
2. Copy the report to Windows so the owner can open it:
```powershell
wsl -e bash -lc "mkdir -p /mnt/c/Users/jmswe/refine_validation && cp ~/refine_validation/REPORT.md /mnt/c/Users/jmswe/refine_validation/REPORT.md && echo copied"
```
3. One line per handoff 01–07: PASS or FAIL.
