# Handoff 00 — set up the perf/refine validation (do this first)

You are validating the branch `perf/refine` of https://github.com/jmsweng/DiffuseCode.
Every later handoff (01–07) checks one commit of that branch against the commit before it.
This handoff builds one `discus_suite` binary per commit. You do **not** edit any code.

## Rules for all perf/refine handoffs
- Work only inside WSL under `~/refine_validation`. Never under `/mnt/c` (Windows line endings break the build).
- Do not edit, commit or push anything in the repository. Do not "fix" failures.
- Run the commands exactly as written, from PowerShell. They are single `wsl` calls.
- If a command fails or a check says FAIL: stop, copy the full output into the report, and end the handoff.
- Append each handoff's output to `~/refine_validation/REPORT.md` (step 3 below shows how).

## Step 1 — get the source (once)
```powershell
wsl -e bash -lc "mkdir -p ~/refine_validation && cd ~/refine_validation && (test -d src || git clone --branch perf/refine-handoffs https://github.com/jmsweng/DiffuseCode src) && git -C src fetch origin && git -C src checkout -q perf/refine-handoffs && git -C src reset -q --hard origin/perf/refine-handoffs && git -C src log -1 --oneline"
```
Expected: one line starting with the newest commit of `perf/refine-handoffs`.

## Step 2 — build the 8 binaries
This takes 10–30 minutes the first time. If your tool call times out, run the same command again: finished binaries are skipped.
```powershell
wsl -e bash -lc "bash ~/refine_validation/src/handoffs/perf-refine/kit/setup.sh 2>&1 | tail -n 12"
```
Expected last line:
```
SETUP DONE: 8 binaries in /home/<user>/refine_validation/bins
```
If it says `cmake failed` or `build of cN failed`, show the end of that log and stop:
```powershell
wsl -e bash -lc "tail -n 40 ~/refine_validation/build/cmake.log ~/refine_validation/build/make_c0.log 2>/dev/null"
```
(Known requirement from earlier DiffuseCode builds: packages `gfortran cmake pgplot5 libhdf5-dev libfftw3-dev libreadline-dev`. Do not install anything without asking the owner.)

## Step 3 — start the report
```powershell
wsl -e bash -lc "cd ~/refine_validation && { echo '# perf/refine validation report'; echo; echo '## 00 setup'; echo; uname -a; gfortran --version | head -1; cat bins/c*.sha | nl -v0; echo; } > REPORT.md && cat REPORT.md"
```
Report back: the output of step 3. Then continue with handoff 01.
