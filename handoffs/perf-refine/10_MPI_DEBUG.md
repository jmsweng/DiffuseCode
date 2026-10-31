# Handoff 10 — why does the MPI master see zero derivatives? (information gathering)

Follow-up to handoff 09, which ended `INCONCLUSIVE`: with MPI, both the unfixed and the fixed build stop with `All derivatives are zero; check user macro`. Requires 09 to have run (it creates the MPI build in `~/refine_validation/mpi`). This handoff **does not fix anything and does not change the repository**. It only gathers evidence.

## What is already known (from the cloud session)
- The workers do their job: 8 worker log files, each with a different `Range y`, and 8 different result files `DISCUS_SUITE_DERIVATIVES/data.nnnn` (valid HDF5, different content).
- The master then loads each result file with `do_load('h5, DISCUS_SUITE_DERIVATIVES/data.nnnn')` and reads "the last KUPLOT data set" (`refine_load_calc`). A debug print shows that **after the load KUPLOT has not gained a data set** (`iz` stays 1) and the master reads the original data again for every derivative run. All derivative runs therefore look identical, the derivatives are zero and the fit stops. The load reports no error (`ier_num = 0`).
- So the problem sits between "worker writes the file" and "master has it in KUPLOT". Why the load does nothing is not known yet.

This handoff checks whether your machine shows exactly the same thing.

## The rules for this handoff
- You may build a debug binary. The script does it in the scratch folder `~/refine_validation/mpi`, applying `kit/mpi_debug.patch` to the scratch copy of the c1 source and restoring it afterwards.
- Do not edit, commit or push anything in `~/refine_validation/src` (the repository) and do not try to fix the problem.
- If a command fails, stop and report the output.

## Run (PowerShell, one command)
The first run recompiles one source file and relinks (a few minutes). Update the kit first so that the script and the patch are present:
```powershell
wsl -e bash -lc "git -C ~/refine_validation/src fetch -q origin && git -C ~/refine_validation/src checkout -q perf/refine-handoffs && git -C ~/refine_validation/src reset -q --hard origin/perf/refine-handoffs && cp ~/refine_validation/src/handoffs/perf-refine/kit/mpi_debug.sh ~/refine_validation/src/handoffs/perf-refine/kit/mpi_debug.patch ~/refine_validation/kit/ && cd ~/refine_validation && { echo; echo '## 10 mpi debug'; echo; echo '~~~'; bash kit/mpi_debug.sh; echo '~~~'; } 2>&1 | tee -a REPORT.md | tail -n 90"
```
It ends with `RESULT: DONE (information only)`. `RESULT: FAIL` means a build or a prerequisite problem; report it and stop.

## What the output contains
- **A.** `DBG` lines. For each of the 8 derivative runs: `DBG after load ... iz=... ndims=... lenc=... y=...` (state of KUPLOT right after the master loaded the result file) and `DBG kid=... sum/max(temp)=... sum(calc)=...` (what the master took from it).
- **B.** The first error lines.
- **C.** The worker result files: sizes, and how many distinct contents there are (expected: 8 distinct).
- **D.** The end of two worker logs.

In the cloud session the output was: `iz=1`, `ndims=0`, `lenc=101`, `y= 4.5117E-01 5.4711E-01 4.6892E-01` (the original data) in every line, identical `sum/max(temp)` for all eight runs, 8 distinct result files.

## Report back
1. The full output of the command.
2. One line saying whether A matches the cloud description (iz stays 1 and the arrays are identical for all eight runs) or differs, and how.
