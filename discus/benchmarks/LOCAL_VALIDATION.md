# Local validation of branch `perf/gemm-fourier`

Instructions for an agent running on the owner's own machine. The cloud session
that wrote the changes could not reach the owner's data, compilers or hardware;
this file lists what still has to be checked there and how.

**Repository:** https://github.com/jmsweng/DiffuseCode (fork of tproffen/DiffuseCode)
**Branch under test:** `perf/gemm-fourier`. **Reference:** `master` of the same fork.

## Ground rules

- Work in a new directory, for example `~/discus_validation/`. Never modify the
  owner's existing macros, structures or data files: copy them into the work directory first.
- Don't push to any branch or open pull requests. If you find a bug, report it,
  with a proposed patch if you have one, and leave the decision to the owner.
- Ask the owner before anything that runs for more than about an hour, or uses all cores on a shared machine.
- Record exact commands, versions and timings: the owner will compare your report with the cloud results below.

## 1. What the branch changes

| # | Change | File(s) | Active |
|---|---|---|---|
| 1 | ZGEMM Fourier kernel, exact phases instead of the 2^16 sin/cos table | `four_strucf_gemm_mod.f90`, `four_strucf_OMP_mod.f90` | Only with `cmake -DDISCUS_BLAS=ON`; layers and volumes, not lines. `DISCUS_FOURIER_KERNEL=table` switches it off at run time |
| 2 | Serial kernel uses a per-plane buffer (cache-friendly) | `four_strucf_OMP_mod.f90` | Always, when OpenMP is off or 1 thread |
| 3 | Debye-Waller table cache fixed (wrong size check, always recomputed, missing cache keys) | `fourier_sup.f90`, `diffuse_mod.f90`, `discus_allocate_appl_mod.f90` | Always (`temp use`) |
| 4 | 3D-PDF FFT: real-to-complex, multithreaded FFTW | `graphic.f90`, top `CMakeLists.txt` | Always for 3D `value 3DPDF`; threads only if FFTW's OpenMP library is found |
| 5 | NUFFT driver rewritten: exact coordinate scaling, any grid | `fourier_sup.f90`, `fourier_finufft.f90`, `appl_dis.hlp` | Always for `set technique:nufft` (not the DISCAMB variant) |

Change 5 fixes a correctness bug. The old NUFFT gave wrong intensities whenever
grid step × unit cells along an axis was not a whole number. For example, step
0.2 on 16 cells put the strongest peak at (−0.8, −0.8, −4) instead of (0, 0, −4).
The error check in master could never trigger. Results from such runs are expected to change.

## 2. What is already validated (cloud, Linux, gfortran 13.3, 2 threads)

`discus/benchmarks/regression/regress.py`: 43 cases, 82 pass/fail checks + 9
informational, comparing every point of full 3D volumes in double precision.
All pass for both the default build and `-DDISCUS_BLAS=ON`:

- With the original kernel selected, the branch reproduces master to ≤1e-12 of
  the maximum (mostly bit-identical) in all of these cases:
  - Grids: odd/even/unequal/off-centre/skewed 3D grids, layers, lines.
  - Radiation: X-ray, neutron, anomalous dispersion.
  - ADPs: iso/aniso with one and two ADP types.
  - Averaging: `set aver`; deterministic `lots` (box/ellipsoid, periodic or not).
  - Crystals: a thermally displaced crystal.
  - Outputs: 3D-PDF and 3D-ΔPDF; the serial kernel.
  - Cache: a 7-step session that changes grid, temperature switch and ADPs, against master and against fresh sessions.
- GEMM agrees with the independent NUFFT to between 1e-8 and 3e-6 of the maximum
  on every grid type. The original table kernel deviates from both:
  - by up to 1e-3 of the maximum in layers;
  - by 2e-3 with `set aver`;
  - by ~1% at single weak diffuse points (3D-ΔPDF: 3e-4 of the maximum).

  So switching to GEMM changes results by about that much. GEMM is the more accurate one.
- New NUFFT agrees with turbo to ≤3e-6 of the maximum on every grid type, and with
  old NUFFT exactly where the old one was valid.

Not covered in the cloud, and the job of this local run:
- the owner's own structures, macros and workflows;
- molecules and vacancies;
- electron radiation;
- the DISCAMB form-factor path;
- `mode:accumulate`;
- other compilers (ifort/ifx, macOS);
- other BLAS libraries;
- realistic core counts and problem sizes.

## 3. Build three binaries

Prerequisites:
- **Compilers:** gfortran or ifort/ifx, plus a C/C++ compiler.
- **Build tools:** CMake ≥ 3.12.
- **Libraries:** FFTW3 (and its OpenMP library `fftw3_omp` if available), HDF5 with Fortran bindings, readline, X11/PNG/zlib (graphics build), **PGPLOT**, and a BLAS (OpenBLAS, MKL or Accelerate).
- **Python 3** with `numpy` and `h5py`.

Two upstream quirks:
- `FC` must be set in the environment.
- PGPLOT is required even with `-DDIFFUSE_GRAPHICS=OFF`, because `suite/prog` requires it unconditionally. On Ubuntu: `apt install pgplot5` and `export PGPLOT_DIR=/usr/lib/pgplot5`.

If the owner already builds DiffuseCode on this machine, reuse their usual
build settings (compiler, flags, library paths) for all three builds, and note them.

```sh
mkdir -p ~/discus_validation && cd ~/discus_validation
git clone https://github.com/jmsweng/DiffuseCode src
git -C src checkout perf/gemm-fourier
git -C src worktree add ../src-master master
export FC=$(which gfortran)            # or ifx/ifort
cmake -S src-master -B build-master                    && cmake --build build-master -j
cmake -S src        -B build-branch                    && cmake --build build-branch -j
cmake -S src        -B build-blas -DDISCUS_BLAS=ON     && cmake --build build-blas   -j
```

Check that the right compile definitions are active, for example:
```sh
grep -o "\-DDISCUS_BLAS\|-DDIFFUSE_FFTW_OMP" build-blas/discus/prog/CMakeFiles/discus_all.dir/flags.make | sort -u
```
`build-blas` should show both, `build-branch` only `-DDIFFUSE_FFTW_OMP` (or
neither if FFTW's OpenMP library is missing), `build-master` neither. Also
record which BLAS CMake found (`grep BLAS build-blas/CMakeCache.txt`).

Any compile failure in `build-branch` or `build-blas`, but not in `build-master`, is a finding.
Report it with the full error text.

## 4. Run the regression suite

```sh
R=src/discus/benchmarks/regression
M=build-master/suite/prog/discus_suite
python3 $R/regress.py --ref $M --test build-branch/suite/prog/discus_suite --work work_branch --threads <ncores> | tee regress_branch.txt
python3 $R/regress.py --ref $M --test build-blas/suite/prog/discus_suite   --work work_blas   --threads <ncores> | tee regress_blas.txt
```

Expected: `ALL CHECKS PASSED` for both. How to read failures:

- **`exact` failure:** the branch differs from master in code other than the GEMM kernel. This is the most important kind; report the case, the numbers, and the `.log` files in `work_*/ref` and `work_*/test_table`.
- **`oracle` failure with `build-blas`:** GEMM or NUFFT is inaccurate on this platform (BLAS or FINUFFT build). Rerun that case with `--only <case>` and `--threads 1` to see whether threading matters.
- **`DISCUS ERRORS` lines:** a macro command failed. Read the `.log`. If it also fails for `ref`, it is a syntax difference in this DISCUS version, not a regression.
- **Timing:** the suite's timings are dominated by start-up and are not a benchmark.

## 5. Run the owner's own workflows (most important)

Ask the owner which macros to use: typical production jobs, ideally including
3D diffuse scattering, 3D-ΔPDF output, any use of `set technique:nufft`, `lots`,
`set aver`, anisotropic ADPs, molecules, DISCAMB, and `mode:accumulate`.
Then, for each one:

1. Copy the macro and every file it reads into its own subdirectory.
2. Separate structure generation from the Fourier part if the macro builds or
   modifies the structure with random steps (`replace`, `mmc`, `ther`,
   `domain`, random `insert`, ...). Run the generation once, `save` the
   structure, and let both binaries `read`/`stru` the saved file. Otherwise
   the inputs differ and nothing can be compared.
3. Random lots (`lots box,...,<n>,...`) are not reproducible even in master
   (`seed` does not reach the lot placement). For a point-by-point comparison,
   temporarily use `lots ..., all, ...` or `lots off`. Otherwise compare only
   statistics: total intensity, and correlation (should be > 0.999).
4. Where the macro writes output, add (or switch to) `form common` so the full
   volume is written in double precision. Keep the original output too if
   the owner relies on it.
5. Run the copy three times, each into its own output directory:
   - `build-master`;
   - `build-branch`;
   - `build-blas` (and once more with `DISCUS_FOURIER_KERNEL=table` if anything looks off).

   Use the same `OMP_NUM_THREADS` for all three.
6. Compare: `python3 $R/compare_outputs.py out_master/ out_branch/` (and `out_blas/`).

Expected outcomes:

| Run | vs master |
|---|---|
| turbo, `build-branch` | ≤1e-12 of the maximum (identical) |
| turbo, `build-blas` with `DISCUS_FOURIER_KERNEL=table` | ≤1e-12 |
| turbo, `build-blas` | about 1e-6 to 1e-3 of the maximum, correlation > 0.9999. Larger with `set aver` because the maximum is small. This is the old table kernel's error |
| nufft, any branch build, step × cells integer on every axis | ≤1e-6 |
| nufft, any branch build, step × cells **not** an integer | Large differences expected (old code was wrong). Confirm the new result against `set technique:turbo` on the same grid: should agree to ≤1e-5 (≤1e-3 for `build-branch`) |
| 3D-PDF / 3D-ΔPDF | same tolerance as the intensities it came from |

Anything outside these ranges is a finding: report the macro (or a minimal
version), the numbers, and where in the volume the largest differences sit
(index and hkl). The HDF5 `common` file stores `data_corner` and
`data_increment_vector`; hkl = corner + index · increment.

## 6. Earlier NUFFT results

Ask the owner whether they have published or saved results made with
`set technique:nufft`. For each, compute step × unit cells along every
non-flat axis (step = increment vector length in r.l.u.; cells = crystal size
`cr_icc`). If any product is not a whole number, the old result was wrong:
rerun with `build-branch` and with `technique:turbo` and report both.

## 7. Performance on this machine

Use one or two of the owner's real jobs (large crystal, 3D grid). Time the
Fourier step (DISCUS prints `Elapsed time` after `run` in the fourier menu) and
the 3D-PDF output step (wall clock around `output ... run`). Do this for:

| Build | Settings |
|---|---|
| `build-master` | `OMP_NUM_THREADS=<ncores>` |
| `build-branch` | same |
| `build-blas` | `OMP_NUM_THREADS=<ncores>` and `OPENBLAS_NUM_THREADS=<ncores>` (or `MKL_NUM_THREADS`) |
| `build-master` and `build-branch` | `OMP_NUM_THREADS=1` (serial kernel) |
| `technique:nufft`, if the owner uses it | master vs branch |

Report a table with core count, CPU model, crystal size (atoms), grid size,
and times. For reference, cloud results on 2 cores:
- **GEMM, 16k atoms, 101³:** 7.2 s → 2.3 s.
- **Serial, same job:** 17.3 s → 6.0 s.
- **3D-PDF FFT, 301³ (standalone):** 2.6 s → 0.5 s.

If GEMM is slower than master on many cores, try
`OPENBLAS_NUM_THREADS` equal to physical (not logical) cores, and report both.

## 8. Report back

Write `~/discus_validation/REPORT.md` with:
1. **Machine:** OS, CPU and core count, compiler and version, BLAS and FFTW versions, CMake options used.
2. **Build results** for the three builds.
3. **Suite results:** the summary lines of `regress_branch.txt` and `regress_blas.txt` (full files attached), and every non-`ok` line.
4. **Owner workflows:** for each, what it does, the comparison table from `compare_outputs.py`, and your verdict against the expected outcomes in section 5.
5. **NUFFT re-runs** from section 6.
6. **Performance table** from section 7.
7. **Findings** (anything outside expectations), each with a minimal reproduction.

Don't change the branch to make a check pass. Report instead.
