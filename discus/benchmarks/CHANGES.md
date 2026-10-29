# DiffuseCode: faster 3D diffuse scattering and 3D-ΔPDF

**Branch:** [`perf/gemm-fourier`](https://github.com/jmsweng/DiffuseCode/tree/perf/gemm-fourier) on `jmsweng/DiffuseCode`, a fork of `tproffen/DiffuseCode`
**Base:** `master` @ `19a214a` (DiffuseCode 6.23.01)
**Status (29 Oct 2026):** validated on a cloud VM and on a Ryzen 7 4800HS workstation. All 82 regression checks pass for both build variants.

## Summary

Computing a 3D diffuse scattering volume in DISCUS is dominated by one routine,
`four_strucf`. For every grid point it sums the phase factors of all atoms. The
3D-ΔPDF is then a single FFT of that volume, which takes seconds.

This branch:

1. adds a matrix-product (ZGEMM) version of that sum. It is opt-in, uses exact phases and is up to ~15× faster;
2. makes the single-threaded fallback ~5× faster;
3. fixes a **correctness bug** in the NUFFT technique, which silently gave wrong intensities on many grids;
4. fixes the Debye-Waller table cache;
5. makes the 3D-PDF FFT real-to-complex and multithreaded;
6. adds a full-volume regression suite and validation tooling.

Without `-DDISCUS_BLAS=ON`, results from the turbo technique are identical to `master`.
Results from the NUFFT technique change only where the old code was wrong.

## Changes

| # | Change | Files | When active | Effect |
|---|---|---|---|---|
| 1 | **ZGEMM Fourier kernel.** On a regular grid the phase factor separates, exp(2πi[p₀+(j−1)u·x+(i−1)v·x+(h−1)w·x]) = U·V·W<sub>h</sub>, so each layer of the structure factor is one complex matrix product with the atoms as the summed index. Phases come from exact sin/cos instead of the 2¹⁶-entry lookup table. Layers are spread over OpenMP threads with single-threaded BLAS | `four_strucf_gemm_mod.f90`, `four_blas_threads.c`, `four_strucf_OMP_mod.f90`, `CMakeLists.txt` | Only with `cmake -DDISCUS_BLAS=ON`. Layers and volumes; lines keep the original kernel. `DISCUS_FOURIER_KERNEL=table` switches back at run time | ~15× faster on 1 thread, ~2.2× at 8 threads; error vs exact sum 1.8e-2 → 1.2e-11 |
| 2 | **Serial kernel.** Accumulates into a per-plane buffer instead of writing `tcsf(j,i,h)` with a stride of `num(1)·num(2)` in the innermost loop. The inner loop now runs to `fnum(3)` instead of `num(3)` | `four_strucf_OMP_mod.f90` | Always, when OpenMP is off or limited to 1 thread | ~5× faster, bit-identical |
| 3 | **NUFFT driver rewritten** (bug fix). The old code scaled atom coordinates by an integer rounded up from step × unit cells, so whenever that product was not a whole number the grid had the wrong step. Example: step 0.2 with 16 cells (3.2, rounded to 4) put the strongest peak at (−0.8, −0.8, −4) instead of (0, 0, −4). Its error check could never trigger. New `four_strucf_nufft` uses coordinates 2π·vᵢ·x wrapped into [−π, π] and puts the corner phase into the source strengths, which is exact for any step, crystal size, grid position and skewed axes. The old oversample-then-discard step is gone | `fourier_sup.f90`, `fourier_finufft.f90`, `appl_dis.hlp` | Always, for `set technique:nufft` (not the DISCAMB variant) | Correct on every grid; agrees with turbo to ≤3e-6 of the maximum |
| 4 | **Debye-Waller table cache.** Fixed a size check comparing `ubound(four_dbw,2)` with `num(1)`, and a test that compared the ADPs with a tolerance instead of with the previous ADPs, so the table was rebuilt on every run. The cache key now includes the grid, lattice and `temp` setting; `four_dbw_valid` is cleared whenever the table is reallocated | `fourier_sup.f90`, `diffuse_mod.f90`, `discus_allocate_appl_mod.f90` | Always, with `temp use` | Unchanged setups reuse the table; changed ones recompute (bit-identical to a fresh session) |
| 5 | **3D-PDF FFT.** Real-to-complex transform instead of copying the intensity into a complex array; the full real part is rebuilt from the half spectrum. Multithreaded when FFTW's OpenMP library is found | `graphic.f90`, `CMakeLists.txt` | Always, for 3D `value 3DPDF` | Output identical; 4× less memory; ~5× faster at 301³ |
| 6 | **Validation tooling.** `regress.py`: 43 cases with full 3D volumes in double precision, master vs branch. `compare_outputs.py` compares any two DISCUS outputs. Also standalone kernel benchmarks, a GEMM-vs-exact-sum test, a README, and a step-by-step guide for validating on another machine | `discus/benchmarks/` | Not part of the DISCUS build | — |

Also corrected: the help text for `set technique:nufft`, which now describes when the old restrictions still apply (only with `set aver`).

### Commits

| Commit | Description |
|---|---|
| `80922b3` | Add ZGEMM-based Fourier kernel and performance benchmarks |
| `3d9598f` | Wire ZGEMM Fourier kernel into DISCUS behind `DISCUS_BLAS` option |
| `87c8346` | `four_dbwtab`: fix cache check so the Debye-Waller table is reused correctly |
| `80d0980` | `four_strucf_serial`: accumulate into a per-plane scratch array |
| `ba1cf16` | 3D-PDF: real-to-complex FFT, multithreaded FFTW |
| `20c16cb` | NUFFT: exact coordinate scaling for any grid, no oversampling |
| `edc8011` | Add full-volume regression suite |
| `bd4b738` | Add local validation instructions |
| `b1761c0` | `regress.py`: `n_line` oracle uses table-kernel tolerance |
| `507b2e6` | GEMM kernel: distribute layers over OpenMP threads, single-threaded BLAS |
| `41bd9d0` | `regress.py`: table-kernel tolerance 5e-3; always run the GEMM detection case |

19 files changed, +1784 / −250 lines. About half of the additions are tests and benchmarks.

## Timing benchmarks

> **How DISCUS reports time.** The `Elapsed time` DISCUS prints after `fourier run`
> is `CPU_TIME`: CPU time summed over all threads, not wall-clock time. It never
> drops with more threads and grows when idle threads spin. All tables below are
> wall-clock time unless marked otherwise.

### 1. Full DISCUS run: owner's workstation

AMD Ryzen 7 4800HS (8 cores / 16 threads), WSL Ubuntu 22.04, gfortran 11.4,
OpenBLAS 0.3.28. Job: 16,384-atom Cu/Au random alloy (16³ fcc cells), X-ray,
101³ grid from −5 to 5, `lots off`, turbo technique. Median of 3 runs, measured
with `/usr/bin/time`, including ~0.4 s start-up and structure reading.

| Build | 1 thread | 2 threads | 4 threads | 8 threads | 16 threads |
|---|---|---|---|---|---|
| `master` | 34.97 s | 4.21 s | 2.98 s | 2.21 s | 2.44 s |
| branch, default build | 7.45 s | 4.45 s | 2.80 s | 2.08 s | 2.12 s |
| branch, `-DDISCUS_BLAS=ON` | **2.65 s** | **1.80 s** | **1.36 s** | **1.23 s** | **1.26 s** |

Speedup over `master` at the same thread count, start-up (0.4 s) subtracted:

| Build | 1 thread | 2 threads | 4 threads | 8 threads | 16 threads |
|---|---|---|---|---|---|
| branch, default build | 4.9× | 0.9× | 1.1× | 1.1× | 1.2× |
| branch, `-DDISCUS_BLAS=ON` | **15.4×** | **2.7×** | **2.7×** | **2.2×** | **2.4×** |

- **Default build with more than one thread:** it runs the same parallel kernel as `master`, so it matches `master`. Its gain is the single-threaded case.
- **BLAS build with OpenBLAS forced to one thread** (16 OpenMP threads, `OPENBLAS_NUM_THREADS=1`): 1.23 s, the same as above. The kernel no longer depends on the BLAS library's own threads.
- **Where the BLAS build levels off:** at about 0.8 s past start-up for this job size. The remaining time is the single-threaded work DISCUS does around the kernel; larger crystals or grids benefit more.

### 2. Full DISCUS run: cloud VM

2-core Xeon (AVX-512), Ubuntu 24.04, gfortran 13.3, OpenBLAS. Same job as
above, median of 3 runs, including ~0.7 s start-up.

| Build | 1 thread | 2 threads |
|---|---|---|
| `master` | 24.7 s | 4.4 s |
| branch, default build | 8.1 s | 4.4 s |
| branch, `-DDISCUS_BLAS=ON` | **2.7 s** | **1.6 s** |

### 3. Fourier kernel alone (standalone replicas)

2-core cloud VM, 2 threads, random atom positions, gfortran `-O3 -ffast-math`.

| Kernel | 8,000 atoms, 61³ | 15,625 atoms, 101³ |
|---|---|---|
| `four_strucf_serial` (original) | 3.7 s | 48.2 s |
| `four_strucf_omp` (original) | 1.5 s | 12.9 s |
| Phase recurrence, SIMD over atoms (pure Fortran, not adopted) | 0.76 s | 3.7 s |
| Separable phases + ZGEMM (adopted) | **0.35 s** | **2.3 s** |

Accuracy against an exact sum (20,000 atoms, skewed 81×71×61 grid, max |F| = 324):
the lookup-table kernel's maximum error is 1.8e-2, the ZGEMM kernel's 1.2e-11.

### 4. 3D-PDF FFT (standalone, FFTW 3.3.10, 2-core VM)

| Grid | Original: complex-to-complex, 1 thread | Real-to-complex, 1 thread | Real-to-complex, 2 threads | Memory, original → new |
|---|---|---|---|---|
| 201³ | 0.77 s | 0.43 s | **0.21 s** | 0.26 → 0.07 GB |
| 301³ | 2.56 s | 1.14 s | **0.52 s** | 0.87 → 0.22 GB |

### 5. NUFFT (DISCUS, 2-core VM, CPU time from `Elapsed time`)

| Grid | Old NUFFT | New NUFFT | Old result |
|---|---|---|---|
| 41³, step 0.2, 16 cells | 0.54 s | 0.13 s | wrong (step × cells = 3.2) |
| 33³, step 0.25, 16 cells | 0.40 s | 0.11 s | correct (identical to new) |

On the owner's workstation, on a grid where the old code did not oversample, NUFFT time was unchanged (~1.0 s).

## Accuracy of the original kernel

Every grid tested gave the same picture. The ZGEMM kernel and the independent
NUFFT agree with each other to 1e-8 – 3e-6 of the maximum intensity. The
original lookup-table kernel deviates from both by:

| Case | Max deviation (fraction of max intensity) |
|---|---|
| 3D volumes | 1e-6 – 4e-5 |
| 2D layers | up to 1.5e-3 |
| with `set aver` (diffuse only) | up to 2.4e-3 |
| 3D-ΔPDF with `set aver` | ~3e-4 |
| single weak diffuse points | up to ~1% |

So switching to `-DDISCUS_BLAS=ON` changes results by about this much, and the
new values are the more accurate ones.

## Validation

- **Regression suite** (`regress.py`, 43 cases, 82 checks plus 9 informational):
  - **Coverage:**
    - grids: odd/even/unequal/off-centre/skewed 3D grids, layers, lines;
    - radiation: X-ray, neutron, anomalous dispersion;
    - ADPs: isotropic and anisotropic, with one or two ADP types;
    - averaging: `set aver`, and `lots` (box/ellipsoid, periodic or not);
    - crystals: thermally displaced, and one cell thick;
    - outputs: 3D-PDF and 3D-ΔPDF;
    - the serial kernel, and a 7-step Debye-Waller cache session.
  - **Results:** all checks pass for both build variants on both machines. With the original kernel selected, the branch reproduces `master` to ≤1e-12 of the maximum (mostly bit for bit).
- **Not tested** (no data available): molecules, electron radiation, DISCAMB form factors, `mode:accumulate`, compilers other than gfortran.
- **Upstream quirk found:** random `lots` are not reproducible in `master`, even after `seed`.

## Building

```sh
export FC=$(which gfortran)
cmake -S DiffuseCode -B build -DDISCUS_BLAS=ON    # omit -DDISCUS_BLAS=ON for the default build
cmake --build build -j
```

- **Libraries:** BLAS (OpenBLAS, MKL or Accelerate) for the GEMM kernel; FFTW with `fftw3_omp` for the threaded 3D-PDF.
- **PGPLOT** is required even with `-DDIFFUSE_GRAPHICS=OFF` (an upstream quirk).

Run the regression suite:

```sh
python3 discus/benchmarks/regression/regress.py \
    --ref  /path/to/master-build/suite/prog/discus_suite \
    --test /path/to/branch-build/suite/prog/discus_suite --threads 8
```

## Still open

- `four_run_nufft_discamb` (DISCAMB form factors) has the same NUFFT scaling bug; untested for lack of DISCAMB data.
- `-march=native` is set for C++ but not for Fortran.
- The Friedel mirror is reapplied for every lot.
- The CUDA kernel launches once per atom.
- The NUFFT tolerance is hard-coded at 1e-7.
- `four_dbw` uses num³ × nanis × 8 bytes of memory.
- `do_fft_2d_cos` is dead code with index errors.
