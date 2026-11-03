# Changes on this fork (jmsweng/DiffuseCode)

This fork of [tproffen/DiffuseCode](https://github.com/tproffen/DiffuseCode) adds two
sets of changes to DiffuseCode 6.23.01. Each set was developed and validated on its own
branch and then merged into `master`. Both branches are kept unchanged for reference:

| Branch | Area | Detailed description | Compare with upstream |
|---|---|---|---|
| [`perf/gemm-fourier`](https://github.com/jmsweng/DiffuseCode/tree/perf/gemm-fourier) | DISCUS: 3D diffuse scattering, 3D-PDF / 3D-ΔPDF | [discus/benchmarks/CHANGES.md](discus/benchmarks/CHANGES.md) | [19a214a…perf/gemm-fourier](https://github.com/jmsweng/DiffuseCode/compare/19a214a...perf/gemm-fourier) |
| [`perf/refine`](https://github.com/jmsweng/DiffuseCode/tree/perf/refine) | REFINE, RMC, KUPLOT, HDF5: structure refinement | [refine/CHANGES.md](refine/CHANGES.md) | [19a214a…perf/refine](https://github.com/jmsweng/DiffuseCode/compare/19a214a...perf/refine) |

Upstream base: `19a214a` (DiffuseCode 6.23.01). `perf/refine` also contains the
upstream commit `6fd5f33` (refine_data_mod.f90 moved to lib_f90), which came in with the merge.

Every fix and speed-up is a separate commit on its branch, so each one can be read,
tested or reverted on its own. The two branches touch different source files (only
`README.md` overlapped), so the merge changed no code.

## 1. 3D diffuse scattering and 3D-ΔPDF (`perf/gemm-fourier`)

| # | Change | Kind | When active | Commit(s) |
|---|---|---|---|---|
| 1 | ZGEMM Fourier kernel: each layer of the structure factor as one complex matrix product, exact phases instead of the lookup table, layers over OpenMP threads | speed-up | only with `cmake -DDISCUS_BLAS=ON`; `DISCUS_FOURIER_KERNEL=table` switches back at run time | `80922b3`, `3d9598f`, `507b2e6` |
| 2 | Serial Fourier kernel accumulates into a per-plane buffer (~5× faster on one thread, bit-identical) | speed-up | OpenMP off or 1 thread | `80d0980` |
| 3 | NUFFT driver rewritten: the old coordinate scaling gave wrong intensities whenever step × cells was not a whole number | **bug fix** | `set technique:nufft` (not DISCAMB) | `20c16cb` |
| 4 | Debye-Waller table cache: wrong size check, table rebuilt every run, incomplete cache key | bug fix | `temp use` | `87c8346` |
| 5 | 3D-PDF FFT real-to-complex and multithreaded (output identical, 4× less memory) | speed-up | `value 3DPDF` | `ba1cf16` |
| 6 | Full-volume regression suite, benchmarks and validation guide (`discus/benchmarks/`) | tooling | not part of the build | `80922b3`, `edc8011`, `bd4b738`, `b1761c0`, `41bd9d0` |

Measured on an 8-core Ryzen, 16 384 atoms, 101³ grid (wall clock): 35.0 s → 2.65 s on
one thread and 2.21 s → 1.23 s at 8 threads with the BLAS build. Without
`-DDISCUS_BLAS=ON`, turbo results are identical to upstream.

## 2. Structure refinement (`perf/refine`)

| # | Change | Kind | When active | Commit |
|---|---|---|---|---|
| 1 | MPI 5-point derivatives used the last parameter's step for all parameters | **bug fix** | `points:5` under MPI | `2e95b93` |
| 2 | Seed restore did not reset the Gaussian generator's cache, so derivatives were noisy | **bug fix** | models with Gaussian random numbers | `df2b814` |
| 3 | RMC with lots: no moves possible, wrong lot updated, lot 1 counted twice, wrong average bounds | **bug fix** | RMC with `set lots` | `a299cfa` |
| 4 | `MAX(dP/sig)` divided by the curvature instead of the parameter sigma; convergence criterion 1 could fire early | bug fix | always | `368094c` |
| 5 | Skip the repeated centre run after the derivatives (bit-identical) | speed-up | except with a plot macro or MPI | `5d65bcb` |
| 6 | Compute derivatives only for accepted steps (bit-identical) | speed-up | always | `863bb4f` |
| 7 | New opt-in one-sided derivatives, `points:2` | feature | `newparam ..., points:2` | `f40b2b5` |
| 8 | KUPLOT returned a single point to the MPI master, so MPI refinements never finished | **bug fix** | MPI refinement of KUPLOT data | `d1edc81` |
| 9 | HDF5 write tested an uninitialised error flag (false "Error opening file") | bug fix | any HDF5 write after the first | `0d73677` |

A 4-parameter test fit needs 46–97 macro runs instead of 52–162; with 25 fixed cycles
131 instead of 262. An RMC change that removes the amplitude copies per move was
tried (`80ea6bb`) and reverted (`41b21e2`); see refine/CHANGES.md for why.

## Behaviour changes to be aware of

- NUFFT (`perf/gemm-fourier` #3): results change where the old code was wrong.
- Seed restore (`perf/refine` #2): any refinement or MC run that draws Gaussian random
  numbers after a seed restore gives different (correct) numbers.
- `MAX(dP/sig)` (`perf/refine` #4): the printed column changes; a fit that stopped on
  criterion 1 may run more cycles.

Validation details, timing tables and what was not tested are in the two detailed
documents linked above.
