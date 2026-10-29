# Fourier / 3D-PDF performance benchmarks

Standalone programs used to evaluate faster alternatives to the DISCUS
"turbo" Fourier kernel (`four_strucf_*` in `discus/prog/four_strucf_OMP_mod.f90`)
and the 3D-PDF FFT (`out_prep_3dpdf_3d` in `discus/prog/graphic.f90`).
None of these files are part of the DISCUS build.

| File | What it does |
|---|---|
| `bench_strucf.f90` | Replicas of `four_strucf_serial` and `four_strucf_omp`, plus a phase-recurrence kernel (pure Fortran, SIMD over atoms) and a separable-phase ZGEMM kernel; checks all against an exact sum |
| `test_gemm.f90` | Validates `discus/prog/four_strucf_gemm_mod.f90` against the repo kernel and an exact sum on a skewed grid with non-zero `lmn` / `off_shift` |
| `bench_pdf.c` | 3D FFT as done now (complex-to-complex, 1 thread, `FFTW_ESTIMATE`) vs real-to-complex and threaded FFTW |

## Using the GEMM kernel in DISCUS

Configure with BLAS enabled (OpenBLAS or MKL; found via CMake's FindBLAS):

```sh
cmake -S DiffuseCode -B build -DDISCUS_BLAS=ON
cmake --build build -j
```

`four_strucf` then uses `four_strucf_gemm` for every grid with more than one
point along both abscissa and ordinate (layers and volumes); lines keep the
original kernel. To compare against the original lookup-table kernel in the
same binary:

```sh
DISCUS_FOURIER_KERNEL=table discus_suite -macro my.mac
```

Threads: for volumes the layers are distributed over the OpenMP threads
(`OMP_NUM_THREADS`) and the BLAS library is switched to one thread meanwhile
(OpenBLAS and MKL on Linux, see `four_blas_threads.c`; elsewhere set
`OPENBLAS_NUM_THREADS=1` / `MKL_NUM_THREADS=1` / `VECLIB_MAXIMUM_THREADS=1`).
A single layer is one GEMM and uses the BLAS library's own threads.

**Timing note.** The `Elapsed time` DISCUS prints after `fourier run` is
`CPU_TIME`, i.e. CPU time summed over all threads, not wall-clock time. It does
not go down with more threads, and spinning BLAS/OpenMP threads inflate it.
Always measure wall-clock time (e.g. `/usr/bin/time`, or `date` around the run).

Wall-clock check inside DISCUS (16,384 atoms, Cu/Au random alloy, X-ray,
`lots off`, 101³ grid -5..5, median of 3 runs, 2-core VM, includes ~0.7 s
start-up and structure reading):

| Build | 1 thread | 2 threads |
|---|---|---|
| master | 24.7 s | 4.4 s |
| branch, default (table kernel) | 8.1 s | 4.4 s |
| branch, `-DDISCUS_BLAS=ON` | 2.7 s | 1.6 s |

## Standalone benchmarks: build and run

```sh
# needs gfortran, OpenBLAS, FFTW3 (with the OpenMP threads library)
gfortran -O3 -ffast-math -fopenmp bench_strucf.f90 -lopenblas -o bench_strucf
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 ./bench_strucf 25 101   # cells per edge, grid points per axis

gcc -O2 -c ../prog/four_blas_threads.c
gfortran -O3 -ffast-math -fopenmp ../prog/four_strucf_gemm_mod.f90 test_gemm.f90 four_blas_threads.o -lopenblas -o test_gemm
./test_gemm

gcc -O2 bench_pdf.c -lfftw3_omp -lfftw3 -fopenmp -lm -o bench_pdf
./bench_pdf 301
```

## Results (2-core Xeon, AVX-512, gfortran -O3 -ffast-math, 2 threads)

| Kernel | 8k atoms, 61³ | 15.6k atoms, 101³ |
|---|---|---|
| `four_strucf_serial` (repo) | 3.7 s | 48.2 s |
| `four_strucf_omp` (repo) | 1.5 s | 12.9 s |
| Phase recurrence, SIMD over atoms | 0.76 s | 3.7 s |
| Separable phases + ZGEMM | 0.35 s | 2.3 s |

`test_gemm` (20k atoms, 81×71×61 skewed grid): 4.72 s → 0.92 s; max error vs
exact sum 1.8e-2 (2^16 lookup table) → 1.1e-11.

3D-PDF FFT, 301³: 2.56 s / 0.87 GB (current) → 0.52 s / 0.22 GB (r2c, threaded).
