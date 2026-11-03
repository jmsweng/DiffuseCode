# Regression tests for the Fourier / 3D-PDF changes

Two scripts, both need Python 3 with `numpy` and `h5py`.

## `regress.py`: built-in test cases

Runs 43 DISCUS cases (91 checks) with a binary built from `master` (reference) and a
binary built from this branch (test). It writes every result as an HDF5 `common`
file (full volume, double precision) and compares all points.

```sh
python3 regress.py --ref  /path/to/master-build/suite/prog/discus_suite \
                   --test /path/to/branch-build/suite/prog/discus_suite \
                   --work /tmp/discus_regress --threads 8
# --only 'n_*'   run a subset (glob on case names)
# --keep         reuse outputs already in --work
```

Exit status 0 means every check passed. It takes about 1 minute on 2 cores.

| Check | Compares | Tolerance (of max) | What it shows |
|---|---|---|---|
| `exact` | master vs branch with `DISCUS_FOURIER_KERNEL=table` | 1e-12 | Every change other than the GEMM kernel reproduces master: serial kernel, Debye-Waller cache, 3D-PDF FFT, NUFFT with `set aver` on valid grids |
| `kernel` | master vs branch default kernel | 1e-2 | Gross-error guard only. See below |
| `oracle` | branch NUFFT vs branch turbo | 1e-5 (1e-3 without BLAS) | Two independent exact-phase methods agree |
| `mask` | branch NUFFT with `set aver` vs without | 1e-12 | `aver` only zeroes the integer-hkl points |
| `info` | reported only | none | e.g. how wrong the old NUFFT was |

The cases cover:
- **Grids:** odd and even sizes (with and without the Friedel shortcut), unequal axis lengths, off-centre and skewed grids, layers, lines.
- **Radiation:** X-ray, neutron, anomalous dispersion.
- **ADPs:** isotropic and anisotropic, with one or two ADP types.
- **Averaging:** `set aver`, and `lots` (box/ellipsoid, periodic or not).
- **Crystals:** a thermally displaced crystal and a one-cell-thick crystal.
- **Debye-Waller cache:** one session with seven setup changes, compared against the same session in master and against fresh sessions.

### The `kernel` check and the original table kernel

The original turbo kernel looks up sin/cos in a 2^16-entry table. Compared with
GEMM and NUFFT (which agree to 1e-8 .. 3e-6), its error reaches:
- 1e-3 of the maximum in layers;
- 2e-3 with `set aver`;
- up to ~1% at single weak diffuse points (3D-ΔPDF: 3e-4 of the maximum).

The `kernel` check only guards against gross errors. Correctness of GEMM is shown
by the `oracle` checks and by `four_strucf_gemm` vs an exact sum (`../test_gemm.f90`).

### Known non-reproducibility in master

Random lots (`lots box, ..., <n>, ...`) pick different lot origins on every run,
even after `seed`. Random lots therefore can't be compared point by point; the
cases use `lots ..., all, ...`, which is deterministic.

## `compare_outputs.py`: compare your own runs

```sh
python3 compare_outputs.py master_out/ branch_out/       # files with equal names
python3 compare_outputs.py master.h5 branch.h5 --tol 1e-5
```

It reads HDF5 (`common`, `hdf5`, `nexus`) and text formats (`gnu`, `xplor`, ...).
Text formats round to ~5 significant digits, so use `form common` where possible.

See `../LOCAL_VALIDATION.md` for the full procedure on a new machine.
