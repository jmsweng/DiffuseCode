> **This fork (`jmsweng/DiffuseCode`, branch `perf/refine`)** fixes six bugs in
> structure refinement (REFINE, RMC, KUPLOT, HDF5 output) and makes REFINE need
> fewer runs of the user macro, which is what dominates a real refinement. The
> branch starts at upstream `master` (`6fd5f33`) and holds one self-contained
> commit per item, fixes first, so it can be merged whole or item by item. The
> original project is [tproffen/DiffuseCode](https://github.com/tproffen/DiffuseCode).
>
> | Commit | Type | Change |
> |---|---|---|
> | `2e95b93` | fix | MPI: 5-point derivatives used another parameter's step size |
> | `df2b814` | fix | Seed restore did not reset the cached Gaussian number, so models with Gaussian randomness got noisy derivatives |
> | `a299cfa` | fix | RMC with `set lots`: no move possible, wrong lot updated, lot 1 counted twice |
> | `368094c` | fix | Printed and used `MAX(dP/sig)` was computed from the curvature matrix, not the sigmas |
> | `5d65bcb` | speed-up | Skip the second macro run at the centre point each cycle (1 run per cycle) |
> | `863bb4f` | speed-up | Derivatives only for accepted steps (a rejected step costs 1 run instead of 2N+1) |
> | `f40b2b5` | speed-up | New opt-in `points:2` one-sided derivatives (N instead of 2N runs) |
> | `d1edc81` | fix | KUPLOT returned only one point of its curve to the MPI master (`All derivatives are zero`) |
> | `0d73677` | fix | HDF5 write: uninitialised `hdferr` made a worker's second write fail |
>
> Tested with fixed-seed 4-parameter Gaussian fits (serial and `mpirun -n 5`) and
> an RMC test on Cu 6x6x6, on gfortran 13.3 and 11.4. Not tested: MPI with DISCUS
> (3D) models, multi-data-set REFINE, RMC with ellipsoidal lots. The test kit and
> handoff documents are kept at the tag `perf-refine-handoffs-final`.

Diffuse program collection
================================================================

Welcome to the 'diffuse program' collection.
The website can be found at http://tproffen.github.io/DiffuseCode

Depending on the package you have downloaded you will find one or all of the
following program directories after unpacking the archive.

DISCUS_SUITE :  The Main program for DISCUS/KUPLOT/DIFFEV/REFINE
                This is the standard program that should be used
                and includes all other sections under a common roof top
AUTHOR  : R.B. Neder  (reinhard.neder@fau.de)

DISCUS  : Diffuse Scattering & Defect Structure Simulation
AUTHORS : R.B. Neder  (reinhard.neder@fau.de)
          Th. Proffen (tproffen@ornl.gov)

DIFFEV  : Generic refinement program using evolutionary algoritm
AUTHOR  : R.B. Neder  (reinhard.neder@fau.de)

KUPLOT  : General plotting program (well suited for DISCUS output)
AUTHOR  : Th. Proffen (tproffen@ornl.gov)
          R.B. Neder  (reinhard.neder@fau.de)

MIXSCAT : Program to generate differential PDFs from n/X data
AUTHOR  : C. Wurden K. Page A. Llobet
          Th. Proffen (tproffen@ornl.gov)

INSTALLATION
================================================================

Binary
------

Binary distributions as well as the source code can be found 
at https://github.com/tproffen/DiffuseCode/releases

Linux
-----

Instructions to build from source please check file [AAA_INSTALL_DISCUS_LINUX.pdf](Manual/AAA_INSTALL_DISCUS_LINUX.pdf).

OSX - Using homebrew
---

Instructions to build from source please check file [AAA_INSTALL_DISCUS_MACOS.pdf](Manual/AAA_INSTALL_DISCUS_MACOS.pdf).

Windows
-------

Recommend you install the binary found [here](https://github.com/tproffen/DiffuseCode/releases).

CYGWIN
------

Instructions to build from source please check file [AAA_INSTALL_DISCUS_CYGWIN.pdf](Manual/AAA_INSTALL_DISCUS_CYGWIN.pdf).


[![Build Status](https://travis-ci.org/tproffen/DiffuseCode.svg?branch=master)](https://travis-ci.org/tproffen/DiffuseCode)
