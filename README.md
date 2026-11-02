> **This fork (`jmsweng/DiffuseCode`, branch `perf/refine`)** fixes six bugs in
> structure refinement (REFINE, RMC, KUPLOT, HDF5 output) and makes REFINE need
> fewer runs of the user macro, which dominate the cost of a real refinement.
> Each fix and each speed-up is in its own commit on top of upstream `master`, so
> they can be reviewed and merged one at a time.
> See **[refine/CHANGES.md](refine/CHANGES.md)** for the changes, their commits,
> validation and macro-run counts. The original project is
> [tproffen/DiffuseCode](https://github.com/tproffen/DiffuseCode).

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
