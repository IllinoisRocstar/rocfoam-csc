RocFOAM
-------
This project is a coupling module for OpenFOAM using IMPACT multiphysics coupling tool. 

## Version ##
Version 1.0.2

RocFOAM follows semantic versioning. The versions will be major.minor.patch.
We will:
* Increase the patch version for bug fixes, security fixes, and code
documentation. Backwards compatible; no breaking changes.
* Increase the minor version for new features and additions to the library’s
interface. Backwards compatible; no breaking changes.
* Increase the major version for breaking changes to the library’s interface or
breaking changes to behavior.

## Options ##
The following table contains all available CMake options to configure the project. The necessary third-party libraries are listed in the notes section.

| Option name            | Option description              | Default | Notes                            |
|------------------------|---------------------------------|---------|----------------------------------|
| ENABLE_TESTING         | Enable Testing system           | ON      |                                  |

## Linux Build Instructions ##

### Dependencies ###
You will need the following dependencies:
* MPI compiler
* OpenFOAM v2006 (available on [official website](https://www.openfoam.com/download/). Choose an option that includes the `dev` (Debian) or `devel` (RPM) package.)
* IMPACT (available on [GitHub](https://github.com/IllinoisRocstar/IMPACT))
* Boost v1.45.0 or higher (Must include development headers)

### Build ###

First, load the OpenFOAM environment. On a default installation, this is done by:
```
$ . /usr/lib/openfoam/openfoam2006/etc/bashrc
```
See OpenFOAM wiki for [alternatives](https://develop.openfoam.com/Development/openfoam/-/wikis/running).

We assume the following environment variables:
* `IMPACT_INSTALL_PATH` is the install location of IMPACT
* `ROCFOAM_INSTALL_PATH` is the desired installation location

Recommended configuration is as follows:
```
$ mkdir build && cd build
$ cmake .. \
        -DBUILD_SHARED_LIBS=ON \
        -DCMAKE_INSTALL_PREFIX=${ROCFOAM_INSTALL_PATH} \
        -DIMPACT_DIR=${IMPACT_INSTALL_PATH}/lib/cmake/IMPACT \
        -DCMAKE_CXX_COMPILER=mpicxx
```

To build, test, and install:
```
$ make -j
$ make test
$ make install -j
```
