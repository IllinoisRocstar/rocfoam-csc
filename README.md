RocFOAM
-------
This project is a coupling module for OpenFOAM using IMPACT multiphysics coupling tool. 

## Version ##
Version 1.0.0

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
* OpenFOAM v2006 (available on [official website](https://www.openfoam.com/download/))
* IMPACT (available on [GitHub](https://github.com/IllinoisRocstar/IMPACT))
* Boost 

### Build ###

First, load the OpenFOAM environment. See OpenFOAM wiki for [detailed instructions](https://develop.openfoam.com/Development/openfoam/-/wikis/running).

We assume the following environment variables:
* `IMPACT_INSTALL_PATH` is the install location of IMPACT
* `ROCFOAM_INSTALL_PATH` is the desired installation location
Recommended configuration is as follows:
```
$ cmake .. \
        -DBUILD_SHARED_LIBS=ON \
        -DCMAKE_INSTALL_PREFIX=${ROCFOAM_INSTALL_PATH} \
        -DCMAKE_IMPACT_DIR=${IMPACT_INSTALL_PATH}/lib/cmake/IMPACT \
        -DCMAKE_CXX_COMPILER=mpicxx
```

To build, test, and install:
```
$ make -j
$ make test
$ make install -j
```

## Using RocFOAM ##

### Running Pre-processor ###
To run, there are two standalone drivers called `rocRhoCentral` (density-based solver) and `rocRhoPimple` (pressure-based solver). They can be used in serial and parallel modes to preprocess a simulation case.
* Serial:
  ```
  rocRhoCentral -preprocess
  ```
  OR
  ```
  rocRhoPimple -preprocess
  ```
* Parallel:
  ```
  mpirun -np <nproc> rocRhoCentral -parallel -prerocess
  ```
  OR
  ```
  mpirun -np <nproc> rocRhoPimple -parallel -prerocess
  ```
where `<nproc>` needs to be replaced with number of processors.

### Validating Pre-processor Output ###
Optionally, to test if the pre-processor has been successfully executed, one could run the commands mentioned above, but replace `-preprocess` with `-production`, that is:
* Serial:
  ```
  rocRhoCentral -production
  ```
  OR
  ```
  rocRhoPimple -production
  ```
* Parallel:
  ```
  mpirun -np <nproc> rocRhoCentral -parallel -production
  ```
  OR
  ```
  mpirun -np <nproc> rocRhoPimple -parallel -production
  ```
where `<nproc>` needs to be replaced with number of processors.
