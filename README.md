# Prototyping ISDF Implementation in Fortran

CMake configure, build and test commands:

## OMP

```shell
cmake --fresh -B serial-cmake-build-debug -DCMAKE_BUILD_TYPE=Debug
cmake --build serial-cmake-build-debug
ctest --test-dir ./serial-cmake-build-debug --output-on-failure
```

## MPI

MPI is OFF by default

```shell
cmake --fresh -B cmake-build-debug -DCMAKE_BUILD_TYPE=Debug -DMPI=On
cmake --build cmake-build-debug
ctest --test-dir ./cmake-build-debug --output-on-failure
```

## TODOs

* Make my k-means clustering repo a submodule of this
* Create some mock states - output from my python module

### Serial

* Write construction of `P` matrix
* Write construction of `(CC^T)^-1` matrix
* Write `construct_interpolation_vectors`, using `P` and `(CC^T)^-1` matrices

### MPI

* Consideration of algorithms for state and domain distribution
* Use of mocked, batched data structures
* Trial of SCALAPACK
