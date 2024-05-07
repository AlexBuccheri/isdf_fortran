# Prototyping ISDF Implementation in Fortran

CMake configure, build and test commands:

## OMP

```shell
cmake --fresh -B serial-cmake-build-debug -DCMAKE_BUILD_TYPE=Debug -DMPI=Off
cmake --build serial-cmake-build-debug
ctest --test-dir ./serial-cmake-build-debug --output-on-failure
```

## MPI

MPI is OFF by default

```shell
cmake --fresh -B cmake-build-debug -DCMAKE_BUILD_TYPE=Debug  
cmake --build cmake-build-debug
ctest --test-dir ./cmake-build-debug --output-on-failure
```

## TODOs

* Make my k-means clustering repo a submodule of this
  
### Serial

* Construction of $P^phi$
* Construction of (CC^T)^-1
* Element-wise product $P^phi$ . $P^psi$ . (CC^T)^-1
* Face-splitting product
* Mocking of some one-particle states
  
### MPI

* Consideration of algorithms for state and domain distribution
* Use of mocked, batched data structures
* Trial of SCALAPACK
