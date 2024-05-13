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
  * Read importing/exporting [guide](https://cmake.org/cmake/help/latest/guide/importing-exporting/index.html)
* Create some mock states - output from my python module

### Serial

* Write unit and/or regression tests

### MPI

* Write some distribution routines for states and domain.
  * Implement the strategy below
  * Note that [Hu et. al.](J.Chem.TheoryComput.2017, 13, 5420-5431) uses 2d block-cyclic distribution for everything from KS states to $(ZC^T)$, but uses 1D column cyclic for the final interpolation vectors
* Use of mocked, batched data structures
* Trial of SCALAPACK

### Parallelisation Strategy 

This assumes domain distribution and state distribution, so parallelisation over indices $i$ and $\mathbf{r}_k$.
Batching of states, or use of scalapack would change the suggested strategy.

#### `construct_quasi_density_matrix_R_intr`

* $P^{\varphi}(\mathbf{r}, \mathbf{r}_\mu) = \sum_{i=1}^{m}  \varphi_i(\mathbf{r}) \varphi_i(\mathbf{r}_\mu)$

When state-parallel (over $i$), this means $P^{\varphi}(\mathbf{r}, \mathbf{r}_\mu)$ has partial contributions
to all elements, on each process. After construction, perform `allreduce` on $P$ using the states communicator, 
such that it is fully-defined on all processes.

If domain-parallel (over $\mathbf{r}$, which we can index with $k$), $P$ stays distributed w.r.t. row index.
However, one still needs to define $\varphi_i(\mathbf{r}_\mu)$ for all interpolation grid points, on each process.
This changes the implementation from the serial API:

* Construct $\varphi_i(\mathbf{r}_\mu)$ and hold in memory on all processes, for all interpolation points $(\mu)$.
* As such, one does not pass `indices` of interpolation grid points to various routines.

#### `select_interpolation_points_of_first_dim`

Given the implementation points in the prior section, it's easiest to make a new routine that takes 
$P^{\varphi}(\mathbf{r}_\nu, \mathbf{r}_\mu)$ and $P^{\psi}(\mathbf{r}_\nu, \mathbf{r}_\mu)$, and computes the element-wise
product. Define on all processes. This is basically the same operations in `construct_quasi_density_matrix_R_intr`.

#### `construct_inverse_coefficients_contraction`

This stays the same, $P^(\mathbf{r}_\nu, \mathbf{r}_\mu)$ are defined on all processes, SVD is performed and the resulting
inverted matrix, $(CC^T)^{-1}$ is defined on all processes.

#### Contraction of $ZC^T$ with $(CC^T)^{-1}$

\[
\begin{bmatrix}
a_{11} & a_{12} \\
a_{21} & a_{22} \\
a_{31} & a_{32} \\
a_{41} & a_{42} \\
\end{bmatrix}
\times
\begin{bmatrix}
b_{11} & b_{12} \\
b_{21} & b_{22} \\
\end{bmatrix} 
\]

If rows 1-4 are on $p_0$ and rows 5-8 are on $p_1$, then the resulting product on $p_0$ is:

\[
\begin{bmatrix}
a_{11}b_{11} + a_{12}b_{21} & a_{11}b_{12} + a_{12}b_{22} \\
a_{21}b_{11} + a_{22}b_{21} & a_{21}b_{12} + a_{22}b_{22} \\
a_{31}b_{11} + a_{32}b_{21} & a_{31}b_{12} + a_{32}b_{22} \\
a_{41}b_{11} + a_{42}b_{21} & a_{41}b_{12} + a_{42}b_{22} \\
\end{bmatrix}
\]

To map the local row index to the global index, one could use an expression like:

`irow_global = irow + (iproc * np)`

where `np` is the number of grid points (hence rows) local to process `iproc`. 

The problem with row distribution is that one cannot use lapack and pass a sub-section of the output array to fill into, because
the memory will be noncontiguous, resulting in a copy.

Would be better if the full grid ran over the column index of Z... but this would also require changing the current implementation.

### Resulting Code Changes

Assuming both state and domain distributions:

* $\varphi_i(\mathbf{r})$ is distributed w.r.t. $i$ and $\mathbf{r}$. 
* $\varphi_i(\mathbf{r}_\mu)$ should be constructed such that each process has all $\{\mathbf{r}_\mu\}$, and remain distributed in $i$
* This allows $P^\varphi(\mathbf{r}_\nu, \mathbf{r}_\mu)$ to be constructed on all processes with no comms over domain.
  * Use allreduce on $P$ to sum the contributions from states across all processes.
* $P^\varphi(\mathbf{r}, \mathbf{r}_\mu)$ is built from $\varphi_i(\mathbf{r})$ distributed w.r.t. $(i, \mathbf{r})$ and $\varphi_i(\mathbf{r}_\mu)$ distributed w.r.t. $i$.
  * Use allreduce on $P$ to get sum the contributions from states across all processes.
  * $P^\varphi(\mathbf{r}, \mathbf{r}_\mu)$ stays distributed w.r.t. row index (i.e. over domains)
* $(CC^T)^{-1}$ constructed from $P(\mathbf{r}_\nu, \mathbf{r}_\mu)$, making it available on all processes
* $(ZC^T)$ is constructed from $P(\mathbf{r}, \mathbf{r}_\mu)$, and retains the same distribution as $P^\varphi(\mathbf{r}, \mathbf{r}_\mu)$.
  * Distributed over rows.
* $\Theta$ is constructed from the contraction of $(ZC^T)$ and $(CC^T)^{-1}$
  * $\Theta$ is required on all processes, so one should initialise $\Theta=0$ and compute the product for each block of rows, per process. Allreduce $\Theta$ will collate contributions from all processes.
  * Open question as to how to map local index of $(ZC^T)$ to the global row index, if one is using metis.
  