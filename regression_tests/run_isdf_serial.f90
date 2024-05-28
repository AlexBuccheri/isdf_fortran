! Run the ISDF code in serial
program run_isdf_serial
    use, intrinsic :: iso_fortran_env
    ! External libs
    ! use mpi, only: mpi_init, MPI_Finalize, mpi_comm_world

    ! TODO(Alex) For some reason, init and finalize cannot be found
    use mpi_m,    only: mpi_t  !, mpi_init, MPI_Finalize
    use grids_m,  only: discretise_values_to_grid
    use kmeans_m, only: weighted_kmeans

    ! TODO(Alex) Cannot use maths_m from this file as the module name clashes with 
    ! maths_m imported from the kmeans library. Extremely annoying
    ! Internal libs
    use face_splitting_m,       only: face_splitting_product
    use parse,                  only: write_to_xyz, output_cube, parse_grid_1d_from_c, parse_grid_2d_from_c
    use sampling_m,             only: choose_initial_centroids_simple
    use interpolation_points_m, only: subsample_transpose_product_matrix
    use lapack_drivers_m,       only: qr_decomposition_with_pivot
    use isdf_serial_m,          only: construct_interpolation_vectors, construct_approximation_product_states
    implicit none

    type(mpi_t),  allocatable :: comm
    real(real64), parameter :: bohr_to_ang =  0.529177249_real64
    real(real64), parameter :: ang_to_bohr =  1._real64 / bohr_to_ang

    ! Inputs
    integer, parameter :: n_states = 22 !< 22 occupied for benzene in minimal basis
    character(len=100) :: root = "/Users/alexanderbuccheri/Codes/isdf_fortran"
    real(real64), allocatable :: grid(:, :)      !< (ndim, np)
    real(real64), allocatable :: phi(:, :)       !< (nstates, np)
    real(real64), allocatable :: rho(:)          !< (np)
    real(real64), allocatable :: centroids(:, :) 
    character(len=1)          :: species(12)
    integer                   :: an(12)
    real(real64)              :: atomic_pos(3, 12), spacings(3, 3)

    ! Centroid Generation
    integer                       :: niter, n_centroid, min_seed_size
    real(real64)                  :: centroid_tol = 1.e-6_real64
    integer,          allocatable :: init_centroid_indices(:), centroid_indices(:), seed(:)
    character(len=1), allocatable :: dummy_species(:)
    logical                       :: use_centroids

    ! QR decomposition
    integer,      allocatable :: interpolation_indices(:)
    real(real64), allocatable :: zt_subspace(:, :)

    ! Interpolation vectors
    real(real64), allocatable :: theta(:, :), product_exact(:, :), product_isdf(:, :)
    real(real64), allocatable :: error(:)
    real(real64)              :: dv 
    character(len=3)          :: ij_char

    integer :: i, np, ir, ic, ij, imax, imin

    ! -----------------------------
    ! Start of code
    ! -----------------------------

    ! call mpi_init(ierr)
    ! comm = mpi_t(mpi_comm_world)
    comm = mpi_t()

    use_centroids = .true.

    ! benzene
    species = ['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']
    an =      [ 6,   6,   6,   6,   6,   6,   1,   1,   1,   1,   1,   1 ]

    ! Angstrom
    atomic_pos = reshape(&
                  [-0.65914,  -1.21034,  3.98683, & ! Atom 1, accessed as atomic_pos(:, 1)
                    0.73798,  -1.21034,  4.02059, &
                   -1.35771,  -0.00006,  3.96990, &
                    1.43653,  -0.00004,  4.03741, &
                   -0.65915,   1.21024,  3.98685, &
                    0.73797,   1.21024,  4.02061, &
                   -1.20447,  -2.15520,  3.97369, &
                    1.28332,  -2.15517,  4.03382, &
                   -2.44839,  -0.00006,  3.94342, &
                    2.52722,  -0.00004,  4.06369, &
                   -1.20448,   2.15509,  3.97373, &
                    1.28330,   2.15508,  4.03386], [3, 12])

    ! Parse grid and density from benzene example
    ! Units of Bohr and (e/Bohr^3), respectively
    call parse_grid_2d_from_c(trim(root) // "/regression_tests/input/grid.out", .true., grid)
    call parse_grid_1d_from_c(trim(root) // "/regression_tests/input/density.out", rho)

    ! Output to visualise
    np = size(grid, 2)
    allocate(dummy_species(np), source='P')
    call write_to_xyz('grid.xyz', dummy_species, grid * bohr_to_ang)
    deallocate(dummy_species)

    ! Order of fastest index converted when reading from file
    spacings = 0._real64
    spacings(1, 1) = grid(1, 2)  - grid(1, 1) 
    spacings(2, 2) = grid(2, 11) - grid(2, 10) 
    spacings(3, 3) = grid(3, 101) - grid(3, 100) 
    dv = spacings(1, 1) * spacings(2, 2) * spacings(3, 3)
    call output_cube('density', an, atomic_pos * ang_to_bohr, [10, 10, 10], spacings, grid(:, 1), rho)

    ! --------------------------------------------------------------------
    ! Compute centroids 
    ! --------------------------------------------------------------------
    if (use_centroids) then
        niter = 100
        ! 32 = sqrt(np) 
        n_centroid = 32   
        allocate(centroids(3, n_centroid), init_centroid_indices(n_centroid))

        ! Set seed to fix random number generation, for testing    
        ! TODO(Alex) Consider other ways of sampling a non-uniform distribution
        ! a) Implement myself
        ! b) Interface to GSL C code (note, not MPI-enabled)                                                              
        call random_seed(size=min_seed_size)
        allocate(seed(min_seed_size), source=[(i, i=1, min_seed_size)])
        call choose_initial_centroids_simple(grid, rho, centroids, fixed_seed=seed)

        ! Output initial centroids
        allocate(dummy_species(n_centroid), source='P')
        call write_to_xyz('initial_centroids.xyz', dummy_species, centroids * bohr_to_ang)
        deallocate(dummy_species)

        write(*,*) 'Computing centroids with kmeans'
        call weighted_kmeans(comm, grid, rho, centroids, niter, centroid_tol, verbose=.false.)
        deallocate(rho)
        allocate(centroid_indices(n_centroid))
        call discretise_values_to_grid(centroids, grid, centroid_indices)

        ! Output final centroids for visualisation
        allocate(dummy_species(n_centroid), source='P')
        call write_to_xyz('final_centroids.xyz', dummy_species, centroids * bohr_to_ang)
        deallocate(dummy_species)

        ! Might be useful to code some measure of centroid choice.
        deallocate(centroids)
    endif
    ! --------------------------------------------------------------------
    ! Compute optimal interpolation points with QR decomposition 
    ! --------------------------------------------------------------------

    if (.not. use_centroids) then
        ! Parse wave functions in packed form (n_states, np)
        call parse_grid_2d_from_c(trim(root) // "/regression_tests/input/wfs.out", .true., phi)

        write(*, *) 'Compute optimal interpolation points using QR decomposition'
        ! 32 = sqrt(np) 
        n_centroid = 32   

        ! Fix for testing
        call random_seed(size=min_seed_size)
        allocate(seed(min_seed_size), source=[(i, i=1, min_seed_size)])

        call subsample_transpose_product_matrix(phi, n_centroid, zt_subspace, random_seed=seed)
        allocate(interpolation_indices(n_centroid))
        call qr_decomposition_with_pivot(zt_subspace, interpolation_indices, preserve_A=.false.)
        deallocate(zt_subspace)

        ! Interpolation points
        allocate(centroids(3, n_centroid))
        do ic = 1, n_centroid
          ir = interpolation_indices(ic)
          centroids(:, ic) = grid(:, ir)
        enddo

        ! Output final centroids for visualisation
        allocate(dummy_species(n_centroid), source='P')
        call write_to_xyz('qr_final_centroids.xyz', dummy_species, centroids * bohr_to_ang)
        deallocate(dummy_species)
        deallocate(centroids)
        deallocate(phi)
    endif

    ! --------------------------------------------------------------------
    ! Compute interpolation vectors
    ! --------------------------------------------------------------------

    ! Parse wave functions in packed form (np, n_states)
    call parse_grid_2d_from_c(trim(root) // "/regression_tests/input/wfs.out", .false., phi)
    if (all(shape(phi) /= [np, n_states])) then
      write(*, *) 'Expected phi in unpacked form for use in ISDF routines'
      error stop
    endif

    ! ! Use interpolation vectors to expand the wave functions
    write(*, *) 'Computing ISDF vectors'
    call construct_interpolation_vectors(phi, centroid_indices, theta)

    if (all(shape(theta) /= [np, n_centroid])) then
      write(*, *) "Error with shape of returned Theta"
      error stop
    endif

    write(*, *) 'Constructing approximate product of KS states using ISDF vectors'
    call construct_approximation_product_states(theta, phi, phi, centroid_indices, product_isdf)
    if (all(shape(product_isdf) /= [np, n_states * n_states])) then
      write(*, *) "Error with shape of returned product_isdf"
      error stop
    endif
    deallocate(theta)
    deallocate(centroid_indices)

    ! Compute face-splitting product
    write(*, *) 'Constructing exact product of KS functions using face-splitting product'
    call face_splitting_product(phi, product_exact)
    if (all(shape(product_exact) /= [np, n_states * n_states])) then
      write(*, *) "Error with shape of returned product_exact"
      error stop
    endif
    deallocate(phi)

    ! Output both for plotting
    write(ij_char, '(I3)') n_states * n_states
    write(*, *) 'Outputting '// ij_char //' product functions in cube format'

    do ij = 1, n_states * n_states
      write(ij_char, '(I3)') ij 
      call output_cube('product_cubes/product_isdf_'//trim(adjustl(ij_char)), an, atomic_pos * ang_to_bohr, [10, 10, 10], &
        spacings, grid(:, 1), product_isdf(:, ij))
      call output_cube('product_cubes/product_exact'//trim(adjustl(ij_char)), an, atomic_pos * ang_to_bohr, [10, 10, 10], &
        spacings, grid(:, 1), product_exact(:, ij))
    enddo

    ! Compare both numerically
    ! RMSE
    allocate(error(n_states * n_states))
    do ij = 1, n_states * n_states
      error(ij) = sqrt(mean_square_error(product_isdf(:, ij), product_exact(:, ij), dv))
    enddo

    imax = maxloc(error, dim=1)
    imin = minloc(error, dim=1)
    write(*, *) 'Mean RMSE: ', sum(error) / real(n_states * n_states, real64)
    write(*, *) 'Min RMSE found for state :', imin, error(imin)
    write(*, *) 'Max RMSE found for state :', imax, error(imax)

    ! Percentage error
    error = 0._real64
    do ij = 1, n_states * n_states
      error(ij) = relative_error(product_isdf(:, ij), product_exact(:, ij))
    enddo

    imax = maxloc(error, dim=1)
    imin = minloc(error, dim=1)
    write(*, *) 'Mean relative error: ', sum(error) / real(n_states * n_states, real64)
    write(*, *) 'Min relative error found for state :', imin, error(imin)
    write(*, *) 'Max relative error found for state :', imax, error(imax)

    ! call MPI_Finalize(ierr)
    deallocate(grid)
    deallocate(product_isdf)
    deallocate(product_exact)
    deallocate(error)

contains

function relative_error(x_ref, x) result(rel_err)
  real(real64), intent(in) :: x_ref(:)
  real(real64), intent(in) :: x(:)
  real(real64) :: rel_err

  integer :: n, i

  n = size(x_ref)
  if (n /= size(x)) then
     write(*, *) 'Size of x_ref and x disagree'
     error stop 
  endif

  rel_err = 0._real64
  do i = 1, n
    rel_err = rel_err + abs((x_ref(i) - x(i)) / x_ref(i))
  enddo

end function relative_error

function mean_square_error(x, y, dv) result(mse)
  real(real64), intent(in) :: x(:)
  real(real64), intent(in) :: y(:)
  real(real64), intent(in) :: dv   !< Volume element
  real(real64) :: mse

  integer :: n, i

  n = size(x)
  if (n /= size(y)) then
     write(*, *) 'Size of x and y disagree'
     error stop 
  endif

  mse = 0._real64

  !$omp parallel do simd default(shared) reduction(+:mse)
  do i = 1, n
    mse = mse + (x(i) - y(i))**2._real64
  enddo
  !$omp end parallel do simd

  mse = mse / dv

end function mean_square_error

end program run_isdf_serial
