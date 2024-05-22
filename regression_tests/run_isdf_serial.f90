! Run the ISDF code in serial
program run_isdf_serial
    use, intrinsic :: iso_fortran_env
    ! External libs
    ! use mpi, only: mpi_init, MPI_Finalize, mpi_comm_world
    ! TODO(Alex) For some reason, init and finalize cannot be found
    use mpi_m,    only: mpi_t  !, mpi_init, MPI_Finalize
    use grids_m,  only: discretise_values_to_grid
    use kmeans_m, only: weighted_kmeans
    ! Internal libs
    use parse,       only: write_to_xyz, output_cube, parse_grid_from_c
    use sampling_m,    only: choose_initial_centroids_simple
    implicit none

    character(len=100) :: root = "/Users/alexanderbuccheri/Codes/isdf_fortran"
    real(real64), allocatable :: grid(:, :)      !< (ndim, np)
    real(real64), allocatable :: wfs(:, :)       !< (nstates, np)
    real(real64), allocatable :: rho(:)          !< (np)
    real(real64), allocatable :: centroids(:, :) 
    type(mpi_t),  allocatable :: comm
    character(len=1)          :: species(12)
    integer :: an(12)
    real(real64) :: atomic_pos(3, 12), spacings(3, 3)
    real(real64), parameter :: bohr_to_ang =  0.529177249_real64
    real(real64), parameter :: ang_to_bohr =  1._real64 / bohr_to_ang

    ! Centroid Generation
    integer                       :: niter, n_centroid, min_seed_size
    real(real64)                  :: centroid_tol = 1.e-6_real64
    integer,          allocatable :: init_centroid_indices(:), seed(:)
    character(len=1), allocatable :: dummy_species(:)

    integer :: i, ierr, np

    ! call mpi_init(ierr)
    ! comm = mpi_t(mpi_comm_world)
    comm = mpi_t()

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
    call parse_grid_from_c(trim(root) // "/regression_tests/input/grid.out", grid)
    call parse_grid_from_c(trim(root) // "/regression_tests/input/density.out", rho)

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
    call output_cube('density', an, atomic_pos * ang_to_bohr, [10, 10, 10], spacings, grid(:, 1), rho)

    ! --------------------------------------------------------------------
    ! Compute centroids 
    ! --------------------------------------------------------------------
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

    write(*,*)' Computing centroids with kmeans'
    call weighted_kmeans(comm, grid, rho, centroids, niter, centroid_tol, verbose=.true.)
    deallocate(rho)
    call discretise_values_to_grid(centroids, grid)

    ! Output final centroids for visualisation
    allocate(dummy_species(n_centroid), source='P')
    call write_to_xyz('final_centroids.xyz', dummy_species, centroids * bohr_to_ang)
    deallocate(dummy_species)

    ! Might be useful to code some measure of centroid choice.

    ! TODO(Alex) Should ulimately implement QR decomposition so I can compare the two

    ! Compute interpolation vectors
    ! call parse_text(trim(root) // "/regression_tests/input/wfs.out", wfs)

    ! use interpolation vectors to expand the wave functions

    ! Compute face-splitting product

    ! Output both for plotting

    ! Compare both numerically

    ! Consider comparing to python result

    deallocate(grid)
    ! deallocate(wfs)
    deallocate(centroids)
    ! call MPI_Finalize(ierr)

end program run_isdf_serial
