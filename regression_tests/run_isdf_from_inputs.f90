

! Run ISDF serial tests using interpolation points defined on input,
! in the form of grid indices
program run_isdf_from_inputs
    use, intrinsic :: iso_fortran_env
    use omp_lib

    use face_splitting_m,       only: face_splitting_product
    use parse,                  only: write_to_xyz, output_cube, parse_grid_1d_from_c, parse_grid_2d_from_c
    use isdf_serial_m,          only: construct_interpolation_vectors, construct_approximation_product_states
    implicit none

    ! Inputs
    character(len=100)      :: root = "/Users/alexanderbuccheri/Codes/isdf_fortran"
    real(real64), parameter :: bohr_to_ang =  0.529177249_real64
    real(real64), parameter :: ang_to_bohr =  1._real64 / bohr_to_ang
    integer, parameter      :: n_states = 21                         !< 21 occupied for benzene in minimal basis
    character(len=1)        :: species(12)
    integer                 :: an(12)
    real(real64)            :: atomic_pos(3, 12)
    character(len=1), allocatable :: dummy_species(:)

    ! Grids, wave functions and interpolation points
    real(real64), allocatable :: grid(:, :)      !< (ndim, np)
    real(real64)              :: spacings(3, 3)
    real(real64), allocatable :: phi(:, :)       !< (nstates, np) or (np, nstates) depending on choice
    integer                   :: niter, np
    integer,      allocatable :: interpolation_indices(:, :), interpolation_indices_packed(:)
    integer                   :: limits(3)       !< Hard-coded limits of each grid dimension

    ! For visualisation
    real(real64), allocatable :: centroids(:, :), rho(:)
    logical :: visualise

    ! Interpolation vectors
    real(real64), allocatable :: theta(:, :), product_exact(:, :), product_isdf(:, :)
    real(real64), allocatable :: error(:)
    character(len=3)          :: ij_char, char_niter
    integer :: n_product

    integer :: i, ij, imax, imin, ic, ir

    ! --------------------
    ! Main code
    ! --------------------
    ! Valid choices = 30, 50, 70, 90
    niter = 30
    visualise = .false.

    write(*, *) 'Using ', niter,'interpolation points'

    ! benzene
    species = ['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']
    an =      [ 6,   6,   6,   6,   6,   6,   1,   1,   1,   1,   1,   1 ]

    ! Angstrom (from .xyz)
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

    ! Parse grid from benzene example. Units of Bohr
    call parse_grid_2d_from_c(trim(root) // "/regression_tests/input/grid.out", .true., grid)

    ! Definitions required for cube file
    limits = [10, 10, 10]
    np = size(grid, 2)
    spacings = 0._real64
    spacings(1, 1) = grid(1, 2)  - grid(1, 1) 
    spacings(2, 2) = grid(2, 11) - grid(2, 10) 
    spacings(3, 3) = grid(3, 101) - grid(3, 100) 

    ! Load centroid indices
    write(*, *) 'Load interpolation indices'
    allocate(interpolation_indices(3, niter))
    write(char_niter,  "(i3)") niter
    call parse_centroid_indices(trim(root) // "/inputs/centroid_indices/centroid_indices_" // trim(adjustl(char_niter)) //".txt", &
        interpolation_indices)

    ! Convert from three to one index
    write(*, *) 'Map interpolation indices form (ix, iy, iz) to ir'
    allocate(interpolation_indices_packed(niter))
    call map_indices_to_comp(interpolation_indices, limits, interpolation_indices_packed)
    deallocate(interpolation_indices)    

    if (visualise) then
      ! Visualise centroid points with the density and grid, to confirm the conversions were done correctly
      ! Grid 
      write(*, *) 'Output grid to .xyz'
      allocate(dummy_species(np), source='P')
      call write_to_xyz('grid.xyz', dummy_species, grid * bohr_to_ang)
      deallocate(dummy_species)

      ! Centroids
      write(*, *) 'Output centroids to .xyz'
      allocate(dummy_species(niter), source='P')
      allocate(centroids(3, niter))
      do ic = 1, niter
        ir = interpolation_indices_packed(ic)
        centroids(:, ic) = grid(:, ir)
      enddo
      call write_to_xyz('centroids.xyz', dummy_species, centroids * bohr_to_ang)
      deallocate(dummy_species)
      deallocate(centroids)

      ! Density
      write(*, *) 'Output density to .xyz'
      call parse_grid_1d_from_c(trim(root) // "/regression_tests/input/density.out", rho)
      call output_cube('density', an, atomic_pos * ang_to_bohr, limits, spacings, grid(:, 1), rho)
      deallocate(rho)
    endif

    ! --------------------------------------------------------------------
    ! Compute interpolation vectors
    ! --------------------------------------------------------------------
    
    ! Parse wave functions in form (np, n_states)
    call parse_grid_2d_from_c(trim(root) // "/regression_tests/input/wfs.out", .false., phi)

    ! TODO(Alex) Note, this did not catch when I was using 22 states - somewhat weird
    if (all(shape(phi) /= [np, n_states])) then
      write(*, *) 'Expected phi in unpacked form for use in ISDF routines'
      error stop
    endif

    ! Output KS states as cube files
    do i = 1, n_states
      write(ij_char, '(I3)') i 
      call output_cube('ks_states/benzene_wf_'//trim(adjustl(ij_char)), an, atomic_pos * ang_to_bohr, limits, &
        spacings, grid(:, 1), phi(:, i))
    enddo

    ! Use interpolation vectors to expand the wave functions
    write(*, *) 'Computing ISDF vectors'
    call construct_interpolation_vectors(phi, interpolation_indices_packed, theta)

    if (all(shape(theta) /= [np, niter])) then
      write(*, *) "Error with shape of returned Theta"
      error stop
    endif

    n_product = n_states * n_states

    write(*, *) 'Constructing approximate product of KS states using ISDF vectors'
    call construct_approximation_product_states(theta, phi, phi, interpolation_indices_packed, product_isdf)
    if (all(shape(product_isdf) /= [np, n_product])) then
      write(*, *) "Error with shape of returned product_isdf"
      error stop
    endif

    write(*, *) 'Constructing exact product of KS functions using face-splitting product'
    call face_splitting_product(phi, product_exact)
    if (all(shape(product_exact) /= [np, n_product])) then
      write(*, *) "Error with shape of returned product_exact"
      error stop
    endif
    deallocate(phi)
    

    ! Output both for plotting
    write(ij_char, '(I3)') n_states * n_states
    write(*, *) 'Outputting '// ij_char //' product functions in cube format'

    do ij = 1, n_states * n_states
      write(ij_char, '(I3)') ij 
      call output_cube('product_cubes/product_isdf_'//trim(adjustl(ij_char)), an, atomic_pos * ang_to_bohr, limits, &
        spacings, grid(:, 1), product_isdf(:, ij))
      call output_cube('product_cubes/product_exact'//trim(adjustl(ij_char)), an, atomic_pos * ang_to_bohr, limits, &
        spacings, grid(:, 1), product_exact(:, ij))
    enddo


    ! Quantify error for each product state
    allocate(error(n_product))
    do ij = 1, n_product
      error(ij) = mean_square_error_on_regular_grid(product_isdf(:, ij), product_exact(:, ij))
      write(*, *) 'MSE error:', ij, error(ij)
    enddo

    imax = maxloc(error, dim=1)
    imin = minloc(error, dim=1)
    write(*, *) '(Min, Mean, Max) MSEs: ', error(imin), sum(error) / real(n_product, real64), error(imax)

    ! TODO ADD Output for plotting

    deallocate(theta)
    deallocate(grid)

contains

    !> Parse centroid indices and convert from 0-based indexing to 1-based indexing
    !!
    !! Note, another way to have done this would be to read in the actual centroid points
    !! then map to the closest grid point (which should be exact)
    subroutine parse_centroid_indices(fname, indices)
        character(len=*), intent(in) :: fname
        integer, intent(out) ::  indices(:, :)  !< shape (n_dim, n_interpolation_points)
        integer :: unit = 101
        integer :: i, ierr

        open(unit=unit, file=trim(fname), form='formatted', access='stream', iostat=ierr)

        if (ierr /= 0) then
            write(*, *) "Error opening file: ", trim(fname)
            stop
        end if

        ! Skip the header line (although one could get n interpolation points from here)
        read(unit, *) 

        do i = 1, size(indices, 2)
            read(unit, *) indices(1:3, i)
            ! Convert to 1-based indexing
            indices(1:3, i) = indices(1:3, i) + 1
        enddo

        close(unit)

    end subroutine parse_centroid_indices

    ! Program Hello
    !     integer :: ix, iy, iz, ir, f
    !     integer :: nx, ny, nz
        
        
    !     nx = 2
    !     ny = 3
    !     nz = 4
        
    !     ir = 0
    !     do iz = 1, nz
    !         do iy = 1, ny
    !             do ix = 1, nx
    !                 f = ix + (iy - 1) * nx + (iz - 1) * nx * ny
    !                 ir = ir + 1
    !                 write(*, *) ir, f
    !             enddo
    !         enddo
    !     enddo
        
    
    ! End Program Hello

    subroutine map_indices_to_comp(indices, limits, cmp_indices)
        integer, intent(in) ::  indices(:, :)
        integer, intent(in) ::  limits(:)
        integer, intent(out) :: cmp_indices(:)

        integer :: iindex, nx, ny, nz, ix, iy, iz

        nx = limits(1)
        ny = limits(2)
        nz = limits(3)

        do iindex = 1, size(cmp_indices)
            ix = indices(1, iindex)
            iy = indices(2, iindex)
            iz = indices(3, iindex)
            cmp_indices(iindex) = ix + (iy - 1) * nx + (iz - 1) * nx * ny
        enddo

    end subroutine map_indices_to_comp



    function mean_square_error_on_regular_grid(x, y) result(mse)
        real(real64), intent(in) :: x(:)
        real(real64), intent(in) :: y(:)
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
      
        mse = mse / real(n, real64)
      
      end function mean_square_error_on_regular_grid

end program run_isdf_from_inputs
