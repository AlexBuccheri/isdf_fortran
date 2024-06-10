
program isdf_components
    use, intrinsic :: iso_fortran_env, only: real64
    use helpers,       only: load_numpy_array, write_numpy_array, load_centroid_indices
    use isdf_serial_m, only: construct_zct, select_interpolation_points_of_first_dim, &
        construct_inverse_coefficients_contraction, construct_quasi_density_matrix_R_intr, &
        construct_quasi_density_matrix_R_intr_alt, construct_approximation_product_states, &
        construct_interpolation_vectors
    use lapack_drivers_m, only: pseudo_inv, inversion_lu
    use parse,            only: output_cube
    implicit none

    ! Reference data
    character(len=59) :: ref_root = "/Users/alexanderbuccheri/Codes/isdf_prototypes/mocked_data/"
    character(len=60) :: fortran_data_root = "/Users/alexanderbuccheri/Codes/isdf_prototypes/fortran_data/"

    integer, parameter :: limits(3) = [10, 10, 10]      !< Grid limits per dimension
    integer, parameter :: n_total = 1000                !< Number of grid points for benzene PYSCF example
    integer, parameter :: n_states = 21                 !< Number of occupied states benzene PYSCF example
    integer, parameter :: n_centroids = 60              !< Number of interpolation points

    ! System details required for cube output
    real(real64), parameter :: bohr_to_ang =  0.529177249_real64
    real(real64), parameter :: ang_to_bohr =  1._real64 / bohr_to_ang
    integer, parameter :: n_atoms = 12
    integer            :: an(n_atoms)                   !< Atomic numbers
    real(real64)       :: atomic_pos(3, n_atoms)        !< Atomic positions
    real(real64)       :: spacings(3, 3)                !< Grid spacings
    real(real64)       :: origin(3)                     !< First grid point


    ! Data
    real(real64), allocatable :: phi(:, :)              !< Batch of wave functions
    real(real64), allocatable :: P_rr(:, :)             !< Density matrix defined for all grid points
    real(real64), allocatable :: P_rr_alt(:, :)         !< Density matrix defined for all grid points
    real(real64), allocatable :: P_r_mu(:, :)           !< Density matrix defined for (grid, intrp)
    real(real64), allocatable :: P_r_mu_alt(:, :)       !< Density matrix defined for (grid, intrp)
    real(real64), allocatable :: P_mu_nu(:, :)          !< Density matrix defined for (intrp, intrp)
    real(real64), allocatable :: zct(:, :)              !< Z and C^T product
    real(real64), allocatable :: cct(:, :)              !< C and C^T product
    real(real64), allocatable :: cct_alt(:, :)          !< C and C^T product constructed by sampling ZC^T
    real(real64), allocatable :: cct_inv(:, :)          !< Inverse of C and C^T product
    real(real64), allocatable :: theta(:, :)            !< Contract of (ZC^T) and (CC^T)^{-1}
    real(real64), allocatable :: z_isdf(:, :)           !< Expansion of product functions in ISDF basis

    integer,      allocatable :: centroid_indices(:)    !< Centroid indices as composite index
    character(len=3)          :: char_nc, ij_char       !< Convert n_centroids to char

    integer :: ir, jr, j_mu, ij

    ! -------------------------
    ! Main routine
    ! -------------------------
    ! Read batched wave functions and compute P_rr
    write(*, *) 'Load phi.txt'
    call load_numpy_array(ref_root // 'phi.txt', phi)
    if (all(shape(phi) /= [n_total, n_states])) then
        write(*, *) 'Shape of phi is wrong', shape(phi)
        error stop
    endif

    ! Read density matrix contracted over states {phi}
    ! Assume all density matrix quantities are for phi, so drop that suffix
    write(*, *) 'Load P_rr.txt'
    call load_numpy_array(ref_root // 'P_rr.txt', P_rr)
    if (all(shape(P_rr) /= [n_total, n_total])) then
        write(*, *) 'Shape of P_rr is wrong', shape(P_rr)
        error stop
    endif

    ! Compute P_rr by contracting over the state index, P = phi @ phi^T
    allocate(P_rr_alt(n_total, n_total))
    call dgemm ('N', 'T',   &
        n_total,              & ! size(op(A), 1)
        n_total,              & ! size(C, 2)
        n_states,        & ! size(op(A), 2)
        1._real64,       & ! alpha
        phi, n_total,         & ! A and LDA
        phi, n_total,         & ! B and LDB
        0._real64,       & ! beta
        P_rr_alt, n_total     & ! C and LDC
        )
    call write_numpy_array(fortran_data_root// "f90_P_rr_alt.txt", P_rr_alt)

    ! Read centroid indices
    write(*, *) 'Load centroid indices'
    write(char_nc, '(I3)') n_centroids
    call load_centroid_indices(ref_root // "centroid_indices_"//trim(adjustl(char_nc))//".txt", &
        centroid_indices)
    if (size(centroid_indices) /= n_centroids) then
        write(*, *) 'Shape of centroid_indices is wrong', size(centroid_indices)
        error stop
    endif

    ! Construct P_r_mu from imported phi, rather than imported P_rr
    call construct_quasi_density_matrix_R_intr(phi, centroid_indices, P_r_mu_alt)
    call write_numpy_array(fortran_data_root// "f90_P_r_mu_alt.txt", P_r_mu_alt)

    deallocate(P_r_mu_alt)
    call construct_quasi_density_matrix_R_intr_alt(phi, centroid_indices, P_r_mu_alt)
    call write_numpy_array(fortran_data_root// "f90_P_r_mu_alt_dgemm.txt", P_r_mu_alt)

    ! Construct P_r_mu from masking P_rr, and export
    write(*, *) 'Construct P_r_mu and export'
    allocate(P_r_mu(n_total, n_centroids))
    do j_mu = 1, n_centroids
        jr = centroid_indices(j_mu)
        do ir = 1, n_total
            P_r_mu(ir, j_mu) = P_rr(ir, jr)
        enddo
    enddo
    call write_numpy_array(fortran_data_root// "f90_P_r_mu.txt", P_r_mu)

    ! Construct and export zct
    write(*, *) 'Construct ZC^T and export'
    allocate(zct(n_total, n_centroids))
    call construct_zct(P_r_mu, zct)
    call write_numpy_array(fortran_data_root// "f90_zct.txt", zct)

    ! Construct P_mu_nu
    write(*, *) 'Construct P_mu_nu and export'
    allocate(P_mu_nu(n_centroids, n_centroids))
    call select_interpolation_points_of_first_dim(P_r_mu, centroid_indices, P_mu_nu)
    call write_numpy_array(fortran_data_root// "f90_P_mu_nu.txt", P_mu_nu)

    ! Compute CC^T, and inv(CC^T)
    write(*, *) 'Construct CC^T and export'
    allocate(cct(n_centroids, n_centroids))
    call construct_inverse_coefficients_contraction(P_mu_nu, P_mu_nu, cct)
    call write_numpy_array(fortran_data_root// "f90_cct.txt", cct)

    ! Also test computing CC^T by selecting rows of ZC^T according to centroid_indices
    write(*, *) 'Construct CC^T by selecting rows of ZC^T'
    allocate(cct_alt(n_centroids, n_centroids))
    call select_interpolation_points_of_first_dim(zct, centroid_indices, cct_alt)
    call write_numpy_array(fortran_data_root// "f90_cct_alt.txt", cct_alt)

    ! Inverse of (CC^T)
    ! NOTE: inversion appears to be a problem - Try LU instead
    write(*, *) 'Compute inv(CC^T)'
    allocate(cct_inv(n_centroids, n_centroids))
    ! call pseudo_inv(cct, cct_inv)
    call inversion_lu(cct, cct_inv)
    call write_numpy_array(fortran_data_root// "f90_cct_inv.txt", cct)

    ! a) Allow comparison to continue
    write(*, *) 'Load cct_inv.txt to allow code comparison to continue'
    deallocate(cct_inv)
    call load_numpy_array(ref_root // 'cct_inv.txt', cct_inv)
    if (all(shape(cct_inv) /= [n_centroids, n_centroids])) then
        write(*, *) 'Shape of imported (CC^T)^{-1} is wrong', shape(cct_inv)
        error stop
    endif

    ! Compute Theta = (ZC^T) @ (CC^T)^{-1}
    write(*, *) 'Compute and export Theta'
    allocate(theta(n_total, n_centroids))
    call dgemm ('N', 'N', &
        n_total, &               ! row(ZCT)
        n_centroids, &           ! col(cct_inv)
        n_centroids, &           ! col(ZCT)
        1._real64, &             ! alpha
        zct, n_total, &          !LDA
        cct_inv, n_centroids, &  ! LDB
        0._real64, &             ! beta
        theta, n_total &         ! LDC
        )
    ! Alternatively, use matmul
    !theta = matmul(zct, cct_inv)
    call write_numpy_array(fortran_data_root// "f90_theta.txt", theta)

    ! Approximate ISDF product functions
    write(*, *) 'Compute and export approximate product functions'
    call construct_approximation_product_states(theta, phi, phi, centroid_indices, z_isdf)
    if (all(shape(z_isdf) /= [n_total, n_states * n_states])) then
        write(*, *) 'Shape of z_isdf is wrong', shape(z_isdf)
        error stop
    endif
    call write_numpy_array(fortran_data_root// "f90_z_isdf.txt", z_isdf)

    ! Output for plotting
    ! NOTE: This won't give the correct result because it assumes the data array(ir) -> array(ix, iy, iz)
    ! whereas it was constructed in python, so will go array(iz, iy, ix)
    ! write(ij_char, '(I3)') n_states * n_states
    ! write(*, *) 'Outputting '// ij_char //' approximate product functions in cube format'

    ! ! Benzene
    ! an = [ 6,   6,   6,   6,   6,   6,   1,   1,   1,   1,   1,   1 ]

    ! ! Atomic positions in angstrom (from .xyz)
    ! atomic_pos = reshape(&
    !               [-0.65914,  -1.21034,  3.98683, &
    !                 0.73798,  -1.21034,  4.02059, &
    !                -1.35771,  -0.00006,  3.96990, &
    !                 1.43653,  -0.00004,  4.03741, &
    !                -0.65915,   1.21024,  3.98685, &
    !                 0.73797,   1.21024,  4.02061, &
    !                -1.20447,  -2.15520,  3.97369, &
    !                 1.28332,  -2.15517,  4.03382, &
    !                -2.44839,  -0.00006,  3.94342, &
    !                 2.52722,  -0.00004,  4.06369, &
    !                -1.20448,   2.15509,  3.97373, &
    !                 1.28330,   2.15508,  4.03386], [3, 12])
    ! atomic_pos = atomic_pos * ang_to_bohr

    ! ! Grid spacings and first grid point, in Bohr
    ! spacings = 0._real64
    ! spacings(1, 1) = 1.71139336_real64
    ! spacings(2, 2) = 1.5716964_real64
    ! spacings(3, 3) = 0.69191971_real64
    ! origin = [-7.62678655_real64, -7.07273774_real64, 4.45198379_real64] 

    ! do ij = 1, n_states * n_states
    !   write(ij_char, '(I3)') ij 
    !   call output_cube(fortran_data_root//'z_isdf/z_isdf_'//trim(adjustl(ij_char)), &
    !     an,          &
    !     atomic_pos,  &
    !     limits,      &
    !     spacings,    &
    !     origin,      &
    !     z_isdf(:, ij)&
    !     )
    ! enddo

    ! Test high-level routine for generating theta
    write(*, *) 'Compute and export theta using high-level routine. Expect inversion to be an issue'
    deallocate(theta)
    call construct_interpolation_vectors(phi, centroid_indices, theta)
    call write_numpy_array(fortran_data_root// "f90_theta_alt.txt", theta)



end program isdf_components
