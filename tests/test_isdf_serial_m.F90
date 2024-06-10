program test_isdf_m
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib

    use grids_m, only: linspace, linspace_to_grid
    use fortuno_interface_m, only: execute_cmd_app, test, check, is_equal
    
    ! maths_m is coming from kmeans lib due to lack of sensible namespacing
    use maths_m,             only: all_close
    use parse,               only: output_cube
    use grid_funcs_m,        only: GridParams3D, construct_sin_basis_with_grid
    
    use isdf_serial_m, only: construct_quasi_density_matrix_R_intr, construct_quasi_density_matrix_R_intr_alt, &
        select_interpolation_points_of_first_dim

    implicit none

    ! Register tests
    call execute_cmd_app(testitems=[&
            test("Test quasi-density matrix", test_construct_quasi_density_matrix_R_intr) &
        ])


contains

    ! DONE:
    ! Tested: construct_quasi_density_matrix_R_intr against dgemm implementation. Checked maths.
    ! Tested: select_interpolation_points_of_first_dim. Checked maths.
    ! Tested: construct_inverse_coefficients_contraction indirectly, by testing matrix inversion. Checked maths

    ! TODO
    ! Test: ZCT. 
    !      Validate by inspection. Or compare to doing P_phi_ri * P_phi_ri .. which is equivalent but less efficient.
    ! Contraction of ZCT and inv_CCT. Validate by inspection
    ! Test theta
    !     Not sure how, other than to use it
    ! Test: construct_approximation_product_states
    !   Should converge on exact states, given enough interpolation points that are well-chosen


    ! Create a set of n_states defined on a real-space grid
    subroutine mock_states(grid, sin_integers, phi)
        real(dp), intent(in) :: grid(:, :)  !<(3, np)
        integer,  intent(in) :: sin_integers(:, :) !< (3, n_states)
        
        integer :: np, n_states, i
        real(dp), intent(out) :: phi(:, :)

        np = size(grid, 2)
        n_states = size(sin_integers, 2)

        do i = 1, n_states
            call construct_sin_basis_with_grid(grid, sin_integers(:, i), phi(:, i))
        enddo

    end subroutine mock_states

  
    subroutine test_construct_quasi_density_matrix_R_intr()
        ! Grid
        type(GridParams3D) :: regular
        real(dp) :: origin(3), limits(3), spacings(3, 3), position(3, 1)
        integer  :: np, points(3)
        ! 
        real(dp), allocatable :: phi(:, :)
        integer, parameter :: n_states = 7
        integer :: sin_integers(3, n_states)
        integer, allocatable :: indices(:)
        ! Density matrices
        real(dp), allocatable :: P(:, :), P_alt(:, :), P_mu_nu_ref(:, :), P_mu_nu(:, :)

        integer :: i, j, n_interp, ir

        ! Grid
        points = [10, 10, 10]
        origin = [0.0_dp, 0.0_dp, 0.0_dp]
        limits = [1.0_dp, 1.0_dp, 1.0_dp]
        call regular%init(points, origin, limits)
        np = regular%np

        sin_integers = reshape([1, 1, 1, &
                                1, 1, 2, &
                                1, 1, 3, &
                                2, 2, 2, &
                                1, 1, 4, &
                                1, 2, 3, &
                                1, 2, 4], [3, n_states])

        allocate(phi(np, n_states))
        call mock_states(regular%grid, sin_integers, phi)

        ! For visualisation
        ! position(:, 1) = regular%origin(:)
        !call output_cube("sin_func", [1], position, regular%points, regular%spacings, regular%origin, phi(:, 7))

        ! Define some interpolation points - can be arbtirary
        indices = [451, 114, 800, 687, 238, 99, 572, 869, 330, 714, 185, 382, 660, 961, 529, 493, 81, &
                   211, 780, 405, 923, 644, 30, 565, 756, 152, 256, 854, 366, 44]
        n_interp = size(indices)
    
        call check(maxval(indices) <= np, msg='Max index should not exceed number of grid points')

        call construct_quasi_density_matrix_R_intr(phi, indices, P)

        call construct_quasi_density_matrix_R_intr(phi, indices, P_alt)

        ! Could do k-diff on the files
        ! do j = 1, size(P_alt, 2)
        !     do i = 1, size(P_alt, 1)
        !         write(201, *) P(i, j)
        !         write(202, *) P_alt(i, j)
        !     enddo
        ! enddo

        call check(all_close(P, P_alt), msg='Should be the same for any choice of phi and indices')
        call check(all(shape(P) == [np, n_interp]))

        ! P is more explicitly, P_r_nu
        allocate(P_mu_nu_ref(n_interp, n_interp))
        P_mu_nu_ref = P(indices, :)

        allocate(P_mu_nu(n_interp, n_interp))
        call select_interpolation_points_of_first_dim(P, indices, P_mu_nu)

        call check(all_close(P_mu_nu, P_mu_nu_ref), msg='Should be the same')


    end subroutine test_construct_quasi_density_matrix_R_intr


end program test_isdf_m
