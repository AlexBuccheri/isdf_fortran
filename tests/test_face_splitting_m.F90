program test_face_splitting_m
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib
    use fortuno_interface_m, only: execute_cmd_app, test, check, is_equal
    use maths_m,             only: all_close
    use face_splitting_m,    only: face_splitting_product_two_funcs_rowwise, &
                                   face_splitting_product_two_funcs_colwise, &
                                   face_splitting_product
    implicit none

    call execute_cmd_app(testitems=[&
            test("Face-splitting of two products", test_face_splitting_product_two_funcs), &
            test("Face-splitting of the product of one matrix", test_face_splitting_product_one_func) &
        ])

contains

    subroutine test_face_splitting_product_two_funcs()
        real(dp) :: c(3, 3), d(3, 3), expected_z(3, 9)
        !integer :: i
        real(dp), allocatable :: z(:, :)

        ! Example matrix from wikipedia:
        ! https://en.wikipedia.org/wiki/Khatri–Rao_product#Face-splitting_product
        ! Transposed so memory layout matches visual layout
        c = transpose(reshape(&
            [1, 2, 3, &
             4, 5, 6, &
             7, 8, 9], [3,3]))

        d = transpose(reshape(&
            [1, 4, 7, &
             2, 5, 8, &
             3, 6, 9], [3,3]))

        expected_z = transpose(reshape(&
             [ 1,  4,  7,  2,  8, 14,  3, 12, 21, &
               8, 20, 32, 10, 25, 40, 12, 30, 48, &
              21, 42, 63, 24, 48, 72, 27, 54, 81], [9, 3]))

        ! First routine
        call face_splitting_product_two_funcs_rowwise(c, d, z)
        call check(size(z, 1) == 3, msg='N rows of Z should be consistent with N rows of both c and d')
        call check(size(z, 2) == 9, msg='N cols of Z should equal product of c and d N cols')
        call check(all_close(z, expected_z), msg='Expect Z to be consistent with reference')

        ! For visual inspection
        ! write(*, *) 'Expected'
        ! do i = 1, 3
        !     write(*, *) expected_z(i, :)
        ! enddo

        ! write(*, *) 'Obtained'
        ! do i = 1, 3
        !     write(*, *) z(i, :)
        ! enddo

        deallocate(z)

        ! Second routine. Expects the inputs of shape (n_states, n_points) and (m_states, n_points), respectively
        ! Performs the kronecker product on the first dimension
        call face_splitting_product_two_funcs_colwise(transpose(c), z, transpose(d))
        call check(size(z, 1) == 9, msg='Nnmber of rows of Z = product of N rows of c and d')
        call check(size(z, 2) == 3, msg='Number of columns of Z consistent with N cols of c and d')
        ! As such, expect the returned Z to be the transpose of what would expect if performing the kronecker
        ! product on the second dimension.
        call check(all_close(z, transpose(expected_z)), msg='Expect Z return as the transpose of the prior routine')

    end subroutine test_face_splitting_product_two_funcs


    subroutine test_face_splitting_product_one_func()
        real(dp) :: c(3, 3)
        real(dp), allocatable :: z1(:, :), z2(:, :)

        ! Example matrix from wikipedia:
        ! https://en.wikipedia.org/wiki/Khatri–Rao_product#Face-splitting_product
        ! Transposed so memory layout matches visual layout
        c = transpose(reshape(&
            [1, 2, 3, &
             4, 5, 6, &
             7, 8, 9], [3,3]))

        call face_splitting_product(c, z1)
        call face_splitting_product(c, c, z2)

        call check(size(z1, 1) == 3, msg='N rows of Z should be consistent with N rows of c')
        call check(size(z1, 2) == 9, msg='N cols of Z should equal (N cols)^2 of c')
        call check(all_close(z1, z2), msg='Expect Z to be the same from both routines')
        deallocate(z1)
        deallocate(z2)

    end subroutine test_face_splitting_product_one_func


end program test_face_splitting_m
