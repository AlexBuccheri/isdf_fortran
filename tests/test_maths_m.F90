program test_maths_m
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib
    use fortuno_interface_m, only: execute_cmd_app, test, check, is_equal
    use maths_m

    implicit none

    call execute_cmd_app(testitems=[&
            test("Gram-Schmidt orthogonalisation", test_gram_schmidt) &
        ])

contains

    ! Test classic Gram Schmidt, and its modified version
    subroutine test_gram_schmidt()
        real(dp), allocatable :: y(:, :), y_ref(:, :)

        ! In each case, one could also test matmul(transpose(y), y) == I

        ! ---------------------------------
        ! 1. Input vectors are orthogonal
        ! ---------------------------------
        allocate(y(3, 3), y_ref(3, 3))
        y = reshape(&
            [1._dp, 0._dp, 0._dp, &
             0._dp, 1._dp, 0._dp, &
             0._dp, 0._dp, 1._dp], [3, 3])
        y_ref = y

        call gram_schmidt(y)
        call check(all_close(y, y_ref), msg='Expect y to be unchanged')

        call modified_gram_schmidt(y)
        call check(all_close(y, y_ref), msg='Expect y to be unchanged')

        ! ---------------------------------
        ! 2. Input vectors are non-orthogonal
        ! ---------------------------------
        ! Reference orthogonalised vectors
        y_ref = reshape(&
            [ 0.26726124191242440_dp,  0.53452248382484879_dp,  0.80178372573727330_dp, & ! column vector 1
             -0.87287156094396956_dp, -0.21821789023599247_dp,  0.43643578047198461_dp, & ! column vector 2
             -0.40824829046386363_dp,  0.81649658092772537_dp, -0.40824829046386363_dp  & ! column vector 3
            ],[3, 3])

        y = reshape(&
            [ 1._dp, 2._dp,  3._dp, &   ! column vector 1
             -1._dp, 0._dp,  1._dp, &   ! column vector 2
              2._dp, 1._dp, -1._dp  &   ! column vector 3
             ],[3, 3])
        call gram_schmidt(y)
        call check(all_close(y, y_ref))

        y = reshape(&
            [ 1._dp, 2._dp,  3._dp, &   ! column vector 1
             -1._dp, 0._dp,  1._dp, &   ! column vector 2
              2._dp, 1._dp, -1._dp  &   ! column vector 3
            ],[3, 3])
        call modified_gram_schmidt(y)
        call check(all_close(y, y_ref), msg='Expect y to be unchanged')


        ! --------------------------------------------------------------
        ! 3. Input vectors are all linearly-dependent
        ! and are not naturally handled by Gram-Schmidt implementations
        ! --------------------------------------------------------------
        y_ref = reshape(&
        [0.40824829046386307_dp, 0.81649658092772615_dp, 0.40824829046386307_dp , & ! column vector 1
         0.0_dp, 0.0_dp, 0.0_dp, & ! column vector 2
         0.0_dp, 0.0_dp, 0.0_dp  & ! column vector 3
        ],[3, 3])

        y = reshape(&
        [ 1._dp,  2._dp,  1._dp, &   ! column vector 1
          2._dp,  4._dp,  2._dp, &   ! column vector 2
         -1._dp, -2._dp, -1._dp  &   ! column vector 3
         ],[3, 3])
        call gram_schmidt(y)
        call check(all_close(y, y_ref), msg=&
            & 'First vector is normalised, the other two are zeroed')

        y = reshape(&
            [ 1._dp,  2._dp,  1._dp, &   ! column vector 1
              2._dp,  4._dp,  2._dp, &   ! column vector 2
             -1._dp, -2._dp, -1._dp  &   ! column vector 3
             ],[3, 3])
        call modified_gram_schmidt(y)
        call check(all_close(y, y_ref), msg=&
            & 'First vector is normalised, the other two are zeroed')

        deallocate(y, y_ref)

    end subroutine test_gram_schmidt

end program test_maths_m
