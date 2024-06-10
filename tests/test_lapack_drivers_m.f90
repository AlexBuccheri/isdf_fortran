program test_lapack_drivers
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortuno_interface_m, only: execute_cmd_app, test, check, is_equal

    ! maths_m is coming from kmeans lib due to lack of sensible namespacing
    use maths_m,             only: all_close

    use lapack_drivers_m, only: pseudo_inv, inversion_lu
    implicit none

    ! Register tests
    call execute_cmd_app(testitems=[&
            test("Test pseudo-inverse", test_pseudo_inverse), &
            test("Test LU inverse", test_inversion_lu) &
        ])

        ! Disagrees amongst numpy, pseudo_inv and LU decom - could just be a problematic choice
        ! a = transpose(reshape([ &
        !       1._dp,  2._dp,  3._dp,  4._dp,  5._dp,  6._dp,  7._dp,  8._dp,  9._dp, 10._dp, &
        !      11._dp, 12._dp, 13._dp, 14._dp, 15._dp, 16._dp, 17._dp, 18._dp, 19._dp, 20._dp, &
        !      21._dp, 22._dp, 23._dp, 24._dp, 25._dp, 26._dp, 27._dp, 28._dp, 29._dp, 30._dp, &
        !      31._dp, 32._dp, 33._dp, 34._dp, 35._dp, 36._dp, 37._dp, 38._dp, 39._dp, 40._dp, &
        !      41._dp, 42._dp, 43._dp, 44._dp, 45._dp, 46._dp, 47._dp, 48._dp, 49._dp, 50._dp, &
        !      51._dp, 52._dp, 53._dp, 54._dp, 55._dp, 56._dp, 57._dp, 58._dp, 59._dp, 60._dp, &
        !      61._dp, 62._dp, 63._dp, 64._dp, 65._dp, 66._dp, 67._dp, 68._dp, 69._dp, 70._dp, &
        !      71._dp, 72._dp, 73._dp, 74._dp, 75._dp, 76._dp, 77._dp, 78._dp, 79._dp, 80._dp, &
        !      81._dp, 82._dp, 83._dp, 84._dp, 85._dp, 86._dp, 87._dp, 88._dp, 89._dp, 90._dp, &
        !      91._dp, 92._dp, 93._dp, 94._dp, 95._dp, 96._dp, 97._dp, 98._dp, 99._dp, 100._dp ], &
        ! [10, 10]))

        ! Hilbert matrix
        ! SVD gets it somewhat correct, but still floating point errors. Also differs slightly to LU, 
        ! which also exhibits floating-point errors
        ! b= transpose(reshape(&
        ! [1._dp,          0.5_dp,         0.33333333_dp, 0.25_dp,       0.2_dp,       &
        !  0.5_dp,         0.33333333_dp,  0.25_dp,       0.2_dp,        0.16666667_dp,&
        !  0.33333333_dp,  0.25_dp,        0.2_dp,        0.16666667_dp, 0.14285714_dp,&
        !  0.25_dp,        0.2_dp,         0.16666667_dp, 0.14285714_dp, 0.125_dp,     &
        !  0.2_dp,         0.16666667_dp,  0.14285714_dp, 0.125_dp,      0.11111111_dp], [5, 5]))

contains

    subroutine test_pseudo_inverse()
        real(dp) :: a(10, 10), a_inv(10, 10), a_inv_ref(10, 10)
        real(dp) :: b(5, 5), b_inv(5, 5), b_inv_ref(5, 5)

        integer :: i

        ! Tridiagonal Matrix
        a = transpose(reshape( &
        [2.0_dp, -1.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
        -1.0_dp,  2.0_dp, -1.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
         0.0_dp, -1.0_dp,  2.0_dp, -1.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
         0.0_dp,  0.0_dp, -1.0_dp,  2.0_dp, -1.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
         0.0_dp,  0.0_dp,  0.0_dp, -1.0_dp,  2.0_dp, -1.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
         0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, -1.0_dp,  2.0_dp, -1.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
         0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, -1.0_dp,  2.0_dp, -1.0_dp,  0.0_dp,  0.0_dp, &
         0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, -1.0_dp,  2.0_dp, -1.0_dp,  0.0_dp, &
         0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, -1.0_dp,  2.0_dp, -1.0_dp, &
         0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, -1.0_dp,  2.0_dp], [10, 10]))

         a_inv_ref = transpose(reshape(&
          [0.90909091_dp, 0.81818182_dp, 0.72727273_dp, 0.63636364_dp, 0.54545455_dp, 0.45454545_dp, &
           0.36363636_dp, 0.27272727_dp, 0.18181818_dp, 0.09090909_dp, &
           0.81818182_dp, 1.63636364_dp, 1.45454545_dp, 1.27272727_dp, 1.09090909_dp, 0.90909091_dp, &
            0.72727273_dp, 0.54545455_dp, 0.36363636_dp, 0.18181818_dp, &
           0.72727273_dp, 1.45454545_dp, 2.18181818_dp, 1.90909091_dp, 1.63636364_dp, 1.36363636_dp, &
            1.09090909_dp, 0.81818182_dp, 0.54545455_dp, 0.27272727_dp, &
           0.63636364_dp, 1.27272727_dp, 1.90909091_dp, 2.54545455_dp, 2.18181818_dp, 1.81818182_dp, &
            1.45454545_dp, 1.09090909_dp, 0.72727273_dp, 0.36363636_dp, &
           0.54545455_dp, 1.09090909_dp, 1.63636364_dp, 2.18181818_dp, 2.72727273_dp, 2.27272727_dp, &
            1.81818182_dp, 1.36363636_dp, 0.90909091_dp, 0.45454545_dp, &
           0.45454545_dp, 0.90909091_dp, 1.36363636_dp, 1.81818182_dp, 2.27272727_dp, 2.72727273_dp, &
            2.18181818_dp, 1.63636364_dp, 1.09090909_dp, 0.54545455_dp, &
           0.36363636_dp, 0.72727273_dp, 1.09090909_dp, 1.45454545_dp, 1.81818182_dp, 2.18181818_dp, &
            2.54545455_dp, 1.90909091_dp, 1.27272727_dp, 0.63636364_dp, &
           0.27272727_dp, 0.54545455_dp, 0.81818182_dp, 1.09090909_dp, 1.36363636_dp, 1.63636364_dp, &
            1.90909091_dp, 2.18181818_dp, 1.45454545_dp, 0.72727273_dp, &
           0.18181818_dp, 0.36363636_dp, 0.54545455_dp, 0.72727273_dp, 0.90909091_dp, 1.09090909_dp, &
            1.27272727_dp, 1.45454545_dp, 1.63636364_dp, 0.81818182_dp, &
           0.09090909_dp, 0.18181818_dp, 0.27272727_dp, 0.36363636_dp, 0.45454545_dp, 0.54545455_dp, &
            0.63636364_dp, 0.72727273_dp, 0.81818182_dp, 0.90909091_dp], [10, 10]))

        call pseudo_inv(a, a_inv)

        call check(all_close(a_inv, a_inv_ref), msg='Inv(A) does not agree with the reference')
       
        b= transpose(reshape(&
            [-3._dp,  8._dp,  4._dp,  4._dp,  4._dp, &
             -8._dp,  3._dp,  5._dp,  0._dp,  0._dp, &
             -7._dp,  3._dp,  9._dp, -5._dp, 10._dp, &
             -4._dp, 10._dp, -1._dp,  2._dp, -7._dp, &
             -3._dp,  4._dp,  9._dp,  2._dp,  8._dp], [5, 5]))

        b_inv_ref = transpose(reshape( &
       [-0.17624148003894838_dp,      -0.18461538461538482_dp,      -1.1684518013631861E-002_dp, &
        0.11762414800389490_dp,      0.20564751703992212_dp,       &
         1.8500486854917127E-002_dp,  -0.10769230769230778_dp,       5.6475170399221099E-002_dp, &
          9.8149951314508405E-002_dp,  6.0370009737098435E-003_dp,   &
        -0.29308666017526758_dp,      -3.0769230769230820E-002_dp,  -5.2580331061343757E-002_dp, &
         0.12930866601752666_dp,      0.32541382667964935_dp,       &
         0.15936384290814667_dp,       0.11794871794871806_dp,      -0.12755598831548204_dp,     &
         -0.11593638429081465_dp,     -2.1681272314183708E-002_dp,   &
         0.21454073352807515_dp,      -1.0256410256410222E-002_dp,   5.8422590068159704E-002_dp, &
         -0.12145407335280746_dp,     -0.16157091853294381_dp], [5, 5]))

         call pseudo_inv(b, b_inv)

         call check(all_close(b_inv, b_inv_ref), msg='Inv(B) does not agree with the reference')

         

    end subroutine test_pseudo_inverse

    subroutine test_inversion_lu()
        real(dp) :: a(10, 10), a_inv(10, 10), a_inv_ref(10, 10)
        real(dp) :: b(5, 5), b_inv(5, 5), b_inv_ref(5, 5)
        integer :: i

        ! Tridiagonal Matrix
        a = transpose(reshape( &
        [2.0_dp, -1.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
        -1.0_dp,  2.0_dp, -1.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
         0.0_dp, -1.0_dp,  2.0_dp, -1.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
         0.0_dp,  0.0_dp, -1.0_dp,  2.0_dp, -1.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
         0.0_dp,  0.0_dp,  0.0_dp, -1.0_dp,  2.0_dp, -1.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
         0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, -1.0_dp,  2.0_dp, -1.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
         0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, -1.0_dp,  2.0_dp, -1.0_dp,  0.0_dp,  0.0_dp, &
         0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, -1.0_dp,  2.0_dp, -1.0_dp,  0.0_dp, &
         0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, -1.0_dp,  2.0_dp, -1.0_dp, &
         0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, -1.0_dp,  2.0_dp], [10, 10]))

         a_inv_ref = transpose(reshape(&
          [0.90909091_dp, 0.81818182_dp, 0.72727273_dp, 0.63636364_dp, 0.54545455_dp, 0.45454545_dp, &
           0.36363636_dp, 0.27272727_dp, 0.18181818_dp, 0.09090909_dp, &
           0.81818182_dp, 1.63636364_dp, 1.45454545_dp, 1.27272727_dp, 1.09090909_dp, 0.90909091_dp, &
            0.72727273_dp, 0.54545455_dp, 0.36363636_dp, 0.18181818_dp, &
           0.72727273_dp, 1.45454545_dp, 2.18181818_dp, 1.90909091_dp, 1.63636364_dp, 1.36363636_dp, &
            1.09090909_dp, 0.81818182_dp, 0.54545455_dp, 0.27272727_dp, &
           0.63636364_dp, 1.27272727_dp, 1.90909091_dp, 2.54545455_dp, 2.18181818_dp, 1.81818182_dp, &
            1.45454545_dp, 1.09090909_dp, 0.72727273_dp, 0.36363636_dp, &
           0.54545455_dp, 1.09090909_dp, 1.63636364_dp, 2.18181818_dp, 2.72727273_dp, 2.27272727_dp, &
            1.81818182_dp, 1.36363636_dp, 0.90909091_dp, 0.45454545_dp, &
           0.45454545_dp, 0.90909091_dp, 1.36363636_dp, 1.81818182_dp, 2.27272727_dp, 2.72727273_dp, &
            2.18181818_dp, 1.63636364_dp, 1.09090909_dp, 0.54545455_dp, &
           0.36363636_dp, 0.72727273_dp, 1.09090909_dp, 1.45454545_dp, 1.81818182_dp, 2.18181818_dp, &
            2.54545455_dp, 1.90909091_dp, 1.27272727_dp, 0.63636364_dp, &
           0.27272727_dp, 0.54545455_dp, 0.81818182_dp, 1.09090909_dp, 1.36363636_dp, 1.63636364_dp, &
            1.90909091_dp, 2.18181818_dp, 1.45454545_dp, 0.72727273_dp, &
           0.18181818_dp, 0.36363636_dp, 0.54545455_dp, 0.72727273_dp, 0.90909091_dp, 1.09090909_dp, &
            1.27272727_dp, 1.45454545_dp, 1.63636364_dp, 0.81818182_dp, &
           0.09090909_dp, 0.18181818_dp, 0.27272727_dp, 0.36363636_dp, 0.45454545_dp, 0.54545455_dp, &
            0.63636364_dp, 0.72727273_dp, 0.81818182_dp, 0.90909091_dp], [10, 10]))

        call inversion_lu(a, a_inv)

        call check(all_close(a_inv, a_inv_ref), msg='Inv(A) does not agree with the reference')
        
        b= transpose(reshape(&
        [-3._dp,  8._dp,  4._dp,  4._dp,  4._dp, &
         -8._dp,  3._dp,  5._dp,  0._dp,  0._dp, &
         -7._dp,  3._dp,  9._dp, -5._dp, 10._dp, &
         -4._dp, 10._dp, -1._dp,  2._dp, -7._dp, &
         -3._dp,  4._dp,  9._dp,  2._dp,  8._dp], [5, 5]))

        b_inv_ref = transpose(reshape( &
      [-0.17624148003894838_dp,      -0.18461538461538482_dp,      -1.1684518013631861E-002_dp, &
        0.11762414800389490_dp,      0.20564751703992212_dp,       &
        1.8500486854917127E-002_dp,  -0.10769230769230778_dp,       5.6475170399221099E-002_dp, &
        9.8149951314508405E-002_dp,  6.0370009737098435E-003_dp,   &
        -0.29308666017526758_dp,      -3.0769230769230820E-002_dp,  -5.2580331061343757E-002_dp, &
        0.12930866601752666_dp,      0.32541382667964935_dp,       &
        0.15936384290814667_dp,       0.11794871794871806_dp,      -0.12755598831548204_dp,     &
        -0.11593638429081465_dp,     -2.1681272314183708E-002_dp,   &
        0.21454073352807515_dp,      -1.0256410256410222E-002_dp,   5.8422590068159704E-002_dp, &
        -0.12145407335280746_dp,     -0.16157091853294381_dp], [5, 5]))
        
        call inversion_lu(b, b_inv)

        call check(all_close(b_inv, b_inv_ref), msg='Inv(B) does not agree with the reference')

    end subroutine test_inversion_lu

end program test_lapack_drivers
