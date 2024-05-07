program test_isdf_m
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib

    use fortuno_interface_m, only: execute_cmd_app, test, check, is_equal
    use maths_m,             only: all_close

    implicit none

    ! Register tests
    call execute_cmd_app(testitems=[&
            test("Placeholder", test_stub) &
        ])

contains

    subroutine test_stub()
    end subroutine test_stub

end program test_isdf_m
