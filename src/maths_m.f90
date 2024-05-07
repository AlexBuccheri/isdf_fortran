module maths_m
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    ! Exposed routines and functions
    public :: close, all_close

    interface all_close
        module procedure :: all_close_real64_1d, all_close_real64_2d
    end interface all_close

contains


    !> @brief Are \f$x\f$ and \f$y\f$ equal within a tolerance.
    elemental logical function close(x, y, rtol, atol)
        real(dp), intent(in) :: x, y
        real(dp), optional, intent(in) :: rtol
        real(dp), optional, intent(in) :: atol
        real(dp) :: atol_, rtol_

        if(present(rtol)) then
            rtol_ = rtol
        else
            rtol_ = 1.e-5_dp
        endif

        if(present(atol)) then
            atol_ = atol
        else
            atol_ = 1.e-8_dp
        endif

        close = abs(x - y) <= (atol_ + rtol_ * abs(y))
    end function close


    ! Cannot be elemental when shape of inputs differ from shape of output
    logical function all_close_real64_1d(x, y, rtol, atol)
        real(dp), intent(in) :: x(:), y(:)
        real(dp), optional, intent(in) :: rtol
        real(dp), optional, intent(in) :: atol
        real(dp) :: atol_, rtol_
        logical, allocatable :: values_close(:)

        if(present(rtol)) then
            rtol_ = rtol
        else
            rtol_ = 1.e-5_dp
        endif

        if(present(atol)) then
            atol_ = atol
        else
            atol_ = 1.e-8_dp
        endif

        ! Explicitly allocate else I can the compiler warning:
        ! `'values_close.offset' is used uninitialized [-Wuninitialized]`
        allocate(values_close(size(x)))
        values_close = abs(x - y) <= (atol_ + rtol_ * abs(y))
        all_close_real64_1d = all(values_close)

    end function all_close_real64_1d


    logical function all_close_real64_2d(x, y, rtol, atol)
        real(dp), intent(in) :: x(:, :), y(:, :)
        real(dp), optional, intent(in) :: rtol
        real(dp), optional, intent(in) :: atol
        real(dp) :: atol_, rtol_
        logical, allocatable :: values_close(:, :)

        if(present(rtol)) then
            rtol_ = rtol
        else
            rtol_ = 1.e-5_dp
        endif

        if(present(atol)) then
            atol_ = atol
        else
            atol_ = 1.e-8_dp
        endif

        ! Explicitly allocate else I can the compiler warning:
        ! `'values_close.offset' is used uninitialized [-Wuninitialized]`
        allocate(values_close(size(x, 1), size(x, 2)))
        values_close = abs(x - y) <= (atol_ + rtol_ * abs(y))
        all_close_real64_2d = all(values_close)

    end function all_close_real64_2d

end module maths_m
