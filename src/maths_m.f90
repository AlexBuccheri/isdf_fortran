module maths_m
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    ! Exposed routines and functions
    public :: close, all_close, pseudo_inv, svd

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


    !> @brief Compute pseudo-inv with SVD
    ! TODO(Alex) See if A is preserved.
    subroutine pseudo_inv(A, A_inv)
        real(dp),              intent(in)   :: A(:, :)
        real(dp),              intent(out)  :: A_inv(:, :)
    
        real(dp), allocatable :: u(:, :), S(:), vt(:, :) !< Result of SVD on A
        real(dp), allocatable :: s_inv(:, :)             !< inverse of the diagonals of S
        real(dp), allocatable :: VSi(:, :)               !< Contraction of V and inverse(S)
        integer :: i, m, n
        
        call svd(A, U, S, VT)
    
        ! Compute inverse of A: A^+ = V S^+ U^T
        m = size(A, 1)
        n = size(A, 2)

        ! Invert the diagonals of S, to give S^+
        ! TODO. Look into refactoring to see if I can use a vector instead of a diagonal matrix
        allocate(s_inv(m, n), source=0._dp)
        do i = 1, size(s, 1)
           s_inv(i, i) = 1._dp / s(i)
        enddo
        deallocate(S)
    
        ! NOTE, I "think" I have the leading dimensions correct when requesting the OP on a transpose
        ! As the test case is square, it may be masking an error
        ! Contract (V S^+) = VS
        allocate(VSi(n, n))
        call dgemm('T', 'N', size(VT, 1), size(s_inv, 2), size(VT, 2), 1._dp, VT, size(VT, 2), s_inv, size(s_inv, 1), 0._dp, VSi, n)
        deallocate(VT)
    
        ! Contract A_inv = (V S^+) U^T
        call dgemm('N', 'T', size(VSi, 1), size(U, 2), size(VSi, 2), 1._dp, VSi, size(VSi, 1), U, size(U, 2), 0._dp, A_inv, m)
    
        deallocate(VSi)
        deallocate(U)
            
      end subroutine pseudo_inv
    
    
      !> @brief SVD wrapper
      !!
      !! Workspace query size taken from:
      !! https://github.com/numericalalgorithmsgroup/LAPACK_Examples/blob/master/examples/source/dgesvd_example.f90
      subroutine svd(A, U, S, VT)
        real(dp),              intent(in)  :: A(:, :)
        real(dp), allocatable, intent(out) :: u(:, :)  !< U matrix
        real(dp), allocatable, intent(out) :: S(:)     !<  min(m,n) singular values of A
        real(dp), allocatable, intent(out) :: vt(:, :) !< V^T
    
        real(dp), allocatable :: work(:)
        real(dp) :: dummy(1, 1)
        integer, parameter :: nb = 64
        integer :: m, n, lda, ldu, ldvt, info, lwork
    
        !                     A    =     U       S       VT
        ! with shape:       (m,n)  =  (m, m)  (m, n)   (n, n) 
        ! which reduces to: (m,n)  =  (m, m) min(m, n) (n, n)
        ! if S is represented as a vector.
        m = size(A, 1)
        n = size(A, 2)
        lda = m
        ldu = m
        ldvt = n
        
        allocate(s(min(m, n)))
        allocate(u(ldu, m))
        allocate(vt(ldvt, n))
    
        ! Query optimal work space
        lwork = -1
        call dgesvd('A', 'S', m, n, A, lda, s, u, ldu, vt, ldvt, dummy, lwork, info)
    
        ! Compute the singular values and left and right singular vectors of A.
        lwork = max(m+4*n+nb*(m+n), nint(dummy(1,1)))
        allocate(work(lwork))
        call dgesvd('A', 'S', m, n, A, lda, s, u, ldu, vt, ldvt, work, lwork, info)
        
        if (info/=0) then
           write (*, *) 'Failure in DGESVD. INFO =', info
        endif
    
        deallocate(work)
    
      end subroutine svd

end module maths_m
