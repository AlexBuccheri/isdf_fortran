module lapack_drivers_m
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private 
    public :: qr_decomposition_with_pivot

contains

   !>@brief Perform QR decomposition with column-pivoting, on matrix A
   !!
   !! Utilises lapack `dgeqp3` to obtain R and p, and `dorgqr` to obtain Q.
   !!
   !! The matrix Q is represented as a product of elementary reflectors
   !!  Q = H(1) H(2) . . . H(k), where k = min(m,n).
   !! Each H(i) has the form:
   !!  H(i) = I - tau * v * v**T
   !! where tau is a real scalar, and v is a real/complex vector
   !! with v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in
   !! A(i+1:m,i), and tau in TAU(i).
   subroutine qr_decomposition_with_pivot(A, p, Q, R, preserve_A)
      real(dp), intent(inout)         :: A(:, :)     !< Input matrix to decompose as Q and R
      integer,  intent(out)           :: p(:)        !< Maps new column order to old column order, size(m)
      real(dp), intent(out), optional :: Q(:, :)     !< Orthogonal matrix
      real(dp), intent(out), optional :: R(:, :)     !< Upper triangular matrix
      logical,  intent(in),  optional :: preserve_A  !< Preserve A, true by default. Else A is mutated
      !                                                 and contains R in the upper triangle
      logical               :: keep_A
      integer               :: n, m, info, lwork, i, j, cnt
      real(dp), allocatable :: tau(:)                !< The scalar factors of the elementary reflectors.
      real(dp), allocatable :: work(:), tmp_A(:, :)  !< Work arrays

      n = size(A, 1)
      m = size(A, 2)

      if (present(preserve_A)) then
         keep_A = preserve_A
      else
         keep_A = .true.
      end if

      if (keep_A) allocate (tmp_A(n, m), source=A)

      allocate (tau(min(n, m)))

      ! Initialise pivot indices to zero
      p = 0

      ! TODO. I need to do this on the transpose of A, such that the permutation
      ! matrix is shape(np), not ninterp.
      ! Query the workspace
      allocate (work(1))
      lwork = -1
      call dgeqp3(n, m, A, n, p, tau, work, lwork, info)

      ! Allocate optimal work space
      lwork = int(work(1))
      deallocate (work)
      allocate (work(lwork))

      ! A P = Q R
      call dgeqp3(n, m, A, n, p, tau, work, lwork, info)
      deallocate (work)

      if (info /= 0) then
         write (*, *) 'Warning: dgeqp3 returned info: ', info
         error stop 101
      end if

      if (present(R)) then
         ! Extract R from upper triangle of A
         ! Note, memory access is suboptimal
         R = 0._dp
         do j = 1, m
            do i = 1, j
               R(i, j) = A(i, j)
            end do
         end do
      end if

      if (present(Q)) then
         ! Query optimal workspace
         allocate (work(1))
         lwork = -1
         ! min(n, m) == size(tau), the number of elementary reflectors
         call dorgqr(n, m, min(n, m), A, n, tau, work, lwork, info)

         ! Allocate optimal work space
         lwork = int(work(1))
         deallocate (work)
         allocate (work(lwork))

         ! Compute Q from the lower triangle of A, and tau
         call dorgqr(n, m, min(n, m), A, n, tau, work, lwork, info)

         if (info /= 0) then
            write (*, *) 'Warning: dgeqp3returned info: ', info
            error stop 102
         end if

         Q = A
      end if

      if (keep_A) A = tmp_A

   end subroutine qr_decomposition_with_pivot

end module lapack_drivers_m
