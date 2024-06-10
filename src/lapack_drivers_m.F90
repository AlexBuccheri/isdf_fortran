module lapack_drivers_m
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private 
    public :: qr_decomposition_with_pivot, inversion_lu, pseudo_inv, svd

contains

   !> @bried Inversion of a matrix using LU decomposition
   !! 
   !! TODO Add checks for non-square matrices
   !! Make a_inv optional, returnin the inverse in a
   subroutine inversion_lu(a, a_inv)
      real(dp)              :: a(:, :)      !< Matrix to inverse
      real(dp), intent(out) :: a_inv(:, :)  !< Inverse(A), expected to be allocated memory by the caller

      real(dp), allocatable :: tmp_a(:, :), work(:)
      integer,  allocatable :: ipiv(:)
      integer               :: m, n, info, lwork

      m = size(a, 1)
      n = size(a, 2)
      allocate(tmp_a(m, n), source=a)
      allocate(ipiv(min(m, n)))

      ! LU factorization
      call dgetrf(m, n, A, m, ipiv, info)

      if (info /= 0) then
         write(*,*) 'Error with getrf LU factorisation: ', info
         stop
      endif
      
      ! On exit, the factors L and U from the factorization
      ! A = P*L*U; the unit diagonal elements of L are not stored 
      ! Pass outputs (A and ipiv) directly to `getri`

      ! Query workspace
      allocate(work(1))
      lwork = -1
      call dgetri(n, A, m, ipiv, work, lwork, info)
      lwork = int(work(1))
      deallocate(work)

      ! Compute inverse using LU decomposition
      ! TODO See if this has an effect
      !allocate(work(lwork))
      allocate(work(n))
      call dgetri(n, A, m, ipiv, work, size(work), info)
      
      if (info /= 0) then
         write(*,*) 'Error with getri, inverting matrix A: ', info
         stop
      endif

      deallocate(ipiv)

      ! Return data in correct arrays (Should add a warning if not a square matrix)
      a_inv = a
      a = tmp_a
      deallocate(tmp_a)

   end subroutine inversion_lu
   

   !> @brief Compute pseudo-inv with SVD
   ! TODO(Alex) See if A is preserved.
   subroutine pseudo_inv(A, A_inv)
      real(dp), intent(in)   :: A(:, :)
      real(dp), intent(out)  :: A_inv(:, :)

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
      allocate (s_inv(n, m), source=0._dp)
      do i = 1, size(s, 1)
         s_inv(i, i) = 1._dp/s(i)
      end do
      deallocate (S)

      ! Contract (V S^+) = VS
      allocate (VSi(n, n))
      call dgemm('T', 'N', &
                 size(VT, 2), &  ! Rows of op(A). op(A) = V
                 size(s_inv, 2), &  ! Cols of op(B). op(B) = s_inv
                 size(VT, 1), &  ! Cols of op(A)
                 1._dp, &
                 VT, &
                 size(VT, 1), &  ! Rows of A
                 s_inv, &
                 size(s_inv, 1), &  ! Rows of B
                 0._dp, &
                 VSi, n)            ! Rows of C
      deallocate (VT)

      ! Contract A_inv = (V S^+) U^T
      call dgemm('N', 'T', &
                 size(VSi, 1), &    ! Rows of op(A). op(A) = VSi
                 size(U, 1), &    ! Cols of op(B). op(B) = U^T
                 size(VSi, 2), &    ! Cols of op(A)
                 1._dp, &
                 VSi, &
                 size(VSi, 1), &    ! Rows of A
                 U, &
                 size(U, 1), &    ! Rows of B
                 0._dp, &
                 A_inv, m)          ! Rows of C

      deallocate (VSi)
      deallocate (U)

   end subroutine pseudo_inv


   !> @brief SVD wrapper
   !!
   !! Workspace query size taken from:
   !! https://github.com/numericalalgorithmsgroup/LAPACK_Examples/blob/master/examples/source/dgesvd_example.f90
   subroutine svd(A, U, S, VT)
      real(dp), intent(in)  :: A(:, :)
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

      allocate (s(min(m, n)))
      allocate (u(ldu, m))
      allocate (vt(ldvt, n))

      ! Query optimal work space
      lwork = -1
      call dgesvd('A', 'S', m, n, A, lda, s, u, ldu, vt, ldvt, dummy, lwork, info)

      ! Compute the singular values and left and right singular vectors of A.
      lwork = max(m + 4*n + nb*(m + n), nint(dummy(1, 1)))
      allocate (work(lwork))
      call dgesvd('A', 'S', m, n, A, lda, s, u, ldu, vt, ldvt, work, lwork, info)

      if (info /= 0) then
         write (*, *) 'Failure in DGESVD. INFO =', info
      end if

      deallocate (work)

   end subroutine svd


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
      write(*, *) 'Workspace queried'
      call dgeqp3(n, m, A, n, p, tau, work, lwork, info)

      ! Allocate optimal work space
      lwork = int(work(1))
      deallocate (work)
      allocate (work(lwork))

      ! A P = Q R
      write(*, *) 'Perform QR with pivoting'
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
