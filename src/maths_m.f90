module maths_m
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    ! Exposed routines and functions
    public :: close, all_close, pseudo_inv, svd, gram_schmidt, modified_gram_schmidt

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


      !>@brief Sum of the projection of vector v_i onto a set of vectors {u}.
      subroutine summed_projection(v, u, v_index, proj)
        real(dp), intent(in), contiguous :: v(:)
        real(dp), intent(in), contiguous :: u(:, :)
        integer,  intent(in) :: v_index
        real(dp), intent(out) :: proj(:)

        integer :: i

        ! Assumes U in [1, v_index - 1] are normalised in the caller
        ! TODO(Alex) See if I can replace the dot product with blas
        ! TODO(Alex) Consider OMP parallelisation
        proj = dot_product(v, u(:, 1)) * u(:, 1)

        do i = 2, v_index - 1
            proj = proj + dot_product(v, u(:, i)) * u(:, i)
        enddo

      end subroutine summed_projection

      ! Fastest way to test with python
      ! def gram_schmidt_columns(X):
      !    Q, R = np.linalg.qr(X)
      !    return Q  # Q = Orthogonalised X

      !>@brief Orthogonalisation of column vectors with classic Gram-Schmidt.
      !! 
      !! TODO Add expressions
      !!
      !! Note, classic Gram Schmidt cannot handle linearly-dependent vectors.
      !! The routine will silently return zeros for corresponding vectors. 
      !! It is the caller''s responsibility to check for linearly-dependent vectors.
      subroutine gram_schmidt(v)
        real(dp), intent(inout), contiguous :: v(:, :) !< In: Array of column vectors
        !                                                Out: Orthogonalised column vectors
        integer  :: m, n_vectors, i
        real(dp) :: norm
        real(dp), allocatable :: proj(:)

        m = size(v, 1)
        n_vectors = size(v, 2)
        allocate(proj(m))

        ! TODO(Alex) See if I can replace norm call with blas
        norm = norm2(v(:, 1))
        call dscal(m, 1._dp / norm, v(:, 1), 1)

        do i = 2, n_vectors
            call summed_projection(v(:, i), v, i, proj)
            ! v(:, i) = v(:, i) - proj(:)
            call daxpy(m, -1._dp, proj, 1, v(:, i), 1)
            norm = norm2(v(:, i))
            ! Handle linearly-dependent vectors
            if (close(norm, 0._dp)) then
                v(:, i) = 0._dp
            endif
            call dscal(m, 1._dp / norm, v(:, i), 1)
        enddo

      end subroutine gram_schmidt


      !> @brief Orthogonalisation of column vectors with modified Gram-Schmidt.
      !!
      !! The projection of each vector v(:,j) is subtracted from only the remaining vectors 
      !! (v(:,j+1) to v(:,n)). This orthogonalizes the remaining vectors with respect to the 
      !! newly orthogonalized v(:,j), and results in a more numerically stable algorithm than the 
      !! classical version, especially for ill-conditioned matrices.
      !! 
      !! Also see https://laurenthoeltgen.name/post/gram-schmidt/
      subroutine modified_gram_schmidt(v)
        real(dp), intent(inout), contiguous :: v(:, :) !< In: Array of column vectors
        !                                                Out: Orthogonalised column vectors
        integer  :: n_vectors, m, i, j
        real(dp) :: norm

        m = size(v, 1)
        n_vectors = size(v, 2)
    
        do j = 1, n_vectors
            norm = norm2(v(:, j))
            ! Handle linearly-dependent vectors
            if (close(norm, 0._dp)) then
                v(:, j) = 0._dp
                cycle
            endif
            call dscal(m, 1._dp / norm, v(:, j), 1)
            do i = j + 1, n_vectors
                v(:, i) = v(:, i) - dot_product(v(:, j), v(:, i)) * v(:, j)
            enddo
        enddo

      end subroutine modified_gram_schmidt


end module maths_m
