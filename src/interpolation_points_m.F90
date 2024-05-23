! Find the interpolation points via QR decomposition
module interpolation_points_m
    use, intrinsic :: iso_fortran_env

    use face_splitting_m, only: face_splitting_product
    use maths_m,          only: modified_gram_schmidt
    implicit none
    private

contains

    !> @brief Generate an orthogonalised 2D matrix of random numbers, sampled from a 
    !! Gaussian distribution
    subroutine orthogonalisd_gaussian_matrix(gaussian, seed)
        real(real64), intent(out)    :: gaussian(:, :)
        integer, intent(in), optional :: seed(:)

        integer :: i, j, n, m
        real(real64) :: random_num

        n = size(gaussian, 1)
        m = size(gaussian, 1)

        if (present(seed)) then
            call random_seed(put=seed)
        else
            call random_seed()
        endif

        ! I have not used reservoir sampling as I am not worried about
        ! getting the same random number > once
        do j = 1, m
            do i = 1, n
                ! Generate random number between 0 and 1                                                                       
                call random_number(random_num)
                gaussian(i, j) = random_num
            enddo
        enddo

        ! Orthogonalise columns of gaussian
        call modified_gram_schmidt(gaussian)

    end subroutine orthogonalisd_gaussian_matrix

    !> @bried Randomly sample the product matrix using Gaussian test matrices.
    !!
    !!\f[
    !!    \tilde{Z}_{\alpha, \beta} = \left( \sum^m_{i=1} \phi_i(\mathbf{r}) G^{\phi}_{i, \alpha} \right)
    !!                                \left( \sum^m_{i=1} \phi_i(\mathbf{r}) G^{\phi}_{i, \beta} \right)
    !! \f]
    !!
    !! Implemented according to eq. 20 of "Interpolative Separable Density Fitting Decomposition for
    !! Accelerating Hybrid Density Functional Calculations with Applications to Defects in Silicon"
    !! J. Chem. Theory Comput. 2017, 13, 5420-5431
    subroutine randomly_sample_product_matrix(phi, n_interp, z_subspace, random_seed)
        real(real64), intent(in)    :: phi(:, :)
        integer,      intent(inout) :: n_interp    !< Number of interpolation points
        real(real64), intent(out), allocatable :: z_subspace(:, :)
        integer,      intent(in), optional     :: random_seed(:)

        real(real64), allocatable :: G1(:, :), A(:, :)
        integer :: np, m_states, p

        np = size(phi, 1)
        m_states = size(phi, 2)
        p = int(sqrt(real(n_interp, kind=real64)))
        n_interp = p * p

        ! The first index of G1 SHOULD be the state index, however
        ! as the Gram-Schmidt implementation returns orthogonalised columns, set the columns to
        ! be the state index, and use the transpose of G1 in the dgemm call
        allocate(G1(p, m_states))
        call orthogonalisd_gaussian_matrix(G1, random_seed)

        allocate(A(np, p))
        call dgemm('N', 'T', &
            np,           &       ! Row op(phi) i.e. row(phi)
            p,            &       ! Col op(A) i.e. col(A)
            m_states,     &       ! Col op(phi) == row op(G1). Use col(phi)
            1._real64,    &       
            phi, np,      &       ! Rows of phi, as declared
            G1, p,        &       ! Rows of G1, as declared
            0._real64,    &
            A, np         &       ! Rows of A
        )        
        deallocate(G1)

        allocate(z_subspace(np, n_interp))
        call face_splitting_product(A, z_subspace)
        deallocate(A)

    end subroutine randomly_sample_product_matrix

    ! https://www.netlib.org/lapack/lug/node42.html
    ! https://netlib.org/lapack/explore-html/d0/dea/group__geqp3_gae96659507d789266660af28c617ef87f.html#gae96659507d789266660af28c617ef87f
    ! https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.qr.html

    !> @brief Find interpolation points from a subspace of the pair-product matrix,
    !!  using QR decomposition with pivoting.
    !!
    !! pivot determine the rows (check.. columns?) the product-basis matrix, Z, to be used as interpolating points. 
    !! Note, n_rows should be the same for Z and the z_subspace.
    subroutine interpolation_points_via_qrpivot(z_subspace, indices)
        real(real64), intent(inout)  :: z_subspace(:, :)  !<  Sub-sampled pair-product matrix, with shape(n_grid_points, n_interp)
        integer,      intent(out) :: indices(:)           !<  Interpolation point indices. 
        
        integer              :: np, n_interp, info, lwork
        integer, allocatable :: jpvt(:)
        real(real64), allocatable :: tau(:), work(:)

        np = size(z_subspace, 1)
        n_interp = size(z_subspace, 2)

        ! Query the workspace
        allocate(tau(min(np, n_interp)))
        allocate(work(1))
        allocate(jpvt(n_interp), source=0)
        lwork = -1
        call dgeqp3(np, n_interp, z_subspace, np, jpvt, tau, work, lwork, info)

        ! Allocate optimal work space
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        ! A P = Q R
        ! Upper triangle of A will contain R
        ! Below the diagonal, together with tau, represent the orthogonal matrix Q as a 
        ! product of min(M,N) elementary reflectors.
        call dgeqp3(np, n_interp, z_subspace, np, jpvt, tau, work, lwork, info)
        indices = jpvt(1: n_interp)

        ! TODO(Alex) Try and rebuild A from Q R to confirm I understand what the quantities are doing
        ! Extract R from upper triangle of A
        ! Extract lower triangle of A into Q.
        ! Reconstruct Q with something like DORGQR
        ! A_reconstructed = QR
        ! Note, this does not account for the pivoting
    

    end subroutine interpolation_points_via_qrpivot


end module interpolation_points_m
