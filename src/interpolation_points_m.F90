! Find the interpolation points via QR decomposition
module interpolation_points_m
    use, intrinsic :: iso_fortran_env

    use face_splitting_m, only: face_splitting_product, face_splitting_product_two_funcs_colwise
    use maths_m,          only: modified_gram_schmidt
    implicit none
    private
    public :: subsample_transpose_product_matrix

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


    !> @brief Randomly sample the transpose of the product matrix
    !!
    !! Construct random Gaussian matrix, with orthogonalised columns. Shape (p, m)
    !! Pass wave functions in packed form (m, np)
    !! Contract G and phi matrices to have shape (p, np)
    !! Construct Z^T as the face-splitting product of (G phi) (G phi), with shape
    !! (p*p, np), such that one needs to apply products to the rows
    subroutine subsample_transpose_product_matrix(phi, n_interp, zt_subspace, random_seed)
        real(real64), intent(in)    :: phi(:, :)   !< Set of KS states, with shape(m_states, np)
        integer,      intent(inout) :: n_interp    !< Number of interpolation points
        real(real64), intent(out), allocatable :: zt_subspace(:, :)
        integer,      intent(in), optional     :: random_seed(:)

        real(real64), allocatable :: G1(:, :), G_phi(:, :), z_subspace(:, :)
        integer :: np, m_states, p

        ! Dimensions
        m_states = size(phi, 1)
        np = size(phi, 2)
        p = nint(sqrt(real(n_interp, kind=real64)))

        ! Construct random Gaussian matrix, with orthogonalised columns. Shape (p, m)
        allocate(G1(p, m_states))
        call orthogonalisd_gaussian_matrix(G1, random_seed)

        ! Sanity check for shape of phi
        if (np < m_states) then
            write(*, *) 'Number of states is < n grid points. Might have phi allocated incorrectly'
            write(*, *) 'm_states, np:', m_states, np
            error stop
        endif

        ! Contract G and phi matrices to have shape (p, np): 
        ! Arrays:  G_phi =     G     phi
        ! Shapes: (p, np) = (p, m) (m, np)
        allocate(G_phi(np, p))
        call dgemm('N', 'N', &
            p,             &  ! row op(G)
            np,            &  ! col op(phi)    
            m_states,      &  ! col op(G)
            1._real64,     &       
            G1, p,         &
            phi, m_states, &
            0._real64,     &
            G_phi, np      &
        )      

        ! Z^T_subspace = Face-splitting product (G_phi) (G_phi) to give shapes:
        ! (p*p, np) = (p, np) ^ (p, np)
        n_interp = p * p
        allocate(zt_subspace(n_interp, np))
        call face_splitting_product_two_funcs_colwise(phi, zt_subspace)

    end subroutine subsample_transpose_product_matrix


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


end module interpolation_points_m
