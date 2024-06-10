!> ISDF non-MPI implementation
module isdf_serial_m
    use, intrinsic :: iso_fortran_env, only: dp => real64
    
    use face_splitting_m, only: face_splitting_product
    use lapack_drivers_m, only: inversion_lu, pseudo_inv

    implicit none
    private

    public :: construct_interpolation_vectors, construct_approximation_product_states, &
        construct_quasi_density_matrix_R_intr, &
        construct_quasi_density_matrix_R_intr_alt, &
        select_interpolation_points_of_first_dim, &
        construct_approx_product_states_summing, &
        construct_zct, &
        construct_inverse_coefficients_contraction


contains

    !> @brief Construct the quasi-density matrix, with index one defined for all grid points
    !! and index two only defined for interpolation points.
    !!
    !! The quasi density matrix for a set of \f$m\f$ states, \f$\{\varphi\}\f$, is defined as:
    !! \f[
    !!    P^{\varphi}(\mathbf{r}, \mathbf{r}_\mu) = \sum_{i=1}^{m}  \varphi_i(\mathbf{r}) \varphi_i(\mathbf{r}_\mu),
    !! \f]
    !! where \f\mathbf{r}\f$ is defined over all grid points, and \f$\mathbf{r}_\mu\f$ is defined
    !! over interpolation points of the grid.
    subroutine construct_quasi_density_matrix_R_intr(phi, indices, P_phi)
        real(dp), intent(in)               :: phi(:, :)       !< A set of states defined on real-space grid
        !                                                        of shape (np, nstates)
        integer,  intent(in)               :: indices(:)      !< Indices of interpolation grid points
        real(dp), intent(out), allocatable :: P_phi(:, :)     !< Quasi-density matrix

        integer :: np        !< Number of grid points
        integer :: nintr     !< Number of interpolation vectors
        integer :: m_states  !< Number of KS states in the set \f$\{\varphi\}\f$

        integer  :: i, ir, jr, jintr
        real(dp) :: phi_jr

        np = size(phi, 1)
        m_states = size(phi, 2)
        nintr = size(indices)
        allocate(P_phi(np, nintr))

        do i = 1, m_states
            do jintr = 1, nintr
                jr = indices(jintr)
                phi_jr = phi(jr, i)
                do ir = 1, np
                    P_phi(ir, jintr) = P_phi(ir, jintr) + phi(ir, i) * phi_jr
                enddo
            enddo
        enddo

    end subroutine construct_quasi_density_matrix_R_intr


    !> This routine instead constructs P_imu,i at the expense of storing ninterpolation * m_state elements
    !! then performs a contraction over the states to define P^{\varphi}(\mathbf{r}, \mathbf{r}_\mu)
    subroutine construct_quasi_density_matrix_R_intr_alt(phi, indices, P_phi)
        real(dp), intent(in)               :: phi(:, :)       !< A set of states defined on real-space grid
        !                                                        of shape (np, m_states)
        integer,  intent(in)               :: indices(:)      !< Indices of interpolation grid points
        real(dp), intent(out), allocatable :: P_phi(:, :)     !< Quasi-density matrix

        real(dp), allocatable :: phi_intr(:, :)  !< KS states, only defined at interpolation points
        integer :: np                             !< Number of grid points
        integer :: nintr                          !< Number of interpolation vectors
        integer :: m_states                       !< Number of KS states in the set \f$\{\varphi\}\f$

        integer :: i, iintr, ir

        np = size(phi, 1)
        m_states = size(phi, 2)
        nintr = size(indices)

        ! Construct the density matrix, defined only for interpolation points
        allocate(phi_intr(nintr, m_states))
        do i = 1, m_states
            do iintr = 1, nintr
                ir = indices(iintr)
                phi_intr(iintr, i) = phi(ir, i)
            enddo
        enddo

        ! Contract over the state index, P = phi @ phi_intr^T
        allocate(P_phi(np, nintr))
        call dgemm ('N', 'T',   &
            size(phi, 1),       & ! size(op(A), 1)
            nintr,              & ! size(C, 2)
            size(phi, 2),       & ! size(op(A), 2)
            1._dp,              & ! alpha
            phi, np,            & ! A and LDA
            phi_intr, nintr,    & ! B and LDB
            0._dp,              & ! beta
            P_phi, np           & ! C and LDC
            )
        deallocate(phi_intr)

    end subroutine construct_quasi_density_matrix_R_intr_alt


    !> @brief Construct the product of quasi-density matrices, with index one defined for all grid points
    !! and index two only defined for interpolation points.
    ! This should just be ZC^T = (P_phi * P_psi)


    ! @brief Convert P matrix from P(r, r_u) to P(r_v, r_u)
    subroutine select_interpolation_points_of_first_dim(P_phi_ri, indices, P_phi_ii)
        real(dp),  intent(in)  :: P_phi_ri(:, :)
        integer,   intent(in)  :: indices(:)
        real(dp),  intent(out) :: P_phi_ii(:, :)

        integer :: i, j, ir, nintr

        nintr = size(indices)

        if (size(P_phi_ri, 2) /= nintr) then
            write(*, *) 'Number Interpolation point indices /= size(P_phi_ri, 2)'
            error stop 101
        endif

        if (all(shape(P_phi_ii) /= [nintr, nintr])) then
            write(*, *) 'Size of density matrix should be (nintr, nintr), but is', shape(P_phi_ii)
            error stop 102
        endif

        !$omp parallel do simd collapse(2) default(shared) private(ir)
        do j = 1, nintr
            do i = 1, nintr
                ir = indices(i)
                P_phi_ii(i, j) = P_phi_ri(ir, j)
            enddo
        enddo
        !$omp end parallel do simd

    end subroutine select_interpolation_points_of_first_dim


    !> @brief Compute the inverse of the contraction of coefficients \f$\mathbf{C}\f$ and \f$\mathbf{C}^T\f$.
    !! 
    !! \f[
    !!    \mathbf{C} \mathbf{C}^T = \mathbf{P}^{\varphi} \odot \mathbf{P}^{\psi},
    !! \f]    
    !! where both \f$\mathbf{P}\f$ matrices are square matrices, defined over interpolation points,
    !! \f$(\mathbf{r}_\nu, \mathbf{r}_\mu)\f$.
    subroutine construct_inverse_coefficients_contraction(P_phi, P_psi, cct_inv)
        real(dp), intent(in)  :: P_phi(:, :)
        real(dp), intent(in)  :: P_psi(:, :)
        real(dp), intent(out) :: cct_inv(:, :)

        integer  :: nintr, i, j
        real(dp), allocatable :: cct(:, :)

        if (all(shape(P_phi) /= shape(P_psi))) then
            write(*, *) 'Both quasi-density matrices should be the same shape: (nintr, nintr)'
            error stop 101
        endif

        if (all(shape(P_phi) /= shape(cct_inv))) then
            write(*, *) 'CCT matrix should be the same shape as the quasi-density matrices: (nintr, nintr)'
            error stop 102
        endif

        ! Element wise product [CC^T] = P_phi \odot P_psi
        nintr = size(P_phi, 1)
        allocate(cct(nintr, nintr))

       !$omp parallel do simd collapse(2) default(shared)
        do j = 1, nintr
            do i = 1, nintr
                cct(i, j) = P_phi(i, j) * P_psi(i, j)
                ! TEMP: Return CC^T instead of the inverse of CC^T
                cct_inv(i, j) = cct(i, j)
            enddo
        enddo
        !$omp end parallel do simd

        !write(*, *) 'Inverting CC^T using LU decomposition'
        !call inversion_lu(cct, cct_inv)
        ! write(*, *) 'Inverting CC^T using SVD'
        ! call pseudo_inv(cct, cct_inv)

        ! deallocate(cct)

    end subroutine construct_inverse_coefficients_contraction


    !> @brief Construct the product of Z and coefficient matrices
    !! from the element-wise product of density matrices.
    subroutine construct_zct(P_phi_r_mu, ZCT, P_psi_r_mu)
        real(dp), intent(in), target :: P_phi_r_mu(:, :)  !< Density matrix for states {\phi}
        real(dp), intent(out) :: ZCT(:, :)
        real(dp), intent(in), optional, target :: P_psi_r_mu(:, :)  !< Density matrix for states {\psi}

        real(dp), pointer :: P_psi(:, :)

        integer :: i_mu, ir, np, nintr

        if (all(shape(P_phi_r_mu) /= shape(ZCT))) then
            write(*, *) 'Shape of P^phi_{r,mu} /= shape(ZC^T)'
            error stop
        endif

        if (present(P_psi_r_mu)) then
            if (all(shape(P_psi_r_mu) /= shape(ZCT))) then
                write(*, *) 'Shape of P^psii_{r,mu} /= shape(ZC^T)'
                error stop
            endif
            P_psi => P_psi_r_mu
        else
            P_psi => P_phi_r_mu
        endif

        np = size(ZCT, 1)
        nintr = size(ZCT, 2)

        !$omp parallel do simd collapse(2) default(shared)
        do i_mu = 1, nintr
            do ir = 1, np
                ZCT(ir, i_mu) = P_phi_r_mu(ir, i_mu) * P_psi(ir, i_mu)
            end do
        enddo
        !$omp end parallel do simd

        nullify(P_psi)

    end subroutine construct_zct


    !> @brief Construct a matrix of ISDF interpolation vectors, \f$\mathbf{\Theta}\f$.
    !! 
    !! This API assumes that one only has one set of KS states
    !! Address assigning optional second set of states to a pointer once first implementation tested.
    !! (see commented code below)
    !!
    !! The implementation is as follows:
    !!
    !! 1. Construct the quasi-density matrix \f$ P^{\varphi}(\mathbf{r}, \mathbf{r}_\mu) \f$
    !!    where the first variable is defined over the whole grid and the second variable runs
    !!    only over interpolation points.
    !! 2. Construct the inverse of the contracted coefficients matrix:
    !!    a) Restrict the first variable of \f$ P^{\varphi}(\mathbf{r}, \mathbf{r}_\mu) \f$ to
    !!       interpolation points
    !!    b) Perform the inversion on \f$P^{\varphi}(\mathbf{r}_\nu, \mathbf{r}_\mu) \odot P^{\psi}(\mathbf{r}_\nu, \mathbf{r}_\mu)\f$
    !! 3. Construct the interpolation vectors, \f$\theta \f$
    !!    a) Construct \f$(\mathbf{Z}\mathbf{C}^T)_{(k, \mu)} = P^{\varphi}(\mathbf{r}_k, \mathbf{r}_\mu) \odot P^{\psi}(\mathbf{r}_k, \mathbf{r}_\mu)\f$ 
    !!      where \f$k\f$ runs over all grid points and \f$\mu\f$ runs over all interpolation points.
    !!    b) Contract \f$\mathbf{Z}\mathbf{C}^T\f$ with  \f$(\mathbf{C}\mathbf{C}^T)^{-1}\f$
    !!
    !! For additional details, consult Fig 3. in: [Interpolative Separable Density Fitting Decomposition for Accelerating Hybrid 
    !! Density Functional Calculations with Applications to Defects in Silicon](10.1021/acs.jctc.7b00807)
    !! 
    subroutine construct_interpolation_vectors(phi, indices, theta)
        real(dp),              intent(in)  :: phi(:, :)      !< A set of states defined on real-space grid
        !                                                       of shape (np, mstates)
        integer,               intent(in)  :: indices(:)     !< Indices of interpolation grid points (nintr)
        real(dp), allocatable, intent(out) :: theta(:, :)    !< Matrix of interpolating vectors, stored columnwise (np, nintr)

        ! Work arrays, where r = realspace point and i=interpolation point (k amd mu/nu in the paper)
        integer                :: np, nintr, ir, iintr
        real(dp),  allocatable :: P_phi_ri(:, :), P_phi_ii(:, :)
        real(dp),  allocatable :: ZCT(:, :), cct(:, :), cct_inv(:, :)

        ! Allocate interpolation vectors matrix
        np = size(phi, 1)
        nintr = size(indices)

        ! 1. Quasi density matrix
        write(*, *) 'Constructing density matrix'
        call construct_quasi_density_matrix_R_intr(phi, indices, P_phi_ri)
 
        ! 2. Construct CC^T
        ! a) Convert P(r, r_u) -> P(r_v, r_u)
        ! TODO(Alex) Should change some of this notation from i(nterpolation) to mu, nu
        write(*, *) 'Convert P to interpolation points only'
        allocate(P_phi_ii(nintr, nintr))
        call select_interpolation_points_of_first_dim(P_phi_ri, indices, P_phi_ii)
        
        ! b) Contract CC^T (element-wise multiply P_phi_ii and P_phi_ii), and invert
        write(*, *) 'Construct and invert CC^T'
        allocate(cct(nintr, nintr), cct_inv(nintr, nintr))
        call construct_inverse_coefficients_contraction(P_phi_ii, P_phi_ii, cct)
        ! TODO Alex. Note that inversion is commented out in `construct_inverse_coefficients_contraction`
        ! so returning cct here
        write(*, *) 'Inverting CC^T using LU decomposition'
        call inversion_lu(cct, cct_inv)
        deallocate(cct)
        deallocate(P_phi_ii)

        ! 3. Construct theta = (ZC^T)(CC^T)^{-1} == [P_phi_ri * P_psi_ri] @ [P_phi_ii * P_psi_ii]^{-1}
        ! a) (ZC^T) == [P_phi_ri * P_psi_ri] i.e. Element-wise multiplication
        write(*, *) 'Construct ZC^T'
        allocate(ZCT(np, nintr))
        call construct_zct(P_phi_ri, ZCT)
        deallocate(P_phi_ri)

        ! b) Contract (ZC^T) with (CC^T)^{-1}
        write(*, *) 'Contract (ZC^T) (CC^T)^{-1}'
        allocate(theta(np, nintr))
        call dgemm ('N', 'N', &
            np, &    ! row(ZCT)
            nintr, & ! col(cct_inv)
            nintr, & ! col(ZCT)
            1._dp, & ! alpha
            ZCT, np, & !LDA
            cct_inv, nintr, &  ! LDB
            0._dp, &  ! beta
            theta, np &  ! LDC
            )
        deallocate(ZCT)
        deallocate(cct_inv)

    end subroutine construct_interpolation_vectors


    ! TODO(Alex) Could test against a loop implementation
    !> @brief Approximate the product basis using ISDF vectors.
    !!
    !! \f[
    !!     \varphi_i(\mathbf{r}) \psi_j(\mathbf{r}) 
    !!       \approx \sum_{\mu=1}^{N_\mu} \zeta_\mu(\mathbf{r})\varphi_i(\mathbf{r}_\mu) \psi_j(\mathbf{r}_\mu)    
    !! \f]
    subroutine construct_approximation_product_states(theta, phi, psi, indices, z_approx)
        real(dp), intent(in) :: theta(:, :)                  !< ISDF interpolation vectors of shape(Np, Nintr)
        real(dp), intent(in) :: phi(:, :)                    !< KS states of shape (Np, m_states)
        real(dp), intent(in) :: psi(:, :)                    !< KS states of shape (Np, n_states)
        integer, intent(in)  :: indices(:)                   !< Indices of interpolation points
        real(dp), allocatable, intent(out) :: z_approx(:, :) !< Approximate product basis 

        integer :: np, m, n, nintr, i, ir, imu
        real(dp), allocatable :: phi_intr(:, :), psi_intr(:, :), z_intr(:, :)

        np = size(phi, 1)
        m = size(phi, 2)
        n = size(psi, 2)
        nintr = size(theta, 2)

        if (size(psi, 1) /= np) then
            write(*, *) 'First dimension of psi inconsistent with phi'
            error stop
        endif

        if (size(theta, 1) /= np) then
            write(*, *) 'First dimension of theta inconsistent with phi'
            error stop
        endif

        ! 1. Construct z_{ij}(r_\mu) by face-splitting \varphi_i(\mathbf{r}_\mu) and \psi_j(\mathbf{r}_\mu)    
        ! a) Construct phi and psi at interpolation points (Nmu, n or m)
        allocate(phi_intr(nintr, m))
!        !$omp parallel do simd collapse(2) default(shared) private(ir)
        do i = 1, m
            do imu = 1, nintr
                ir = indices(imu)
                phi_intr(imu, i) = phi(ir, i)
            enddo
        enddo
 !       !$omp end parallel do simd

        allocate(psi_intr(nintr, n))
!        !$omp parallel do simd collapse(2) default(shared) private(ir)
        do i = 1, n
            do imu = 1, nintr
                ir = indices(imu)
                psi_intr(imu, i) = psi(ir, i)
            enddo
        enddo
 !       !$omp end parallel do simd

        ! b) Face split
        allocate(z_intr(nintr, n * m))
        call face_splitting_product(phi_intr, psi_intr, z_intr)
        deallocate(phi_intr)
        deallocate(psi_intr)

        ! 2. Contract Theta and Z_intr of shapes (np, Nmu) (Nmu, n * m) -> (np, n * m)
        allocate(z_approx(np, n * m))
        call dgemm ('N', 'N', &
            np, &    !row(theta)
            n * m, & ! col(z_intr)
            nintr, & ! col(theta)
            1._dp, &
            theta, np, &
            z_intr, nintr, &
            0._dp, &
            z_approx, np &
            )
        deallocate(z_intr)

    end subroutine construct_approximation_product_states


    subroutine construct_approx_product_states_summing(theta, phi, indices, z_approx)
        real(dp), intent(in) :: theta(:, :)                  !< ISDF interpolation vectors of shape(Np, Nintr)
        real(dp), intent(in) :: phi(:, :)                    !< KS states of shape (Np, m_states)
        integer, intent(in)  :: indices(:)                   !< Indices of interpolation points
        real(dp), allocatable, intent(out) :: z_approx(:, :) !< Approximate product basis 

        integer :: np, m, nintr, r, i, j, ij, ir, i_mu
        real(dp) :: phi_mu_i, phi_mu_j

        np = size(phi, 1)
        m = size(phi, 2)
        nintr = size(theta, 2)

        allocate(z_approx(np, m * m))

        ! NOTE, the order of i and j matters iff phi /= psi
        ij = 0
        do i = 1, m
            do j = 1, m
                ij = ij + 1
                do i_mu = 1, nintr
                    r = indices(i_mu)
                    phi_mu_i = phi(r, i)
                    phi_mu_j = phi(r, j)
                    do ir = 1, np
                        z_approx(ir, ij) = z_approx(ir, ij) + theta(ir, i_mu) * phi_mu_i * phi_mu_j
                    enddo
                enddo
            enddo
        enddo

    end subroutine construct_approx_product_states_summing


    ! ! @brief Convert P matrix from P(r, r_u) to P(r_v, r_u)
    ! subroutine select_interpolation_points_of_first_dim(P_phi_ri, indices, P_phi_ii)
    !     real(dp),  intent(in)  :: P_phi_ri(:, :)
    !     integer,   intent(in)  :: indices(:)
    !     real(dp),  intent(out) :: P_phi_ii(:, :)

    !     integer :: i, j, ir, nintr

    !     nintr = size(indices)

    !     if (size(P_phi_ri, 2) /= nintr) then
    !         write(*, *) 'Number Interpolation point indices /= size(P_phi_ri, 2)'
    !         error stop 101
    !     endif

    !     if (all(shape(P_phi_ii) /= [nintr, nintr])) then
    !         write(*, *) 'Size of density matrix should be (nintr, nintr)'
    !         error stop 102
    !     endif

    !     !$omp parallel do simd collapse(2) default(shared) private(ir)
    !     do j = 1, nintr
    !         do i = 1, nintr
    !             ir = indices(i)
    !             P_phi_ii(i, j) = P_phi_ri(ir, j)
    !         enddo
    !     enddo
    !     !$omp end parallel do simd

    ! end subroutine select_interpolation_points_of_first_dim


    ! ! Address assigning optional arg to a pointer once first implementation tested.
    ! !> @brief Construct a matrix of ISDF interpolation vectors
    ! subroutine construct_interpolation_vectors(phi, indices, theta, psi)
    !     real(dp),          intent(in)  :: phi(:, :)      !< A set of states defined on real-space grid
    !     !                                                   of shape (np, mstates)
    !     integer,           intent(in)  :: indices(:)     !< Indices of interpolation grid points

    !     real(dp), allocatable, intent(out) :: theta(:, :)    !< Matrix of interpolating vectors
    !     integer, optional, intent(in)  :: psi(:)         !< A second set of states defined on real-space grid
    !     !                                                   of shape (np, nstates)

    !     integer                        :: np, nintr, ir, iintr, i, j
    !     ! Work arrays, where r = realspace point and i=interpolation point
    !     real(dp),  target, allocatable :: P_phi_ri(:, :), P_phi_ii(:, :)
    !     real(dp),  target, allocatable :: P_psi_ri_work(:, :), P_psi_ii_work(:, :)
    !     real(dp),          pointer     :: P_psi_ri(:, :), P_psi_ii(:, :)

    !     ! Allocate interpolation vectors matrix
    !     np = size(phi, 1)
    !     nintr = size(indices)
    !     allocate(theta(np, nintr))

    !     ! 1. Quasi density matrices
    !     call construct_quasi_density_matrix_R_intr(phi, indices, P_phi_ri)
    !     if (present(psi)) then
    !         call construct_quasi_density_matrix_R_intr(psi, indices, P_psi_ri_work)
    !         P_psi_ri => P_psi_ri_work
    !     else
    !         P_psi_ri => P_phi_ri
    !     endif
        
    !     ! 2. Construct CC^T
    !     ! a) Convert P(r, r_u) -> P(r_v, r_u)
    !     allocate(P_phi_ii(nintr, nintr))
    !     call select_interpolation_points_of_first_dim(P_phi_ri, indices, P_phi_ii)
    !     if (present(psi)) then
    !         allocate(P_psi_ii_work(nintr, nintr))
    !         call select_interpolation_points_of_first_dim(psi, indices, P_psi_ii_work)
    !         P_psi_ii => P_psi_ii_work
    !     else
    !         P_psi_ii => P_phi_ii
    !     endif

    !     ! b) Contract and inverse


    !     ! TODO Check how I get the consistency between (ZC^T) and (CC^T)^{-1} dimensions
    !     ! hence which dimension to contract over.

    !     ! 3. Construct theta = (ZC^T)(CC^T)^{-1} == P_phi_ri * P_psi_ri * [P_phi_ii * P_psi_ii]^{-1}
    !     do iintr = 1, nintr
    !         do ir = 1, np
    !             theta(ir, iintr) = P_phi_ri(ir, iintr) * P_psi_ri(ir, iintr)
    !         end do
    !     enddo

    !     deallocate(P_phi_ri)
    !     if(allocated(P_psi_work)) deallocate(P_psi_work)
    !     nullify(P_psi_ri)
    !     deallocate(cc_inv)

    ! end subroutine construct_interpolation_vectors


end module isdf_serial_m
