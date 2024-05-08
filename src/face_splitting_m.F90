! TODO(alex). Look at kronecker products offered by lapack
! Look at scalapack
module face_splitting_m
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: face_splitting_product, &
              face_splitting_product_two_funcs_rowwise, &
              face_splitting_product_two_funcs_colwise

    interface face_splitting_product
        module procedure :: face_splitting_product_one_func, face_splitting_product_two_funcs_rowwise
    end interface face_splitting_product

contains

    !> @brief Face-splitting Product. Create products of KS wave functions on a real-space grid.
    !!
    !! See [wikipedia](https://en.wikipedia.org/wiki/Khatri–Rao_product#Face-splitting_product) entry and associated
    !! image for a clear description of the algorithm
    !!
    !! This type of operation is a row-by-row Kronecker products of two matrices.
    subroutine face_splitting_product_one_func(phi, z)
        real(dp),              intent(in)  :: phi(:, :)   !< Set of wave functions on a real-space grid
        real(dp), allocatable, intent(out) :: z(:, :)     !< Face-split product

        integer  :: nrow, ncol, icol, jcol, ij

        nrow = size(phi, 1)
        ncol = size(phi, 2)
        allocate(z(nrow, ncol * ncol))

        !$omp parallel do simd default(shared) private(icol, jcol)
        do ij = 1, ncol * ncol
            icol = ((ij - 1) / ncol) + 1
            jcol = ij - (icol - 1) * ncol
            z(:, ij) = phi(:, icol) * phi(:, jcol)
        enddo
        !$omp end parallel do simd

    end subroutine face_splitting_product_one_func


    !> @brief Face-splitting Product. Create products of KS wave functions on a real-space grid.
    !!
    !! This type of operation is a row-by-row Kronecker products of two matrices.
    !! Use this routine if phi and psi have shape (n_grid, n_states) and (n_grid, m_states), respectively, 
    !! and one wishes to return a product of shape (n_grid, n_states * m_states).
    !!
    !! See [wikipedia](https://en.wikipedia.org/wiki/Khatri–Rao_product#Face-splitting_product) entry and associated
    !! image for a clear description of the algorithm.
    !!
    subroutine face_splitting_product_two_funcs_rowwise(phi, psi, z)
        real(dp), intent(in)  :: phi(:, :)                !< Set of wave functions on a real-space grid
        real(dp), intent(in)  :: psi(:, :)                !< Second set of wave functions on a real-space grid
        real(dp), allocatable, intent(out) :: z(:, :)     !< Face-split product

        integer  :: nrow, ncol, mcol, icol, jcol, ij

        ! Declarations specific to 1st and 2nd implementations
        ! integer  :: irow, , ij
        ! real(dp), allocatable :: phi_column(:)

        if (all(shape(psi) /= shape(phi))) then
            write(*, *) 'Shape of first argument does not matcht the shape of the second'
            error stop
        endif

        nrow = size(phi, 1)
        ncol = size(phi, 2)
        mcol = size(psi, 2)
        allocate(z(nrow, ncol * mcol))

        ! Implementation with 3 loops
        ! allocate(phi_column(nrow))
        ! ij = 0
        ! do icol = 1, ncol
        !     phi_column = phi(:, icol)
        !     do jcol = 1, mcol
        !         ij = ij + 1
        !         do irow = 1, nrow
        !             z(irow, ij) = phi_column(irow) * psi(irow, jcol)
        !         enddo
        !     enddo
        ! enddo

        ! Implementation with 2 loops.
        ! allocate(phi_column(nrow))
        ! ij = 0
        ! do icol = 1, ncol
        !     phi_column = phi(:, icol)
        !     do jcol = 1, mcol
        !         ij = ij + 1
        !         z(:, ij) = phi_column(:) * psi(:, jcol)
        !     enddo
        ! enddo

        ! Single loop OMP implementation
        !$omp parallel do simd default(shared) private(icol, jcol)
        do ij = 1, ncol * mcol
            ! Loop unrolling assumes outer loop is icol, and inner loop is jcol (consistent with above)
            icol = ((ij - 1) / mcol) + 1
            jcol = ij - (icol - 1) * mcol
            z(:, ij) = phi(:, icol) * psi(:, jcol)
        enddo
        !$omp end parallel do simd

    end subroutine face_splitting_product_two_funcs_rowwise


    !> @brief Face-splitting Product. Create products of KS wave functions on a real-space grid.
    !!
    !! This type of operation is a column-by-column Kronecker products of two matrices.
    !! Use this routine if phi and psi have shape (n_states, n_grid) and (m_states, n_grid), respectively, 
    !! and one wishes to return a product of shape (n_states * m_states, n_grid).
    subroutine face_splitting_product_two_funcs_colwise(phi, psi, z)
        real(dp), intent(in)  :: phi(:, :)                !< Set of wave functions on a real-space grid
        real(dp), intent(in)  :: psi(:, :)                !< Second set of wave functions on a real-space grid
        real(dp), allocatable, intent(out) :: z(:, :)     !< Face-split product

        integer  :: icol, irow, ncol, nrow, mrow

        ! Declarations specific to first implementation
        ! integer  :: jrow, ij
        ! real(dp) :: phi_element

        ! Declarations specific to second implementation
        integer  :: ij_start, ij_end
        real(dp), allocatable :: psi_column(:)

        if (all(shape(psi) /= shape(phi))) then
            write(*, *) 'Shape of first argument does not matcht the shape of the second'
            error stop
        endif

        nrow = size(phi, 1)
        mrow = size(psi, 1)
        ncol = size(phi, 2)
        allocate(z(nrow * mrow, ncol))

        ! Implementation with 3 loops
        ! do icol = 1, ncol
        !     ! Kroneckor product
        !     ij = 0
        !     do irow = 1, nrow
        !         phi_element = phi(irow, icol)
        !         do jrow = 1, mrow
        !             ij = ij + 1
        !             z(ij, icol) = phi_element * psi(jrow, icol)
        !         enddo
        !     enddo
        ! enddo

        ! Implementation with 2 loops
        allocate(psi_column(mrow))
        do icol = 1, ncol
            psi_column = psi(:, icol)
            ! Kroneckor product on two column vectors
            do irow = 1, nrow
                ij_start = 1 + (irow - 1) * mrow
                ij_end = irow * mrow
                z(ij_start:ij_end, icol) = phi(irow, icol) * psi_column
            enddo
        enddo

    end subroutine face_splitting_product_two_funcs_colwise

end module face_splitting_m
