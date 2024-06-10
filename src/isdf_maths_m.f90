module isdf_maths_m
   use, intrinsic :: iso_fortran_env, only: dp => real64
   implicit none
   private

   ! Exposed routines and functions
   public :: close, all_close, gram_schmidt, modified_gram_schmidt, construct_sin_basis_with_grid

   interface all_close
      module procedure :: all_close_real64_1d, all_close_real64_2d
   end interface all_close

   real(dp), parameter :: pi = 3.14159265359_dp

contains

   !> @brief Are \f$x\f$ and \f$y\f$ equal within a tolerance.
   elemental logical function close (x, y, rtol, atol)
      real(dp), intent(in) :: x, y
      real(dp), optional, intent(in) :: rtol
      real(dp), optional, intent(in) :: atol
      real(dp) :: atol_, rtol_

      if (present(rtol)) then
         rtol_ = rtol
      else
         rtol_ = 1.e-5_dp
      end if

      if (present(atol)) then
         atol_ = atol
      else
         atol_ = 1.e-8_dp
      end if

      close = abs(x - y) <= (atol_ + rtol_*abs(y))
   end function close


   ! Cannot be elemental when shape of inputs differ from shape of output
   logical function all_close_real64_1d(x, y, rtol, atol)
      real(dp), intent(in) :: x(:), y(:)
      real(dp), optional, intent(in) :: rtol
      real(dp), optional, intent(in) :: atol
      real(dp) :: atol_, rtol_
      logical, allocatable :: values_close(:)

      if (present(rtol)) then
         rtol_ = rtol
      else
         rtol_ = 1.e-5_dp
      end if

      if (present(atol)) then
         atol_ = atol
      else
         atol_ = 1.e-8_dp
      end if

      ! Explicitly allocate else I can the compiler warning:
      ! `'values_close.offset' is used uninitialized [-Wuninitialized]`
      allocate (values_close(size(x)))
      values_close = abs(x - y) <= (atol_ + rtol_*abs(y))
      all_close_real64_1d = all(values_close)

   end function all_close_real64_1d


   logical function all_close_real64_2d(x, y, rtol, atol)
      real(dp), intent(in) :: x(:, :), y(:, :)
      real(dp), optional, intent(in) :: rtol
      real(dp), optional, intent(in) :: atol
      real(dp) :: atol_, rtol_
      logical, allocatable :: values_close(:, :)

      if (present(rtol)) then
         rtol_ = rtol
      else
         rtol_ = 1.e-5_dp
      end if

      if (present(atol)) then
         atol_ = atol
      else
         atol_ = 1.e-8_dp
      end if

      ! Explicitly allocate else I can the compiler warning:
      ! `'values_close.offset' is used uninitialized [-Wuninitialized]`
      allocate (values_close(size(x, 1), size(x, 2)))
      values_close = abs(x - y) <= (atol_ + rtol_*abs(y))
      all_close_real64_2d = all(values_close)

   end function all_close_real64_2d


   !>@brief Sum of the projection of vector v_i onto a set of vectors {u}.
   subroutine summed_projection(v, u, v_index, proj)
      real(dp), intent(in), contiguous :: v(:)
      real(dp), intent(in), contiguous :: u(:, :)
      integer, intent(in) :: v_index
      real(dp), intent(out) :: proj(:)

      integer :: i

      ! Assumes U in [1, v_index - 1] are normalised in the caller
      ! TODO(Alex) See if I can replace the dot product with blas
      ! TODO(Alex) Consider OMP parallelisation
      proj = dot_product(v, u(:, 1))*u(:, 1)

      do i = 2, v_index - 1
         proj = proj + dot_product(v, u(:, i))*u(:, i)
      end do

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
      allocate (proj(m))

      ! TODO(Alex) See if I can replace norm call with blas
      norm = norm2(v(:, 1))
      if(close(norm, 0._dp)) then
         write(*, *)'Norm of the first vector is < 1.e-8'
         error stop
      endif
      call dscal(m, 1._dp/norm, v(:, 1), 1)

      do i = 2, n_vectors
         call summed_projection(v(:, i), v, i, proj)
         ! v(:, i) = v(:, i) - proj(:)
         call daxpy(m, -1._dp, proj, 1, v(:, i), 1)
         norm = norm2(v(:, i))
         ! Handle linearly-dependent vectors
         if (close(norm, 0._dp)) then
            v(:, i) = 0._dp
         else
            call dscal(m, 1._dp/norm, v(:, i), 1)
         end if
      end do

   end subroutine gram_schmidt


   !> @brief Orthogonalisation of column vectors with modified Gram-Schmidt.
   !!
   !! The projection of each vector v(:,j) is subtracted from only the remaining vectors
   !! (v(:,j+1) to v(:,n)). This orthogonalizes the remaining vectors with respect to the
   !! newly orthogonalized v(:,j), and results in a more numerically stable algorithm than the
   !! classical version, especially for ill-conditioned matrices.
   !!
   !! Also see https://laurenthoeltgen.name/post/gram-schmidt/
   subroutine modified_gram_schmidt(v, verbose)
      real(dp), intent(inout), contiguous :: v(:, :) !< In: Array of column vectors
      !                                                Out: Orthogonalised column vectors
      logical, intent(in), optional :: verbose
      logical  :: print_out
      integer  :: n_vectors, m, i, j
      real(dp) :: norm

      m = size(v, 1)
      n_vectors = size(v, 2)
      print_out = .false.
      if(present(verbose)) print_out = verbose

      do j = 1, n_vectors
         norm = norm2(v(:, j))
         ! Handle linearly-dependent vectors
         if (close (norm, 0._dp)) then
            if(print_out) write(*, *) 'Vector j is linearly-dependent, hence zeroed', j
            v(:, j) = 0._dp
            cycle
         end if
         call dscal(m, 1._dp/norm, v(:, j), 1)
         do i = j + 1, n_vectors
            v(:, i) = v(:, i) - dot_product(v(:, j), v(:, i))*v(:, j)
         end do
      end do

   end subroutine modified_gram_schmidt


   ! subroutine construct_sin_basis(n_points, origin, end, sin_pns, func)
   !    integer,  intent(in)  :: n_points(:)  !> Points per dimension
   !    real(dp), intent(in)  :: end(:)       !> Limits
   !    real(dp), intent(in)  :: origin(:)    !> Origin
   !    integer,  intent(in)  :: sin_pns(:)   !> Principal numbers for each sin function
   !    real(dp), intent(out) :: func(:)      !> Product of sins on the grid

   !    integer :: nx, ny, nz, ix, iy, iz, n, m, l, ir
   !    real(dp) :: x, y, z, x0, y0, z0, dx, dy, dz, Lx, Ly, Lz

   !    nx = n_points(1)
   !    ny = n_points(2)
   !    nz = n_points(3)

   !    x0 = origin(1)
   !    y0 = origin(2)
   !    z0 = origin(3)

   !    Lx = (end(1) - x0) 
   !    Ly = (end(2) - y0) 
   !    Lz = (end(3) - z0)

   !    ! Assumes end point is included 
   !    dx = Lx / real(nx - 1, dp)
   !    dy = Ly / real(ny - 1, dp)
   !    dz = Lz / real(nz - 1, dp)

   !    n = sin_pns(1)
   !    m = sin_pns(2)
   !    l = sin_pns(3)

   !    ir = 0
   !    do iz = 1, nz
   !       z = z0 + (iz - 1) * dz
   !       do iy = 1, ny
   !          y = y0 + (iy - 1) * dy
   !          do ix = 1, nx
   !             x = x0 + (ix - 1) * dx
   !             ir = ir + 1
   !             func(ir) = sin(n * pi * x / Lx) * sin(m * pi * y / Ly) * sin(l * pi * z / Lz)
   !          enddo
   !       enddo
   !    enddo

   ! end subroutine construct_sin_basis


   ! Assumes that the grid is stored in a cartesian format
   subroutine construct_sin_basis_with_grid(grid, sin_pns, func)
      real(dp), intent(in)  :: grid(:, :)   !> Cartesian grid
      integer,  intent(in)  :: sin_pns(:)   !> Principal numbers for each sin function
      real(dp), intent(out) :: func(:)      !> Product of sins on the grid

      integer :: n, m, l, ir, np
      real(dp) :: x, y, z, Lx, Ly, Lz

      ! Assumes first grid point is the origin corner, and the last grid point is the further corner
      ! from the initial point
      np = size(grid, 2)
      Lx = grid(1, np) - grid(1, 1)
      Ly = grid(2, np) - grid(2, 1)
      Lz = grid(3, np) - grid(3, 1)

      n = sin_pns(1)
      m = sin_pns(2)
      l = sin_pns(3)

      do ir = 1, np
         x = grid(1, ir)
         y = grid(2, ir)
         z = grid(3, ir)
         func(ir) = sin(n * pi * x / Lx) * sin(m * pi * y / Ly) * sin(l * pi * z / Lz)
      enddo

   end subroutine construct_sin_basis_with_grid

end module isdf_maths_m
