module grid_funcs_m
    use, intrinsic :: iso_fortran_env, only: dp => real64

    ! From k-means library
    use grids_m, only: linspace, linspace_to_grid

    implicit none
    private

    type GridParams3D
        integer :: points(3)
        real(dp) :: origin(3)
        real(dp) :: limits(3)
        real(dp) :: spacings(3, 3)
        real(dp), allocatable :: grid(:, :)
        integer :: np
    contains 
        procedure :: init => init_grid_params
    end type

    public :: GridParams3D, construct_sin_basis_with_grid

contains

    subroutine init_grid_params(this, points, origin, limits)
        class(GridParams3D), intent(inout) :: this
        integer , intent(in) :: points(3)
        real(dp), intent(in) :: origin(3)
        real(dp), intent(in) :: limits(3)

        real(dp), allocatable :: x(:), y(:), z(:)

        this%points = points
        this%origin = origin
        this%limits = limits
        this%np = product(this%points)

        allocate(x(this%points(1)), y(this%points(2)), z(this%points(3)))
        allocate(this%grid(3, this%np))

        call linspace(this%origin(1), this%limits(1), this%points(1), x)
        call linspace(this%origin(2), this%limits(2), this%points(2), y)
        call linspace(this%origin(3), this%limits(3), this%points(3), z)
        
        call linspace_to_grid(x, y, z, this%grid)

        this%spacings = 0._dp
        this%spacings(1, 1) = x(points(1)) - x(1)
        this%spacings(2, 2) = y(points(2)) - y(1)
        this%spacings(3, 3) = z(points(3)) - z(1)

    end subroutine init_grid_params


   ! TODO Move this to a module within this library
   ! Assumes that the grid is stored in a cartesian format
    subroutine construct_sin_basis_with_grid(grid, sin_pns, func)
        real(dp), intent(in)  :: grid(:, :)   !> Cartesian grid
        integer,  intent(in)  :: sin_pns(:)   !> Principal numbers for each sin function
        real(dp), intent(out) :: func(:)      !> Product of sins on the grid

        integer :: n, m, l, ir, np
        real(dp) :: x, y, z, Lx, Ly, Lz, norm
        real(dp), parameter :: pi = 3.14159265359_dp

        ! Assumes first grid point is the origin corner, and the last grid point is the further corner
        ! from the initial point
        np = size(grid, 2)
        Lx = grid(1, np) - grid(1, 1)
        Ly = grid(2, np) - grid(2, 1)
        Lz = grid(3, np) - grid(3, 1)

        n = sin_pns(1)
        m = sin_pns(2)
        l = sin_pns(3)
        norm = sqrt(2._dp / Lx) * sqrt(2._dp / Ly) * sqrt(2._dp / Lz)

        do ir = 1, np
            x = grid(1, ir)
            y = grid(2, ir)
            z = grid(3, ir)
            ! Feel like I'm missing a factor 2 from args, but plot looks ok
            func(ir) = norm * sin(n * pi * x / Lx) * sin(m * pi * y / Ly) * sin(l * pi * z / Lz)
        enddo

    end subroutine construct_sin_basis_with_grid

end module grid_funcs_m