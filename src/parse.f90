module parse
    use, intrinsic :: iso_fortran_env
    implicit none
    private

    integer, parameter :: new_line_char = 10  !< ASCII

    public :: write_to_xyz, output_cube, parse_grid_1d_from_c, parse_grid_2d_from_c, parse_from_cube

contains
    
    ! subroutine parse_text_1d(fname, data)
    !     character(len=*), intent(in) :: fname
    !     real(real64), allocatable, intent(out) :: data(:)

    !     character(len=1) :: dummy
    !     integer :: n, ierr, i

    !     open(unit=101, file=trim(fname), form='formatted', access='stream', iostat=ierr)

    !     if (ierr /= 0) then
    !         write(*, *) "Error opening file."
    !         stop
    !     end if

    !     ! Parse shape from the header
    !     read(101, *) dummy, n

    !     allocate(data(n))
    !     do i = 1, n
    !         read(101, *) data(i)
    !     enddo

    !     close(101)

    ! end subroutine parse_text_1d

    subroutine parse_grid_1d_from_c(fname, data)
        character(len=*), intent(in) :: fname
        real(real64), allocatable, intent(out) :: data(:)

        character(len=1) :: dummy
        integer :: nx, ny, nz, ierr, ix, iy, iz
        real(real64), allocatable :: tmp(:, :, :)

        open(unit=101, file=trim(fname), form='formatted', access='stream', iostat=ierr)

        if (ierr /= 0) then
            write(*, *) "Error opening file."
            stop
        end if

        ! Parse shape from the header
        ! Expect 3D grid
        read(101, *) dummy, nx, ny, nz

        allocate(tmp(nx, ny, nz))

        ! C-style memory access so the order is consistent in fortran
        do ix = 1, nx
            do iy = 1, ny
                do iz = 1, nz
                    read(101, *) tmp(ix, iy, iz)
                enddo
            enddo
        enddo
        close(101)

        data = reshape(tmp, [nx * ny * nz])

    end subroutine parse_grid_1d_from_c


    subroutine parse_grid_2d_from_c(fname, packed, data)
        character(len=*), intent(in) :: fname   !< File name
        logical                      :: packed  !< Octopus convention. If packed, (istate, ip), Not packed (ip, istate)
        real(real64), allocatable, intent(out) :: data(:, :)  !< Function/s on grid

        character(len=1) :: dummy
        integer :: n_dim, nx, ny, nz, ierr, ix, iy, iz
        real(real64), allocatable :: tmp(:, :, :, :)

        open(unit=101, file=trim(fname), form='formatted', access='stream', iostat=ierr)

        if (ierr /= 0) then
            write(*, *) "Error opening file."
            stop
        end if

        ! Parse shape from the header. NOTE, this is order is a hard-coded convention
        ! Expect 3D grid
        read(101, *) dummy, nx, ny, nz, n_dim

        if (.not. packed) then
            allocate(tmp(nx, ny, nz, n_dim))

            ! C-style memory access so the order is consistent in fortran
            do ix = 1, nx
                do iy = 1, ny
                    do iz = 1, nz
                        read(101, *) tmp(ix, iy, iz, :)
                    enddo
                enddo
            enddo
            close(101)

            data = reshape(tmp, [(nx * ny * nz), n_dim])
        else
            allocate(tmp(n_dim, nx, ny, nz))

            ! C-style memory access so the order is consistent in fortran
            do ix = 1, nx
                do iy = 1, ny
                    do iz = 1, nz
                        read(101, *) tmp(:, ix, iy, iz)
                    enddo
                enddo
            enddo
            close(101)

            data = reshape(tmp, [n_dim, (nx * ny * nz)])
        endif
        deallocate(tmp)

    end subroutine parse_grid_2d_from_c


    ! !> @brief Parse a binary file, where the first line contains the shape of the array
    ! !!
    ! !! Expects a file of the form:
    ! !!
    ! !! ```
    ! !! 1000
    ! !! �Uԁ���X�{J�\{����
    ! !!```
    ! subroutine parse_function_on_grid_1d(fname, data)
    !     character(len=*), intent(in) :: fname
    !     real(real64), allocatable, intent(out) :: data(:)

    !     integer :: n, ierr
    !     character(len=1) :: char_val

    !     ! Parse shape from the header
    !     open(unit=101, file=trim(fname), form='formatted', access='stream', iostat=ierr)

    !     if (ierr /= 0) then
    !         write(*, *) "Error opening file."
    !         stop
    !     end if

    !     read(101, *) n
    !     close(101)
    !     allocate(data(n))

    !     ! Read the remaining binary data (now in unformatted mode)
    !     open(unit=101, file=trim(fname), form='unformatted', access='stream', iostat=ierr)

    !     ! Skip the first line of the file
    !     rewind(101)
    !     do while (.true.)
    !         ! Read a single byte in unformatted mode
    !         read(101) char_val
    !         if (ichar(char_val) == new_line_char) exit   ! Check for newline (ASCII 10)
    !     end do

    !     ! Parse binary stream into array
    !     read(101) data
    !     close(101)

    ! end subroutine parse_function_on_grid_1d


    ! subroutine parse_function_on_grid_2d(fname, data)
    !     character(len=*), intent(in) :: fname
    !     real(real64), allocatable, intent(out) :: data(:, :)

    !     integer :: n, m, ierr
    !     character(len=1) :: char_val

    !     ! Parse shape from the header
    !     open(unit=101, file=trim(fname), form='formatted', access='stream', iostat=ierr)

    !     if (ierr /= 0) then
    !         write(*, *) "Error opening file."
    !         stop
    !     end if

    !     read(101, *) n, m
    !     close(101)
    !     allocate(data(n, m))
    !     write(*, *) trim(fname), n, m

    !     ! Read the remaining binary data (now in unformatted mode)
    !     open(unit=101, file=trim(fname), form='unformatted', access='stream', iostat=ierr)

    !     ! Skip the first line of the file
    !     rewind(101)
    !     do while (.true.)
    !         ! Read a single byte in unformatted mode
    !         read(101) char_val
    !         if (ichar(char_val) == new_line_char) exit   ! Check for newline (ASCII 10)
    !     end do

    !     ! Parse binary stream into array
    !     read(101) data
    !     close(101)

    ! end subroutine parse_function_on_grid_2d


    !> @brief Write xyz file for finite systems.
    subroutine write_to_xyz(file, species, positions)
        character(len=*), intent(in) :: file
        character(len=*), intent(in) :: species(:) 
        real(real64),     intent(in) :: positions(:, :)
        integer :: natoms, i

        natoms = size(positions, 2)
        if (size(species) /= natoms) then
            write(*, *) 'size(species) == natoms'
        endif

        open(unit=101, file=trim(file))

        write(101, *) ''
        write(101, *) natoms
        do i = 1, natoms
            write(101, *) species(i), positions(:, i)
        enddo

        close(101)

    end subroutine write_to_xyz



    !> @brief Parse a cube file into a 3D array
    !!
    !! Also return the grid data, for grid reconstruction
    subroutine parse_from_cube(fname, an, position, n_points, spacing, origin, data)
        character(len=*),          intent(in)  :: fname              !< File name (minus extension)
        integer,      allocatable, intent(out) :: an(:)
        real(real64), allocatable, intent(out) :: position(:, :)
        integer,                   intent(out) :: n_points(:)        !< N points in grid
        real(real64),              intent(out) :: spacing(:, :)      !< Grid spacings, stored row-wise
        real(real64),              intent(out) :: origin(:)          !< First grid point
        real(real64), allocatable, intent(out) :: data(:, :, :)      !< (ix, iy, iz)

        integer, parameter        :: row_len=6                       !< Row length fixed by cube file format
        integer                   :: n_atoms
        integer                   :: i, ia, nrow_z, n_remainder, rowz, ix, iy, iz1, iz2, iremainder
        real(real64)              :: rdummy

        open(unit=001, file=trim(adjustl(fname))//'.cube')

        ! Skip lines 1 and 2
        read(001, *)
        read(001, *)

        ! TODO. Should extract the format specifiers from octopus
        !read(001, '(i5, X, 3(F12.6, X))') n_atoms, origin(:)
        read(001, *) n_atoms, origin(:)

        ! Grid points and spacings
        do i = 1, 3
            !read(001,'(i5, X, 3(F12.6, X))') n_points(i), spacing(i, :)
            read(001, *) n_points(i), spacing(i, :)
        enddo

        ! Atomic number, dummy charge and position
        allocate(position(3, n_atoms), an(n_atoms))
        do ia = 1, n_atoms
            !read(001, '(i5, X, 4(F12.6, X))') idummy, rdummy, position(:, ia)
            read(001, *) an(ia), rdummy, position(:, ia)
        enddo

        !Remaining lines = Volumetric data
        ! Number of full rows
        nrow_z = int(n_points(3) / row_len)

        !If not divisible by 6, how many values will be put on the last row:
        n_remainder = mod(n_points(3), row_len)

        allocate(data(n_points(1), n_points(2), n_points(3))) 

        ! No simple way around poor memory access
        do ix = 1, n_points(1)
            do iy = 1, n_points(2)
               ! Full rows 
               do rowz = 1, nrow_z
                    iz1 = 1 + ((rowz-1)* row_len)
                    iz2 = row_len + ((rowz-1)* row_len)
                    read(001, *) data(ix, iy, iz1:iz2)
                    !read(001,'(6E13.5)') data(ix, iy, iz1:iz2) Did not work
               enddo
               !Remainder row, which does not get touched if n_remainder == 0
               do iremainder = 1, n_remainder
                read(001,'(E13.5)', advance="no") data(ix, iy, iz2 + iremainder)
               enddo
               if (n_remainder > 0) read(001, *)
            enddo
        enddo 

        ! NOTE< flattening this index is bad - won't be consistent, hence why using the above
        ! ie, want from inner to outer, ix, iy, iz
        ! do ix = 1, n_points(1)
        !     do iy = 1, n_points(2)
        !         base_index = (iy - 1) * n_points(2) + (ix - 1) * n_points(2) * n_points(1)
        !        ! Full rows 
        !        do rowz = 1, nrow_z
        !             iz1 = 1 + ((rowz-1)* row_len)
        !             iz2 = row_len + ((rowz-1)* row_len)
        !             ir1 = iz1 + base_index
        !             ir2 = iz2 + base_index
        !             read(001,'(6E13.5)') data(ir1:ir2)
        !        enddo
        !        !Remainder row, which does not get touched if n_remainder == 0
        !        do iremainder = 1, n_remainder
        !           read(001,'(E13.5)', advance="no") data(ir2 + iremainder)
        !        enddo
        !        if (n_remainder > 0) read(001, *)
        !     enddo
        ! enddo 

        close(001)

    end subroutine parse_from_cube


 
    !> @brief Output function defined on grid, to Gaussian cube file format.
    !!
    !! Output looks like:
    !! ```cube
    !! Header line 1
    !! Header line 2
    !! natoms origin_x origin_y origin_z
    !! ni  i_x  i_y  i_z                         # N grid points in i direction, plus spacings
    !! nj  j_x  j_y  j_z
    !! nk  k_x  k_y  k_z
    !! AN  0.00000  x_atom1  y_atom1  z_atom1    # Atomic number   dble(0)  atomic position
    !! AN  0.00000  x_atom2  y_atom2  z_atom2
    !! ..
    !! AN  0.00000  x_Natom  y_Natom  z_Natom    # Voxel points
    !! 2.03274E-08  5.08882E-08  1.00820E-07  1.58068E-07  1.96097E-07  1.92495E-07 
    !! 1.49515E-07  9.18925E-08  4.46919E-08  1.72014E-08
    !! ...
    !! ```
    !!
    !! Grid should be in Bohr.
    subroutine output_cube(fname, atomic_nums, position, n_points, spacing, origin, data)
        character(len=*), intent(in) :: fname              !< File name (minus extension)
        integer,          intent(in) :: atomic_nums(:)     !< Atomic numbers of all atoms in cell
        real(real64),     intent(in) :: position(:, :)     !< Position of all atoms in the cell
        integer,          intent(in) :: n_points(:)        !< N points in grid
        real(real64),     intent(in) :: spacing(:, :)      !< Grid spacings, stored row-wise
        real(real64),     intent(in) :: origin(:)          !< First grid point
        real(real64),     intent(in), target :: data(:)    !< Assumed it is flattened (ix, iy, iz)

        integer, parameter :: row_len=6    !< Row length fixed by cube file format
        integer            :: n_atoms, nrow_z, n_remainder
        integer            :: ia, i, iremainder, ix, iy, rowz, iz1, iz2
        real(real64)       :: charge
        real(real64), pointer, dimension(:, :, :) :: data_3d

        if (product(n_points) /= size(data)) then
            write(*, *) 'Cube: Size of function does not match the total number of grid points'
            error stop 101
        endif

        n_atoms = size(atomic_nums)

        open(unit=001, file=trim(adjustl(fname))//'.cube')

        write(001, '(A)') 'Header line 1'
        write(001, '(A)') 'Header line 2'
        write(001, '(i5, X, 3(F12.6, X))') n_atoms, origin

        ! Grid points and spacings
        do i = 1, 3
            write(001,'(i5, X, 3(F12.6, X))') n_points(i), spacing(i, :)
        enddo

        ! Atomic number, dummy charge and position
        charge = 0._real64
        do ia = 1, n_atoms
            write(001, '(i5, X, 4(F12.6, X))') atomic_nums(ia), charge, position(:, ia)
        enddo

        !Remaining lines = Volumetric data
        ! Number of full rows
        nrow_z = int(n_points(3) / row_len)

        !If not divisible by 6, how many values will be put on the last row:
        n_remainder = mod(n_points(3), row_len)

        data_3d(1:n_points(1), 1:n_points(2), 1:n_points(3)) => data 
        ! No simple way around poor memory access
        do ix = 1, n_points(1)
            do iy = 1, n_points(2)
               ! Full rows 
               do rowz = 1, nrow_z
                    iz1 = 1 + ((rowz-1)* row_len)
                    iz2 = row_len + ((rowz-1)* row_len)
                    write(001,'(6E13.5)') data_3d(ix, iy, iz1:iz2)
               enddo
               !Remainder row, which does not get touched if n_remainder == 0
               do iremainder = 1, n_remainder
                    write(001,'(E13.5)', advance="no") data_3d(ix, iy, iz2 + iremainder)
               enddo
               if (n_remainder > 0) write(001, *)
            enddo
        enddo 

        close(001)

        nullify(data_3d)

    end subroutine output_cube


    ! Given any function defined on the grid (ix, iy, iz)
    ! write out (ix, iy, iz) in an order such that when read by C
    ! with a single index ir, ir correctly maps to C ordering of (ix, iy, iz)
    !
    ! ASSUMES first index is the grid index
    ! subroutine write_grid_function_to_c()

    !     do j = 1, m
    !     ! Loop over (iz, iy, ix) according to C-ordering
    !     ! Means that access of a will be slow, but ir will be consistent with what C-code expects
    !     ir  = 0
    !     do ix = 1, nx
    !         do iy = 1, ny
    !             do iz = 1, nz
    !                 ir = ir + 1
    !                 write(*, *) a(ir, j)
    !             enddo
    !         enddo
    !     enddo

    ! end subroutine write_grid_function_to_c
    
end module parse
