module helpers
    use, intrinsic :: iso_fortran_env
    implicit none
    private

    public :: load_numpy_array, mean_square_error, load_centroid_points, discretise_values_to_grid, &
        load_centroid_indices, write_numpy_array

    interface load_numpy_array
        module procedure :: load_numpy_array_2d
    end interface

    interface write_numpy_array
        module procedure :: write_numpy_array_2d
    end interface

contains

function mean_square_error(x, y) result(mse)
    real(real64), intent(in) :: x(:)
    real(real64), intent(in) :: y(:)
    real(real64) :: mse
  
    integer :: n, i
  
    n = size(x)
    if (n /= size(y)) then
       write(*, *) 'Size of x and y disagree'
       error stop 
    endif
  
    mse = 0._real64
    do i = 1, n
      mse = mse + (x(i) - y(i))**2._real64
    enddo
    mse = mse / real(n)
  
  end function mean_square_error
  

    subroutine load_centroid_indices(fname, indices)
        character(len=*), intent(in) :: fname
        integer, allocatable, intent(out) :: indices(:)
        integer :: n_centroids, i

        open(unit=101, file=trim(fname))
        read(101, *) n_centroids
        allocate(indices(n_centroids))

        do i = 1, n_centroids
            read(101, *) indices(i)
            ! Convert from python to fortran indexing
            indices(i) = indices(i) + 1
        enddo

        close(101)

    end subroutine load_centroid_indices


    subroutine load_centroid_points(fname, centroid_points)
        character(len=*), intent(in) :: fname
        real(real64), allocatable, intent(out) :: centroid_points(:, :)
        integer :: n_centroids, i

        open(unit=101, file=trim(fname))
        read(101, *) n_centroids
        allocate(centroid_points(3, n_centroids))

        do i = 1, n_centroids
            read(101, *) centroid_points(:, i)
        enddo

        close(101)

    end subroutine load_centroid_points


    subroutine load_numpy_array_2d(fname, array)
        character(len=*), intent(in) :: fname
        real(real64), allocatable, intent(out) :: array(:, :)

        integer :: n_row, m_col, i

        open(unit=101, file=trim(adjustl(fname)))
        read(101, *) n_row, m_col

        ! Inefficient memory access, but ensures same shape
        ! as python data. One could alternatively parse in col-major,
        ! then transpose
        allocate(array(n_row, m_col))
        do i = 1, n_row
            read(101, *) array(i, :)
        enddo

        close(101)

    end subroutine load_numpy_array_2d


    ! Write row by row - inefficient but allows for simple parsing and comparison
    subroutine write_numpy_array_2d(fname, array)
        character(len=*), intent(in) :: fname
        real(real64), intent(in) :: array(:, :)

        integer :: n_row, m_col, i

        open(unit=101, file=trim(adjustl(fname)))

        n_row = size(array, 1)
        m_col = size(array, 2)

        write(101, *) n_row, m_col

        do i = 1, n_row
            write(101, *) array(i, :)
        enddo

        close(101)

    end subroutine write_numpy_array_2d


    !> @brief Discretise continious values with discrete grid points.
    !!
    !! Perform a linear search through all grid points and replace each continuous
    !! value with the closest discrete point. Note, this routine uses the metric 
    !! ||r - r'||^2 rather than ||r - r'||, as it is faster to evaluate.
    subroutine discretise_values_to_grid(values, grid, indices)
        real(real64), intent(inout) :: values(:, :)     !< In: Continuous values  (n_dim, M) 
        !                                                Out: Discretised values (n_dim, M) 
        real(real64), intent(in)    :: grid(:, :)       !< Discrete grid (n_dim, Nr)
        integer,      intent(out), optional :: indices(:) !< Indices of discretised centroids on grid
        
        integer :: ir, iv, nr, iopt, n_dim
        real(real64) :: norm, new_norm
        real(real64), allocatable :: val(:)

        n_dim = size(grid, 1)
        nr = size(grid, 2)
        allocate(val(n_dim))

        do iv = 1, size(values, 2)
            val = values(:, iv)
            ! Initialise 
            norm = sum((grid(:, 1) - val(:))**2)
            iopt = 1
            do ir = 2, nr
                new_norm = sum((grid(:, ir) - val(:))**2)
                if (new_norm < norm) then
                    norm = new_norm
                    iopt = ir
                endif
            enddo
            values(:, iv) = grid(:, iopt)
            if (present(indices)) indices(iv) = iopt
        enddo

    end subroutine discretise_values_to_grid

end module helpers
