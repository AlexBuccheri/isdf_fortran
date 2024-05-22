module sampling_m
    use, intrinsic :: iso_fortran_env
    implicit none
    private
    public :: choose_initial_centroids_simple

contains


    ! TODO(Alex) This implement expects weight to NOT be distributed spatially
    ! i.e it''s serial-only
    ! This seems like the most straighforward implementation
    subroutine select_finite_weight(weight, indices, threshold)
        real(real64),            intent(in)  :: weight(:)  !< Weights (np)
        integer,    allocatable, intent(out) :: indices(:) !< Indices of grid points with finite weights
        real(real64),  optional, intent(in)  :: threshold  !< Threshold for finite weight

        real(real64), parameter :: default_threshold = 1.e-8_real64
        integer                 :: ir, cnt, np
        real(real64)            :: thres
        integer, allocatable    :: tmp_indices(:)

        if (present(threshold)) then
            thres = threshold
        else
            thres = default_threshold
        endif

        ! Allocate with upper bound
        np = size(weight)
        allocate(tmp_indices(np))

        cnt = 0
        do ir = 1, np
            if (abs(weight(ir)) >= thres) then
                cnt = cnt + 1
                tmp_indices(cnt) = ir
            endif
        enddo

        if (cnt == 0) then
            write(*, *) 'All weight values are zero - check the weight function!'
            error stop 101
        endif

        allocate(indices(cnt), source= tmp_indices(1:cnt))
        deallocate(tmp_indices)

    end subroutine select_finite_weight


    !> @brief Generate an array of random integers in the range [1:max_integer]
    ! TODO(Alex) THIS CAN give the same random integer more than once, which is problematic 
    !
    ! Updates to the standard: https://cyber.dabamos.de/programming/modernfortran/random-numbers.html    
    ! Some points on sampling non-uniform distributions: https://masuday.github.io/fortran_tutorial/random.html                                       
    subroutine generate_random_integers(max_integer, random_integers, seed)
        integer, intent(in)           :: max_integer
        integer, intent(out)          :: random_integers(:)
        integer, intent(in), optional :: seed(:)

        integer :: i
        real(real64) :: random_num 

        if (present(seed)) then
            call random_seed(put=seed)
        else
            call random_seed()
        endif

        do i = 1, size(random_integers)
            ! Generate random number between 0 and 1                                                                       
            call random_number(random_num)
            random_integers(i) = int(random_num * max_integer + 1, kind=int32)
        end do

    end subroutine generate_random_integers


    !> @brief Generate an array of random integers in the range [1:max_integer]
    !! Transcribed from https://gist.github.com/Pseudomanifold/1ed2b3b5b0a1389bdad97b2fdf5af47e
    subroutine reservoir_sampling(n, selected, seed)
        integer, intent(in) :: n  !< Max integer 
        integer, intent(out) :: selected(:)
        integer, intent(in), optional :: seed(:)

        integer :: i, j, cnt, m
        real(real64) :: random_num
      
        m = size(selected)

        if (present(seed)) then
            call random_seed(put=seed)
        else
            call random_seed()
        endif

        cnt = 0
        do i = 1, n
           call random_number(random_num)
          if (real(n - i + 1) * random_num < m - cnt) then
            cnt = cnt + 1
            selected(cnt) = i
          end if
          if (cnt == m) exit
        end do
      
        ! check for duplicates
        ! do i = 1, m
        !    i_element = selected(i)
        !    do j = i + 1, m
        !       if (i_element == selected(j)) then
        !          write(*, *) 'Elements', i, 'and', j, 'are the same:', i_element, selected(j)
        !       endif
        !    enddo
        ! enddo
      
    end subroutine reservoir_sampling


    !> @brief Initialise centroids according to a weight function.
    !!
    !! Simple implementation:
    !! * Creates an index array of grid indices with finite weights
    !! * Generates a set of random indices in the range [1: n_finite_weights)] 
    !! * Samples initial centroids randomly from grid points with finite weights
    subroutine choose_initial_centroids_simple(grid, weight, centroids, fixed_seed)
        real(real64), intent(in)  :: grid(:, :)         !< (ndim, np)
        real(real64), intent(in)  :: weight(:)          !< (np)
        real(real64), intent(out) :: centroids(:, :)    !< (ndim, n_centroid)
        integer,      intent(in), optional :: fixed_seed(:)   !< Optional seed to fix rand nums generated

        integer, allocatable :: finite_weight_indices(:), random_integers(:)
        integer :: i, ic, ir, irand, n_centroid

        ! Avoid initialising centroids at points with no weight, as this breaks kmeans
        ! Points with finite weights
        n_centroid = size(centroids, 2)
        call select_finite_weight(weight, finite_weight_indices)

        ! Randomly sample from points with finite weight
        allocate(random_integers(n_centroid))
        !call generate_random_integers(size(finite_weight_indices), random_integers, fixed_seed)
        call reservoir_sampling(size(finite_weight_indices), random_integers, fixed_seed)

        do ic = 1, n_centroid
            irand = random_integers(ic)
            ir = finite_weight_indices(irand)
            centroids(:, ic) = grid(:, ir)
            ! write(103, *) ic, ir, centroids(:, ic), weight(ir)
        enddo

        deallocate(finite_weight_indices)
        deallocate(random_integers)

    end subroutine choose_initial_centroids_simple

end module sampling_m
