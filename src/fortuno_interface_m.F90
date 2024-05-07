!> @brief Expose Fortuno routines
module fortuno_interface_m
#ifdef USE_MPI    
    use fortuno_mpi,    only : execute_cmd_app => execute_mpi_cmd_app, &
        test => mpi_case_item, &
        check => mpi_check, &
        is_equal, &
        test_item, &
        global_comm, &
        this_rank, &
        as_char
    implicit none
    ! Can only do this within intrinsics
    !! integer :: comm_world = global_comm()

#else
    use fortuno_serial, only : execute_cmd_app => execute_serial_cmd_app, &
                               test => serial_case_item, &
                               check => serial_check, &
                               is_equal, &
                               test_item, &
                               as_char
    implicit none
    ! Note, probably no need to overload this
    integer, parameter :: comm_world = 100
#endif  

end module fortuno_interface_m
