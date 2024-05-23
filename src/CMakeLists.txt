target_sources(isdf_fortran PRIVATE
       isdf_serial_m.F90
       maths_m.f90
       face_splitting_m.F90
       parse.f90
       sampling_m.f90
       interpolation_points_m.F90
       )

target_sources(fortuno_wrapper PRIVATE
       fortuno_interface_m.F90
       )
