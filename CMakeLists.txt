cmake_minimum_required(VERSION 3.20)

project(ISDF_Fortran_Prototyping LANGUAGES Fortran)

# Build directory
set(CMAKE_Fortran_BIN_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_Fortran_LIB_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

# Compiler flags for GCC only
set(GCC_BASE 
   -fopenmp   # OMP always enabled
   -cpp       # Use c-preprocessor
   )

set(GCC_DEBUG
     -O0                                          # 0 optimisation. Probably implicit
     -fbacktrace                                  # symbolic stack traceback
     -ffpe-trap=invalid,zero,overflow,underflow   # control over floating-point exception
     -finit-real=nan                              # All real scalars are initialised to NaN
     -fcheck=all                                  # Enable all run-time test of -fcheck: array-temps, bits, bounds, do, mem, pointer, recursion
     -Wall                                        # Enable all warnings
    )

set(GCC_RELEASE
    -O3
   )

# Convert list to string
list(JOIN GCC_BASE " " GCC_BASE)
list(JOIN GCC_DEBUG " " GCC_DEBUG)
list(JOIN GCC_RELEASE " " GCC_RELEASE)

# Set flags globally - Note, this is a bad strategy as I DO NOT want to get compiler warnings from external libs, etc
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${GCC_BASE}") 
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${GCC_BASE} ${GCC_DEBUG}")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${GCC_BASE} ${GCC_RELEASE}")

# Interface to control scope of compiler options
add_library(common_options INTERFACE)

# External dependencies
# ---------------------------------
include(FetchContent)

find_package(BLAS REQUIRED)
find_package(OpenMP REQUIRED)

# Fortuno unit testing framework
FetchContent_Declare(
  Fortuno
  GIT_REPOSITORY "https://github.com/fortuno-repos/fortuno"
  GIT_TAG "main"
)
FetchContent_MakeAvailable(Fortuno)
set(fortuno_libs Fortuno::fortuno_serial)

# MPI-dependent settings
option(MPI "Use MPI" OFF)
if(MPI)
   message("-- Using MPI")

   # Find MPI library and headers
   find_package(MPI REQUIRED)
   add_definitions(${MPI_Fortran_COMPILE_FLAGS})
   include_directories(${MPI_Fortran_INCLUDE_PATH})

   # Fetch and assign fortuno MPI lib
   FetchContent_Declare(
     FortunoMPI
     GIT_REPOSITORY "https://github.com/fortuno-repos/fortuno-mpi"
     GIT_TAG "main"
   )
   FetchContent_MakeAvailable(FortunoMPI)
   list(APPEND fortuno_libs Fortuno::fortuno_mpi)

   # Set preprocessor variable 
   target_compile_definitions(common_options INTERFACE USE_MPI)

else()
   set(MPI_Fortran_LIBRARIES "")
   message("-- Not using MPI")
endif()

# This library
# -----------------------------
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

add_library(isdf_fortran "")

set_target_properties(isdf_fortran
        PROPERTIES
        VERSION 0.0.1
        SOVERSION 0.0.1
        LIBRARY_OUTPUT_NAME isdff90
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        )

# Linked libraries
target_link_libraries(isdf_fortran PRIVATE OpenMP::OpenMP_Fortran ${MPI_Fortran_LIBRARIES} BLAS::BLAS common_options)

# Source code
target_sources(isdf_fortran PUBLIC
       src/isdf_m.F90
       src/maths_m.f90
       )

# Fortuno wrapper layer
# -----------------------------
add_library(fortuno_wrapper "")
set_target_properties(fortuno_wrapper
        PROPERTIES
        VERSION 0.0.1
        SOVERSION 0.0.1
        LIBRARY_OUTPUT_NAME fortunow
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        )
target_link_libraries(fortuno_wrapper PUBLIC ${fortuno_libs} ${MPI_Fortran_LIBRARIES})
target_sources(fortuno_wrapper PRIVATE
       src/fortuno_interface_m.F90
)

enable_testing()

#  Library unit tests
# -------------------------------

add_executable(test_isdf tests/test_isdf_m.F90)
target_link_libraries(test_isdf PRIVATE isdf_fortran fortuno_wrapper common_options) 
add_test(NAME test_stub COMMAND test_isdf)
