#[==============================================================================================[
#                            CMakeLists.txt for ISDF Project                                  #
#                                                                                               #
#]==============================================================================================]
cmake_minimum_required(VERSION 3.20)
project(ISDF_Fortran_Prototyping VERSION 1.0.0 LANGUAGES Fortran)

# Set sub-directories in the build folder
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

# Load and execute CMake code from provided modules
include(GNUInstallDirs)
include(FetchContent)

# Add custom modules in cmake/
# list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#[==============================================================================================[
#                        Compiler Flags and Options                                             #
# Note, flags "should" be specified in cmake presets rather in CMakeLists.txt                   #
#]==============================================================================================]

set(GCC_BASE 
   -fopenmp
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

# Interface to control scope of compiler flags and options
add_library(common_flags INTERFACE)

# Build-specific flags
target_compile_options(common_flags INTERFACE 
   $<$<AND:$<CONFIG:Debug>,$<Fortran_COMPILER_ID:GNU>>:${GCC_BASE} ${GCC_DEBUG}>
   $<$<AND:$<CONFIG:Release>,$<Fortran_COMPILER_ID:GNU>>:${GCC_BASE} ${GCC_RELEASE}>
   )

#[==============================================================================================[
#                                   External Dependencies                                      #
#]==============================================================================================]
find_package(BLAS REQUIRED)
find_package(OpenMP REQUIRED)
# TODO(Alex) This should be removed beyond personal use, and replaced with FetchContent
set(CMAKE_PREFIX_PATH "/Users/alexanderbuccheri/Codes/kmeans/serial_install")
find_package(kmeans REQUIRED)

# Fortuno unit testing framework
FetchContent_Declare(
  Fortuno
  GIT_REPOSITORY "https://github.com/fortuno-repos/fortuno"
  GIT_TAG "main"
)
FetchContent_MakeAvailable(Fortuno)
set(fortuno_libs Fortuno::fortuno_serial)

# MPI-dependent settings
option(MPI "Use MPI" ON)
if(MPI)
   find_package(MPI REQUIRED Fortran)

   # Fetch and assign fortuno MPI lib
   FetchContent_Declare(
     FortunoMPI
     GIT_REPOSITORY "https://github.com/fortuno-repos/fortuno-mpi"
     GIT_TAG "main"
   )
   FetchContent_MakeAvailable(FortunoMPI)
   list(APPEND fortuno_libs Fortuno::fortuno_mpi)

   # Set preprocessor variable 
   target_compile_definitions(common_flags INTERFACE USE_MPI)

else()
   message("-- Not using MPI")
endif()

#[==============================================================================================[
#                                   Core Library                                                #
#]==============================================================================================]
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# Main library
add_library(isdf_fortran)
add_library(fortuno_wrapper)

set_target_properties(isdf_fortran
        PROPERTIES
        VERSION "${PROJECT_VERSION}"
        SOVERSION "${PROJECT_VERSION_MAJOR}"
        LIBRARY_OUTPUT_NAME isdff90
        LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
        )

# Linked libraries
target_link_libraries(isdf_fortran PRIVATE 
   OpenMP::OpenMP_Fortran 
   $<$<BOOL:${MPI}>:MPI::MPI_Fortran> 
   common_flags 
   BLAS::BLAS 
   )

# Location of source - added after defining library variables
add_subdirectory(src)

set_target_properties(fortuno_wrapper PROPERTIES
        VERSION "${PROJECT_VERSION}"
        SOVERSION "${PROJECT_VERSION_MAJOR}"
        LIBRARY_OUTPUT_NAME fortunow
        LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
        )

target_link_libraries(fortuno_wrapper PUBLIC 
   ${fortuno_libs} 
   $<$<BOOL:${MPI}>:MPI::MPI_Fortran>
   )

#[==============================================================================================[
#                                   Testing                                                     #
#]==============================================================================================]
enable_testing()

# Unit tests
add_executable(test_face_splitting tests/test_face_splitting_m.F90)
target_link_libraries(test_face_splitting PRIVATE isdf_fortran fortuno_wrapper common_flags) 
add_test(NAME test_face_splitting_product_two_funcs COMMAND test_face_splitting "Face-splitting of two products")
add_test(NAME test_face_splitting_product_one_func COMMAND test_face_splitting "Face-splitting of the product of one matrix")

add_executable(test_maths tests/test_maths_m.F90)
target_link_libraries(test_maths PRIVATE isdf_fortran fortuno_wrapper common_flags) 
add_test(NAME test_gram_schmidt COMMAND test_maths "Gram-Schmidt orthogonalisation")


# Application tests
add_executable(run_isdf_serial regression_tests/run_isdf_serial.F90)
target_link_libraries(run_isdf_serial PRIVATE 
   isdf_fortran kmeans::lib_kmeans 
   $<$<BOOL:${MPI}>:MPI::MPI_Fortran>  
   common_flags) 
