cmake_minimum_required(VERSION 3.10)

project(mmf LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 90)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

add_compile_options(-fPIC)

set(SRC_FILES 
  "./src/module_parallel.f90"
  "./src/interp_lib.f90"
  "./src/module_rootdepth.f90"
  "./src/module_wtable.f90"
  "./src/module_nrtype.f90"
  "./src/module_initial.f90"
  "./src/module_io.f90"
  "./src/module_forcings.f90"
)
# message(STATUS "SRC_FILES: ${SRC_FILES}")

find_package(NetCDF REQUIRED)
include_directories(
  ${NetCDF_F90_INCLUDEDIR} 
  # ${CMAKE_BINARY_DIR}/share
  # ${CMAKE_BINARY_DIR}/mksrfdata
)

# # add_subdirectory(../share ${CMAKE_BINARY_DIR}/share)
# add_subdirectory(../mksrfdata ${CMAKE_BINARY_DIR}/mksrfdata)

# # add_library(${PROJECT_NAME} STATIC ${SRC_FILES})
add_library(${PROJECT_NAME} SHARED ${SRC_FILES})

# # add_executable(${PROJECT_NAME} ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} PRIVATE netcdff)
