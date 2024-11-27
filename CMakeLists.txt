cmake_minimum_required(VERSION 3.1)

project(multidim LANGUAGES Fortran)
enable_language(Fortran)

find_package(LAPACK)
if(LAPACK_FOUND)
	set(IS_LAPACK .true.)
	configure_file(${PROJECT_SOURCE_DIR}/src/libraries.f95.in 
		${PROJECT_SOURCE_DIR}/src/libraries.f95)
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -llapack")
else()
	set(IS_LAPACK .false.)
endif()

set(FVERSION "-std=f2008")
set(FLAG "-O2 -fcheck=all -g -fimplicit-none --pedantic -Wall -Wextra")
set(PREPROCESS "-cpp") # enable preprocessor
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FVERSION} ${FLAG} ${PREPROCESS}")

# Source code
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)

include_directories(src)
link_directories(src)

add_subdirectory(src)


