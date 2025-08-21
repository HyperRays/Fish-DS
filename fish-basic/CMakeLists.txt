cmake_minimum_required(VERSION 3.18)

project(FISH LANGUAGES Fortran)

# Require Fortran 2008
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# Default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# --- Source files ---
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)

file(GLOB_RECURSE FORTRAN_SOURCES
    "${SOURCE_DIR}/*.f90"
    "${SOURCE_DIR}/*.f"
)

# Configure paths for the Fortran code
set(DATA_PATH "${CMAKE_CURRENT_SOURCE_DIR}/data/")
set(LOGS_PATH "${CMAKE_CURRENT_SOURCE_DIR}/logs/")

configure_file(
    "${SOURCE_DIR}/mhd_parameter_module.f90"
    "${CMAKE_CURRENT_BINARY_DIR}/mhd_parameter_module_configured.f90"
    @ONLY
)

list(REMOVE_ITEM FORTRAN_SOURCES "${SOURCE_DIR}/mhd_parameter_module.f90")
list(APPEND FORTRAN_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/mhd_parameter_module_configured.f90")

# --- Target ---
add_executable(fish ${FORTRAN_SOURCES})

# Module output directory
set_target_properties(fish PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

target_include_directories(fish PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# --- Dependencies ---

# MPI
find_package(MPI REQUIRED COMPONENTS Fortran)
target_link_libraries(fish PRIVATE MPI::MPI_Fortran)

find_package(MKL REQUIRED)
target_link_libraries(fish PRIVATE MKL::MKL pthread)

# --- Compiler flags (scoped) ---
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    target_compile_options(fish PRIVATE
        -heap-arrays
        -qopenmp
        -mcmodel=medium
        -double-size 64
        -fpp
    )
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        target_compile_options(fish PRIVATE -O2 -shared-intel -qmkl)
    elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(fish PRIVATE 
            -g                    # Debug symbols
            -debug                # Debug info
            -check all            # All runtime checks
            -traceback           # Stack traceback on error
            -warn all            # All warnings
            -O0                  # No optimization
        )
    endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_compile_options(fish PRIVATE
        -fopenmp
        -fdefault-real-8
        -fdefault-double-8
    )
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        target_compile_options(fish PRIVATE -O2)
    elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(fish PRIVATE -g -fcheck=all -fbacktrace -Wall)
    endif()
endif()

# --- Output directories (runtime) ---
file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/logs)

# --- Info ---
message(STATUS "FISH Configuration:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  Fortran compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
message(STATUS "  Module directory: ${CMAKE_Fortran_MODULE_DIRECTORY}")
message(STATUS "  Output directory: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

# --- Convenience targets ---
add_custom_target(clean-data
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_SOURCE_DIR}/data
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/data
    COMMENT "Cleaning data directory"
)

add_custom_target(clean-logs
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_SOURCE_DIR}/logs
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/logs
    COMMENT "Cleaning logs directory"
)
