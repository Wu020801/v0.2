cmake_minimum_required(VERSION 3.12)
project(MCEM CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (UNIX)
    add_compile_options(-Wno-stringop-overflow)
    add_compile_options(-Wno-error)
endif()

message("avs_mesh_compression")
set (ROOT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
# Set dependency path
SET( DEPENDENCY_DIRECTORY  "${CMAKE_CURRENT_SOURCE_DIR}/dependencies" CACHE STRING "Set sepcific folder to clone dependencies")  
get_filename_component( DEPENDENCY_DIRECTORY ${DEPENDENCY_DIRECTORY} REALPATH)
if( NOT EXISTS ${DEPENDENCY_DIRECTORY} )
  file(MAKE_DIRECTORY ${DEPENDENCY_DIRECTORY})
endif()
MESSAGE("DEPENDENCY_DIRECTORY = ${DEPENDENCY_DIRECTORY}" )

add_subdirectory(dependencies/hpm-HPM-15.0)
include( dependencies/cmake/directx.cmake )
include(dependencies/cmake/draco.cmake)
include(dependencies/cmake/uvatlas.cmake)

# Wrapper and core codec libraries
add_subdirectory(source/lib/encoder)
add_subdirectory(source/lib/decoder)
add_subdirectory(source/lib/common)

# Applications
add_subdirectory(source/app)