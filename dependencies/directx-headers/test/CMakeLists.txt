# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.

list(APPEND dxlibs "")
# Check if in WSL and if has DirectX driver and runtime
if(EXISTS "/usr/lib/wsl/lib/")
  find_library(libd3d12 d3d12 HINTS /usr/lib/wsl/lib)
  find_library(libdxcore dxcore HINTS /usr/lib/wsl/lib)
  list(APPEND dxlibs ${libd3d12} ${libdxcore})
else()
# Fallback to default: let CMake look for libs
  list(APPEND dxlibs d3d12 dxcore)
endif()

project(DirectX-Headers-Test
  DESCRIPTION "DirectX-Header test"
  HOMEPAGE_URL "https://github.com/microsoft/DirectX-Headers/"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(TEST_EXES DirectX-Headers-Test DirectX-Headers-Check-Feature-Support-Test)

add_executable(DirectX-Headers-Test test.cpp)

add_executable(DirectX-Headers-Check-Feature-Support-Test feature_check_test.cpp)

foreach(t IN LISTS TEST_EXES)
  target_link_libraries(${t} DirectX-Headers DirectX-Guids ${dxlibs})
endforeach()

if ( CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
    target_compile_options(DirectX-Headers-Check-Feature-Support-Test PRIVATE -Wno-unused-variable)
endif()

if(WIN32)
    foreach(t IN LISTS TEST_EXES)
      target_compile_definitions(${t} PRIVATE _UNICODE UNICODE _WIN32_WINNT=0x0A00)
    endforeach()

    if(WINDOWS_STORE)
      foreach(t IN LISTS TEST_EXES)
        target_compile_definitions(${t} PRIVATE WINAPI_FAMILY=WINAPI_FAMILY_APP)
      endforeach()
    endif()
endif()
