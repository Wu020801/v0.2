include_directories(${CMAKE_SOURCE_DIR}/dependencies
                    ${CMAKE_SOURCE_DIR}/source/lib/encoder
                    ${CMAKE_SOURCE_DIR}/source/lib/decoder
                    ${CMAKE_SOURCE_DIR}/source/lib/common)



add_executable(encode encode.cpp)
if (WIN32)
    target_link_libraries(encode encoder common UVAtlas DirectX-Headers DirectXMath DirectXMesh)
elseif (UNIX)
    target_link_libraries(encode encoder common UVAtlas DirectX-Headers DirectXMath DirectXMesh stdc++fs)
endif()

add_executable(decode decode.cpp)
if (WIN32)
    target_link_libraries(decode decoder common)
elseif (UNIX)
    target_link_libraries(decode decoder common stdc++fs)
endif()
