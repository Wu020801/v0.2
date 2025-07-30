message("source/lib/encoder")
# message(${CMAKE_CURRENT_SOURCE_DIR})

get_filename_component(MYNAME ${CMAKE_CURRENT_LIST_DIR} NAME)
string(REPLACE " " "_" MYNAME ${MYNAME})
project( ${MYNAME} C CXX )

file(GLOB SRC *.h *.cpp *.c)
include_directories(${CMAKE_BINARY_DIR} 
					${ROOT_PATH}/source/lib/common 
					${DEPENDENCY_DIRECTORY} 
					${DEPENDENCY_DIRECTORY}/draco/src 
					${DEPENDENCY_DIRECTORY}/hpm-HPM-15.0/inc 
					${DEPENDENCY_DIRECTORY}/hpm-HPM-15.0/src 
					${DEPENDENCY_DIRECTORY}/hpm-HPM-15.0/app 
					${DEPENDENCY_DIRECTORY}/uvatlas/UVAtlas/inc
                    ${DEPENDENCY_DIRECTORY}/directx-math/Inc
                    ${DEPENDENCY_DIRECTORY}/directx-mesh/DirectXMesh
                    ${DEPENDENCY_DIRECTORY}/directx-headers/include)

add_library(${MYNAME} ${SRC})

set(DRACO_LIB_PATH $<TARGET_FILE:draco>)
# HPM libs
find_library(COMMON_LIB_PATH NAMES CommonLib)
find_library(ENCODER_LIB_PATH NAMES EncoderLib)
link_directories(${DRACO_LIB_PATH} ${COMMON_LIB_PATH} ${ENCODER_LIB_PATH})

target_link_libraries(${MYNAME} PUBLIC ${DRACO_LIB} common CommonLib EncoderLib)
set_target_properties(${MYNAME} PROPERTIES LINKER_LANGUAGE CXX)
