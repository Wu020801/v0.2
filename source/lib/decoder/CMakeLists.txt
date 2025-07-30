message("source/lib/decoder")
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
					${DEPENDENCY_DIRECTORY}/hpm-HPM-15.0/app)

add_library(${MYNAME} ${SRC})

set(DRACO_LIB_PATH $<TARGET_FILE:draco>)
# add HPM libs
find_library(COMMON_LIB_PATH NAMES CommonLib)
find_library(DECODER_LIB_PATH NAMES DecoderLib)
link_directories(${DRACO_LIB_PATH} ${COMMON_LIB_PATH} ${DECODER_LIB_PATH})

target_link_libraries(${MYNAME} PUBLIC ${DRACO_LIB} common CommonLib DecoderLib)
set_target_properties(${MYNAME} PROPERTIES LINKER_LANGUAGE CXX)



