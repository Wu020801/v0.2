message("source/lib/common")
# message(${CMAKE_CURRENT_SOURCE_DIR})

get_filename_component(MYNAME ${CMAKE_CURRENT_LIST_DIR} NAME)
string(REPLACE " " "_" MYNAME ${MYNAME})
project( ${MYNAME} C CXX )

file(GLOB SRC *.h *.cpp *.c)

include_directories(${CMAKE_BINARY_DIR} 
					${DEPENDENCY_DIRECTORY} )

add_library(${MYNAME} ${SRC})

set_target_properties(${MYNAME} PROPERTIES LINKER_LANGUAGE CXX)
