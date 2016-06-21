project (hermelin)
cmake_minimum_required(VERSION 3.0)

set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(CMAKE_BUILD_TYPE Release)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/.." "." "..")
set(CMAKE_BINARY_DIR "${CMAKE_SOURCE_DIR}")


add_custom_command(OUTPUT script.cc tokens.h
	COMMAND bison --defines=${PROJECT_SOURCE_DIR}/tokens.h --output-file=${PROJECT_SOURCE_DIR}/script.cc ${PROJECT_SOURCE_DIR}/lang/script.yy
)


find_package(Boost REQUIRED COMPONENTS iostreams system filesystem REQUIRED)
find_package(Threads REQUIRED)
find_package(Minuit2 REQUIRED)


add_definitions(-std=c++11)
add_definitions(-march=native -s)

add_definitions(-DNO_LAPACK -DNO_QT)
add_definitions(${Boost_CXX_FLAGS})
add_definitions(-DTLIBS_INC_HDR_IMPLS)

add_definitions(-DINSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")



include_directories("${PROJECT_SOURCE_DIR}" "." "..")
include_directories("${Boost_INCLUDE_DIRS}/..")
include_directories("${Minuit2_INCLUDE_DIRS}")

#include_directories("/usr/include/root")
#link_directories("/usr/lib64/root")


file(GLOB lang "lang/*.cpp" "lang/*.h" "${PROJECT_SOURCE_DIR}/script.cc" "${PROJECT_SOURCE_DIR}/tokens.h")
file(GLOB runtime "runtime/*.cpp" "runtime/*.h")
file(GLOB tlibs
	"../math/linalg2.cpp" "../math/rand.cpp" "../math/fourier.cpp"
	"../string/spec_char.cpp"
	"../gfx/gnuplot.cpp"
	"../file/file.cpp"
	"../log/log.cpp" "../helper/flags.cpp"
)


add_executable(hermelin 
	${lang} ${runtime} ${tlibs} "${PROJECT_SOURCE_DIR}/script.cc" "${PROJECT_SOURCE_DIR}/tokens.h"
)

target_link_libraries(hermelin 
	${Boost_LIBRARIES} Threads::Threads ${Minuit2_LIBRARIES}
)


add_custom_command(TARGET hermelin POST_BUILD 
	COMMAND strip hermelin 
	MAIN_DEPENDENCY hermelin
)



install(TARGETS hermelin DESTINATION bin)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/apps DESTINATION share/hermelin)
