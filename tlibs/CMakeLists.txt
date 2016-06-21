#
# tlibs
# @author tweber
# @date 2012-2016
# @license GPLv2 or GPLv3
#

project(tlibs)
cmake_minimum_required(VERSION 3.0)

#set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(CMAKE_BUILD_TYPE Release)


list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}" ".")
set(CMAKE_BINARY_DIR "${CMAKE_SOURCE_DIR}")
#set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

find_package(Boost REQUIRED COMPONENTS iostreams system filesystem regex REQUIRED)
find_package(Threads REQUIRED)
find_package(Qt4 REQUIRED)
find_package(Freetype REQUIRED)
find_package(Minuit2 REQUIRED)

set(CMAKE_AUTOMOC OFF)
set(CMAKE_AUTOUIC OFF)


#message("Project source dir: ${PROJECT_SOURCE_DIR}")
#message("CMake source dir: ${CMAKE_SOURCE_DIR}")
#message("Binary dir: ${CMAKE_BINARY_DIR}")
message("Install path prefix: ${CMAKE_INSTALL_PREFIX}")


add_definitions(-DNDEBUG)
add_definitions(-DUSE_IOSTR -DUSE_BOOST_REX)
add_definitions(-DNO_LAPACK)
#add_definitions(-DUSE_FFTW)
add_definitions(-DINSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")
add_definitions(${Boost_CXX_FLAGS})


include_directories("${PROJECT_SOURCE_DIR}" ".")
include_directories("${Boost_INCLUDE_DIRS}/..")
include_directories("${QT_INCLUDE_DIR}")
include_directories("${FREETYPE_INCLUDE_DIRS}")
include_directories("${Minuit2_INCLUDE_DIRS}")


add_definitions(-std=c++11)
#add_definitions(-march=native -s)



# -----------------------------------------------------------------------------
# build
# -----------------------------------------------------------------------------
add_library(tlibs SHARED
	file/tmp.cpp  file/loadinstr.cpp  file/recent.cpp  file/x3d.cpp
	#fit/minuit.cpp
	log/log.cpp log/debug.cpp
	helper/flags.cpp
	math/fourier.cpp  math/linalg2.cpp  math/rand.cpp
	net/tcp.cpp
	string/spec_char.cpp  string/eval.cpp
	gfx/gl.cpp  gfx/gnuplot.cpp
	version.cpp
)


target_link_libraries(tlibs ${Boost_LIBRARIES} Threads::Threads
	Qt4::QtCore Qt4::QtGui #Qt4::QtOpenGL
	${FREETYPE_LIBRARY} ${Minuit2_LIBRARIES}
)


if(CMAKE_BUILD_TYPE STREQUAL "Release")
	add_custom_command(TARGET tlibs POST_BUILD
		COMMAND strip -v $<TARGET_FILE:tlibs>
		MAIN_DEPENDENCY tlibs
	)
endif()
# -----------------------------------------------------------------------------




# -----------------------------------------------------------------------------
# install
# -----------------------------------------------------------------------------
install(TARGETS tlibs  DESTINATION lib)

file(GLOB tlibs_rt  ${PROJECT_SOURCE_DIR}/*.h)
file(GLOB tlibs_file  ${PROJECT_SOURCE_DIR}/file/*.h)
file(GLOB tlibs_helper  ${PROJECT_SOURCE_DIR}/helper/*.h)
file(GLOB tlibs_math  ${PROJECT_SOURCE_DIR}/math/*.h ${PROJECT_SOURCE_DIR}/math/*.hpp)
file(GLOB tlibs_net  ${PROJECT_SOURCE_DIR}/net/*.h)
file(GLOB tlibs_string  ${PROJECT_SOURCE_DIR}/string/*.h)
file(GLOB tlibs_gfx  ${PROJECT_SOURCE_DIR}/gfx/*.h)
file(GLOB tlibs_time  ${PROJECT_SOURCE_DIR}/time/*.h)
file(GLOB tlibs_cl  ${PROJECT_SOURCE_DIR}/cl/*.h)
file(GLOB tlibs_log  ${PROJECT_SOURCE_DIR}/log/*.h)
file(GLOB tlibs_fit  ${PROJECT_SOURCE_DIR}/fit/*.h)

install(FILES ${tlibs_rt}  DESTINATION include/tlibs)
install(FILES ${tlibs_file}  DESTINATION include/tlibs/file)
install(FILES ${tlibs_helper}  DESTINATION include/tlibs/helper)
install(FILES ${tlibs_math}  DESTINATION include/tlibs/math)
install(FILES ${tlibs_net}  DESTINATION include/tlibs/net)
install(FILES ${tlibs_string}  DESTINATION include/tlibs/string)
install(FILES ${tlibs_gfx}  DESTINATION include/tlibs/gfx)
install(FILES ${tlibs_time}  DESTINATION include/tlibs/time)
install(FILES ${tlibs_cl}  DESTINATION include/tlibs/cl)
install(FILES ${tlibs_log}  DESTINATION include/tlibs/log)
install(FILES ${tlibs_fit}  DESTINATION include/tlibs/fit)

install(FILES ${PROJECT_SOURCE_DIR}/FindMinuit2.cmake DESTINATION share/tlibs)

install(FILES ${PROJECT_SOURCE_DIR}/COPYING.GPLv2 ${PROJECT_SOURCE_DIR}/COPYING.GPLv3 
	${PROJECT_SOURCE_DIR}/AUTHORS ${PROJECT_SOURCE_DIR}/LICENSES
        DESTINATION share/tlibs
)
# -----------------------------------------------------------------------------
