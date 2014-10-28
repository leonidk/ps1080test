cmake_minimum_required(VERSION 2.6)
project(PS1080Test)
if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
endif()

# Use our modified FindSDL2* modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PS1080Test_SOURCE_DIR}/cmake")

# Look up SDL2 and add the include directory to our include path
set(ENV{SDL2} "C:/Users/nid/usr/SDL2")
set(ENV{EIGEN} "C:/Users/nid/usr/eigen")

find_package(SDL2 REQUIRED)
find_package(EIGEN REQUIRED)

include_directories(${SDL2_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR})

if(NOT OPENNI2_ROOT AND ("ON" STREQUAL "ON"))
    get_filename_component(OPENNI2_LIBRARY_HINT "OPENNI_LIBRARY-NOTFOUND" PATH)
endif(NOT OPENNI2_ROOT AND ("ON" STREQUAL "ON"))
set(OPENNI2_SUFFIX)
if(WIN32 AND CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(OPENNI2_SUFFIX 64)
endif(WIN32 AND CMAKE_SIZEOF_VOID_P EQUAL 8)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(PC_OPENNI2 libopenni2)
endif(PKG_CONFIG_FOUND)
find_path(OPENNI2_INCLUDE_DIRS OpenNI.h
    HINTS /usr/include/openni2 /usr/include/ni2
    PATHS "$ENV{OPENNI2_INCLUDE${OPENNI2_SUFFIX}}"
    PATH_SUFFIXES openni openni2 include Include)
find_library(OPENNI2_LIBRARY
    NAMES OpenNI2	# No suffix needed on Win64
    HINTS /usr/lib
    PATHS "$ENV{OPENNI2_LIB${OPENNI2_SUFFIX}}"
    PATH_SUFFIXES lib Lib Lib64)	
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenNI2 DEFAULT_MSG OPENNI2_LIBRARY OPENNI2_INCLUDE_DIRS)
if(OPENNI2_FOUND)
    get_filename_component(OPENNI_LIBRARY_PATH ${OPENNI2_LIBRARY} PATH)
    set(OPENNI2_LIBRARY_DIRS ${OPENNI2_LIBRARY_PATH})
    set(OPENNI2_LIBRARIES "${OPENNI2_LIBRARY}")
    set(OPENNI2_REDIST_DIR $ENV{OPENNI2_REDIST${OPENNI2_SUFFIX}})
endif(OPENNI2_FOUND)
include_directories(${OPENNI2_INCLUDE_DIRS})


# Look in the Lesson0 subdirectory to find its CMakeLists.txt so we can build the executable
add_subdirectory(src)
