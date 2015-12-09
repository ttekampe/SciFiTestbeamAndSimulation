cmake_minimum_required(VERSION 2.8)

project(SciFiTestbeamAndSimulation)
set(PROJECT_VERSION 0.0.0.0)
set(PROJECT_BRIEF "\"No description\"")



# Enable/disable testing
option(TESTS_ENABLED "Enable automatic tests" OFF)
# Enable/disable test coverage
#option(COVERAGE_ENABLED "Enable test coverage" ON)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")

include_directories(${CMAKE_SOURCE_DIR}/src)

#if(COVERAGE_ENABLED)
#    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -fPIC -O0 -ggdb -fprofile-arcs -ftest-coverage")
#else()
#    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -O0 -ggdb")
#endif()

#commented out -Wall
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -Wall")
#set(CMAKE_CXX_FLAGS_RELEASE "-g")
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(FATAL_ERROR "Compiler does not support C++11")
endif()


set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, defaulting to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

# Boost, a collection of useful C++ utilities and idioms
find_package(Boost COMPONENTS
    program_options
    filesystem
    thread
    system
    REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
set(LIBRARIES ${LIBRARIES} ${Boost_LIBRARIES})

# ROOT, a library for HEP data analysis (contains RooFit)
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIR})
set(LIBRARIES ${LIBRARIES} ${ROOT_LIBRARIES})
set(LIBRARIES ${LIBRARIES} ${ROOFIT_LIBRARIES})



include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})



include_directories(${ROOT_INCLUDE_DIR})


set(ALL_LIBRARIES ${ROOT_LIBRARIES} ${ROOFIT_LIBRARIES} ${ADDITIONAL_LIBRARIES} "-lTreePlayer")

add_subdirectory(src)
add_subdirectory(main)
#add_subdirectory(energyCorrection)

#MACRO(SUBDIRLIST result curdir)
#  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
#  SET(dirlist "")
#  FOREACH(child ${children})
#    IF(IS_DIRECTORY ${curdir}/${child})
#        LIST(APPEND dirlist ${child})
#    ENDIF()
#  ENDFOREACH()
#  SET(${result} ${dirlist})
#ENDMACRO()
#
#SUBDIRLIST(SUBDIRS ${MY_CURRENT_DIR})
#
#file(WRITE .clang_complete "")
#
#FOREACH(subdir ${SUBDIRS})
#    file(APPEND .clang_complete "-I"${subdir})
#    ADD_SUBDIRECTORY(${subdir})
#ENDFOREACH()
