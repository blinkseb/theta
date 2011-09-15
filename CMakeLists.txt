cmake_minimum_required(VERSION 2.6)
project(theta)
set(CMAKE_VERBOSE_MAKEFILE, 1)

option(coverage "build for coverage tests" OFF)
option(profile "build for profiling" OFF)
option(optiontest "build tests" OFF)
option(release "build release version" ON)
option(crlibm "use the log function from crlibm" ON)
option(generic "use generic implementation of funtions; not optimized, platform-specific ones" OFF)
option(psql "use postgresql" OFF)

if(coverage)
   message(STATUS "building for coverage test")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage -O0")
elseif(profile)
   message(STATUS "building for profiling")
   set(CMAKE_BUILD_TYPE Release) # to enable optimizations
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
   #set(CMAKE_LD_FLAGS "${CMAKE_LD_FLAGS} -pg")
else()
  if(release)
     set(CMAKE_BUILD_TYPE Release)
     message(STATUS "building release version")
  else()
     message(STATUS "building release version +debug")
     set(CMAKE_BUILD_TYPE RelWithDebInfo)
   endif()
endif()

if(crlibm)
   message(STATUS "using log from crlibm")
   add_definitions(-DUSE_CRLIBM)
endif()


if(NOT generic)
  enable_language(ASM-ATT)
  if(NOT "${CMAKE_SYSTEM_PROCESSOR}" MATCHES "x86_64" OR NOT CMAKE_ASM-ATT_COMPILER_WORKS)
      message(STATUS "non-optimized platform / assembler not working; enabling 'generic' option")
      set(generic ON)
  endif()
endif()
   

if(generic)
   message(STATUS "using generic functions")
   add_definitions(-DGENERIC_ARCH)
else()
   enable_language(ASM-ATT)
endif()

find_package(Boost COMPONENTS regex program_options filesystem system iostreams)

if(NOT Boost_FOUND)
   message(SEND_ERROR "boost not found")
endif()

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(test_dir ${CMAKE_SOURCE_DIR}/test)

add_subdirectory(src)
add_subdirectory(plugins)
add_subdirectory(libconfig)
add_subdirectory(bin)
add_subdirectory(root)

if(optiontest)
   message(STATUS "building tests")
   add_subdirectory(test)
endif()

