cmake_minimum_required(VERSION 3.13)
cmake_policy(SET CMP0076 NEW)

# set the project name
project(SpaSM VERSION 1.3 
	LANGUAGES C CXX
	HOMEPAGE_URL https://github.com/cbouilla/spasm)

# if(NOT CMAKE_BUILD_TYPE)
#   message(VERBOSE "Setting build type to 'Debug' as none was specified.")
#   set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build." FORCE)
# endif()
# 
# message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
# message(STATUS "Configuration type: ${CMAKE_CONFIGURATION_TYPES}")

if(CMAKE_BUILD_TYPE)
    message(WARNING "Non-standard build type: ${CMAKE_BUILD_TYPE}")
    string(TOLOWER ${CMAKE_BUILD_TYPE} build_type)
endif()

# specify the C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(CTest)
find_package(OpenMP REQUIRED)

# pkg-config then libm4ri
find_package(PkgConfig REQUIRED)
pkg_check_modules(GIVARO REQUIRED IMPORTED_TARGET givaro)
pkg_check_modules(FFLAS_FFPACK REQUIRED IMPORTED_TARGET fflas-ffpack)

# set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/config;${CMAKE_MODULE_PATH}")
# find_package(PAPI)

if (build_type STREQUAL cov)
    add_compile_options(-Wall -O0 --coverage -Wno-pragmas -fopenmp)
    add_link_options(--coverage)
    message(WARNING "special coverage options")
else()
    # normal
    add_compile_options(-Wall -Wextra -O3 -g)
endif()

#add_link_options(-z noexecstack)

########################################################################

add_subdirectory(src)
add_subdirectory(tools)

if (BUILD_TESTING)
  add_subdirectory(tests)
endif()
