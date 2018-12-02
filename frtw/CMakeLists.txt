# The name the project is "FRTW". CMakeLists files in this project can
# refer to the root source directory of the project as ${FRTW_SOURCE_DIR} and
# to the root binary directory of the project as ${FRTW_BINARY_DIR}.
cmake_minimum_required (VERSION 2.6)
project (FRTW C)

#Set output paths
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

#Add compiler flags
add_definitions(-DNTTW_MAKEDLL)
#GCC
if(CMAKE_COMPILER_IS_GNUCC)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native -fwrapv")
endif(CMAKE_COMPILER_IS_GNUCC)
#MSVC
if(MSVC)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Ox /Ot /Oi /MP /GS-")
endif(MSVC)

#Need FFTW
find_path(FFTW_INCLUDES fftw3.h)
find_library (FFTW_LIBRARIES NAMES fftw3)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)

mark_as_advanced (FFTW_LIBRARIES FFTW_INCLUDES)

#Need NTTW
IF(NOT NTTW_DIR)
    IF(WIN32)
        IF(CMAKE_SIZEOF_VOID_P EQUAL 8 AND EXISTS $ENV{ProgramW6432}) #Win64 environment
            set(NTTW_DIR $ENV{ProgramW6432}/NTTW CACHE PATH "NTTW Base Directory")
        ELSE(CMAKE_SIZEOF_VOID_P EQUAL 8 AND EXISTS $ENV{ProgramW6432}) # else use Win32 environment
            set(NTTW_DIR $ENV{PROGRAMFILES}/NTTW CACHE PATH "NTTW Base Directory")
        ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 8 AND EXISTS $ENV{ProgramW6432})
    ENDIF(WIN32)
    IF(UNIX)
        set(NTTW_DIR /usr/local CACHE PATH "NTTW Base Directory")
    ENDIF(UNIX)
ENDIF(NOT NTTW_DIR)
#~ MARK_AS_ADVANCED(NTTW_DIR)
message("Using NTTW from ${NTTW_DIR}")
# Add Headers into List
file(GLOB NTTW_HEADERS ${NTTW_DIR}/include/*.h)
include_directories (${NTTW_DIR}/include ${FFTW_INCLUDES})
link_directories(${NTTW_DIR}/lib ${FFTW_LIBRARIES})

# 64-bit Build?
IF(CMAKE_SIZEOF_VOID_P EQUAL 8) #64-bit environment
  set (NTTW_64 ON) #True
  message("Detected 64-bit system. Using appropriate 64-bit integers.")
ELSEIF(CMAKE_SIZEOF_VOID_P EQUAL 8) #64-bit environment
  set (NTTW_64 OFF) #False
ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 8) #64-bit environment
if (NTTW_64)
	add_definitions (-DNTTW_64)
endif (NTTW_64)

#Build options
OPTION(BUILD_APPLICATIONS "Build the applications for FRTW" ON)
OPTION(BUILD_BENCH_APPLICATIONS "Build the benchmark applications for FRTW" OFF)
OPTION(BUILD_TESTING "Build the test applications for FRTW" OFF)

# Recurse into the "apps" subdirectory. This does not actually
# cause another cmake executable to run.
add_subdirectory (src)
IF(BUILD_APPLICATIONS)
  add_subdirectory (apps)
ENDIF(BUILD_APPLICATIONS)
IF(BUILD_BENCH_APPLICATIONS)
  add_subdirectory (bench)
ENDIF(BUILD_BENCH_APPLICATIONS)
IF(BUILD_TESTING)
  add_subdirectory (tests)
ENDIF(BUILD_TESTING)

#-----------------------
#Packaging
include(InstallRequiredSystemLibraries)

#All platforms
set(CPACK_PACKAGE_NAME "libfrtw")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/Copyright.txt")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "50")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")
#set(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")

#Debian
#if(UNIX)
set(CPACK_GENERATOR "DEB")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Shekhar Chandra <shakes76@users.sourceforge.net>")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "FRTW provides the Discrete Radon Transform and its variants.")
#set(CPACK_PACKAGE_DESCRIPTION "NTTs are pure digital transforms computed on integer fields. This allows convolutions to be done exactly and essentially without overflow.")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6 (>= 2.7-18), libnttw (>= 1.24)")
set(CPACK_DEBIAN_PACKAGE_SECTION "Science")
#endif(UNIX)

include(CPack)