# The name the project is "NTTW". CMakeLists files in this project can
# refer to the root source directory of the project as ${NTTW_SOURCE_DIR} and
# to the root binary directory of the project as ${NTTW_BINARY_DIR}.
cmake_minimum_required (VERSION 2.6)
project (NTTW C)

#Set output paths
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# 64-bit Build?
IF(CMAKE_SIZEOF_VOID_P EQUAL 8) #64-bit environment
  set (NTTW_64 ON) #True
  message("Detected 64-bit system. Using appropriate 64-bit integers.")
ELSEIF(CMAKE_SIZEOF_VOID_P EQUAL 8) #64-bit environment
  set (NTTW_64 OFF) #False
ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 8) #64-bit environment

#Add compiler flags
#GCC
if(CMAKE_COMPILER_IS_GNUCC)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native -fwrapv")
endif(CMAKE_COMPILER_IS_GNUCC)
#MSVC
if(MSVC)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Ox /Ot /Oi /MP /GS-")
endif(MSVC)

if(BUILD_SHARED_LIBS)
message("Building shared libraries")
add_definitions (-DNTTW_DLL -DNTTW_MAKEDLL)
endif(BUILD_SHARED_LIBS)
if (NTTW_64)
	add_definitions (-DNTTW_64)
endif (NTTW_64)

#Build options
OPTION(BUILD_APPLICATIONS "Build the applications for NTTW" ON)
OPTION(BUILD_BENCH_APPLICATIONS "Build the benchmark applications for NTTW" OFF)
OPTION(BUILD_TESTING "Build the test applications for NTTW" OFF)
OPTION(BUILD_ITK_SUPPORT "Build the test applications for NTTW that use the Insight Toolkit" OFF)

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
set(CPACK_PACKAGE_NAME "libnttw")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/Copyright.txt")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "41")
set(CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")

#Debian
if(UNIX)
  set(CPACK_GENERATOR "DEB")
  set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Shekhar Chandra <shakes76@users.sourceforge.net>")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "NTTW library is for Number Theoretic Transforms (NTTs).")
  set(CPACK_PACKAGE_DESCRIPTION "NTTs are pure digital transforms computed on integer fields. This allows convolutions to be done exactly and essentially without overflow.")
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6 (>= 2.7-18)")
  set(CPACK_DEBIAN_PACKAGE_SECTION "Science")
endif(UNIX)

include(CPack)