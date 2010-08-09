# The name the project is "NTTW". CMakeLists files in this project can
# refer to the root source directory of the project as ${NTTW_SOURCE_DIR} and
# to the root binary directory of the project as ${NTTW_BINARY_DIR}.
cmake_minimum_required (VERSION 2.6)
project (NTTW C)

#Set output paths
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

#Add compiler flags
#GCC
if(CMAKE_COMPILER_IS_GNUCC)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native -mfpmath=sse -msse2 -funroll-all-loops -fwrapv")
endif(CMAKE_COMPILER_IS_GNUCC)
#MSVC
if(MSVC)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Ox /Ot /Oi /fp:fast /arch:SSE2 /MP /GS-")
endif(MSVC)

# Recurse into the "apps" subdirectory. This does not actually
# cause another cmake executable to run.
add_subdirectory (src)
add_subdirectory (apps)
