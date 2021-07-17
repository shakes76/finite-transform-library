# The name the project is "FTL". It consists of sub-projects NTTW and FRTW.
cmake_minimum_required (VERSION 2.6)
project (FTL C)

set(COMBINED_BUILD ON)

#-----------------------------------------------------------------------------
# Enable/Disable Shared Libraries
OPTION(BUILD_SHARED_LIBS "Build FTL with shared libraries." ON)
SET(FTL_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})

#Set output paths
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

subdirs(nttw)

set(NTTW_DIR ${PROJECT_SOURCE_DIR}/nttw CACHE PATH "NTTW Base Directory")
subdirs(frtw)