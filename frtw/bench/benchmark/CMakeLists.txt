# Make sure the compiler can find include files from our nttw library.
include_directories (${NTTW_SOURCE_DIR}/include)
include_directories (${FRTW_SOURCE_DIR}/include)

project(benchmark C)

#Additional Platform specific libs
if(CMAKE_COMPILER_IS_GNUCC)
	set (PLATFORM_LIBS rt)
endif(CMAKE_COMPILER_IS_GNUCC)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Make sure the linker can find the nttw library once it is built.
link_directories (${LIBRARY_OUTPUT_PATH})

# Add executable called "benchmark" that is built from the source files
# "benchmark.c". The extensions are automatically found.
add_executable (benchmark benchmark.c)
# Link the executable to the NTTW library.
target_link_libraries (benchmark frtw nttw ${PLATFORM_LIBS}) 

