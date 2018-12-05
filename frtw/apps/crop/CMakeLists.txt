# Make sure the compiler can find include files from our nttw library.
include_directories (${FRTW_SOURCE_DIR}/include ${NTTW_SOURCE_DIR}/include)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Make sure the linker can find the nttw library once it is built.
link_directories (${LIBRARY_OUTPUT_PATH})

# Add executable called "ntt" that is built from the source files
# "ntt.c". The extensions are automatically found.
add_executable (crop crop.c)

#Additional Platform specific libs
if(CMAKE_COMPILER_IS_GNUCC)
	set (PLATFORM_LIBS rt)
endif(CMAKE_COMPILER_IS_GNUCC)

# Link the executable to the NTTW library.
target_link_libraries (crop nttw frtw ${PLATFORM_LIBS}) 

# Installation of the program
install(TARGETS crop
        DESTINATION "bin"
        PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)