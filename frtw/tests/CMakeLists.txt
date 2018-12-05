# Make sure the compiler can find include files from our nttw library.
include_directories (${FRTW_SOURCE_DIR}/include ${NTTW_SOURCE_DIR}/include)

#Additional Platform specific libs
if(CMAKE_COMPILER_IS_GNUCC)
	set (PLATFORM_LIBS rt)
endif(CMAKE_COMPILER_IS_GNUCC)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Make sure the linker can find the nttw library once it is built.
link_directories (${LIBRARY_OUTPUT_PATH})

# Add executable called "array_test" that is built from the source files
# "array_test.c". The extensions are automatically found.
add_executable (fft_test fft_test/fft_test.c)
# Link the executable to the NTTW library.
target_link_libraries (fft_test nttw frtw ${PLATFORM_LIBS}) 

#do the same for the remaining tests
add_executable (gauss_noise_test gauss_noise_test/gauss_noise_test.c)
target_link_libraries (gauss_noise_test frtw nttw ${PLATFORM_LIBS})
