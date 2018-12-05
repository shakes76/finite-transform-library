# Make sure the compiler can find include files from our nttw library.
include_directories (${NTTW_SOURCE_DIR}/include)

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
add_executable (array_test array_test/array_test.c)
# Link the executable to the NTTW library.
target_link_libraries (array_test nttw ${PLATFORM_LIBS}) 

#do the same for the remaining tests
add_executable (timing_test timing_test/timing_test.c)
target_link_libraries (timing_test nttw ${PLATFORM_LIBS})

add_executable (pgm_test pgm_test/pgm_test.c)
target_link_libraries (pgm_test nttw ${PLATFORM_LIBS}) 

add_executable (primitiveRoot_test primitiveRoot_test/primitiveRoot_test.c)
target_link_libraries (primitiveRoot_test nttw ${PLATFORM_LIBS}) 

add_executable (primeDash_test primeDash_test/primeDash_test.c)
target_link_libraries (primeDash_test nttw ${PLATFORM_LIBS}) 

#~ add_executable (ntt_test ntt_test/ntt_test.c)
#~ target_link_libraries (ntt_test nttw ${PLATFORM_LIBS})

add_executable (fntt_test fntt_test/fntt_test.c)
target_link_libraries (fntt_test nttw ${PLATFORM_LIBS})

add_executable (fntt_prime_test fntt_prime_test/fntt_prime_test.c)
target_link_libraries (fntt_prime_test nttw ${PLATFORM_LIBS}) 

add_executable (fntt_2D_test fntt_2D_test/fntt_2D_test.c)
target_link_libraries (fntt_2D_test nttw ${PLATFORM_LIBS}) 

add_executable (fntt_2D_prime_test fntt_2D_prime_test/fntt_2D_prime_test.c)
target_link_libraries (fntt_2D_prime_test nttw ${PLATFORM_LIBS}) 

add_executable (convolve_test convolve_test/convolve_test.c)
target_link_libraries (convolve_test nttw ${PLATFORM_LIBS}) 
