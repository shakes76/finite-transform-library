# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Add Sources
set(FRTW_SOURCES array_complex vector mojette noise fourier radon ghosts)
# Add Headers into List
file(GLOB FRTW_HEADERS ${FRTW_SOURCE_DIR}/include/*.h) 

# Create a library called "FRTW" which includes the source files.
# The extension is already found. Any number of sources could be listed here.
include_directories (${FRTW_SOURCE_DIR}/include ${NTTW_SOURCE_DIR}/include) 
add_library (frtw ${FRTW_SOURCES}) 

#MSVC
if(MSVC)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS) #remove "unsafe functions" warnings
endif(MSVC)

# Set the build version (VERSION) and the API version (SOVERSION)
SET_TARGET_PROPERTIES(frtw
                      PROPERTIES
                      VERSION 1.50
                      SOVERSION 1)
		      
#Additional Platform specific libs
if(CMAKE_COMPILER_IS_GNUCC)
	set (PLATFORM_LIBS rt)
endif(CMAKE_COMPILER_IS_GNUCC)

# Link the executable to the NTTW library.
target_link_libraries (frtw nttw ${FFTW_LIBRARIES} ${PLATFORM_LIBS}) 

# Installation of the Library
# Installs the header files
install(FILES ${FRTW_HEADERS} DESTINATION "include")
# Installs the library
install(TARGETS frtw
        DESTINATION "lib"
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)