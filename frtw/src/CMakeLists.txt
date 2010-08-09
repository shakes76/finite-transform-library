# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Add Sources
set(FRTW_SOURCES array_complex vector mojette noise fourier radon)
# Add Headers into List
file(GLOB FRTW_HEADERS ${FRTW_SOURCE_DIR}/include/*.h) 

# 64-bit Build?
set (NTTW_64 ON) #True
#set (NTTW_64 OFF) #False

#Add compiler flags
add_definitions (-DNTTW_DLL -DNTTW_MAKEDLL)
if (NTTW_64)
	add_definitions (-DNTTW_64)
endif (NTTW_64)

# Create a library called "FRTW" which includes the source files.
# The extension is already found. Any number of sources could be listed here.
include_directories (${FRTW_SOURCE_DIR}/include) 
add_library (frtw SHARED ${FRTW_SOURCES}) 

# Set the build version (VERSION) and the API version (SOVERSION)
SET_TARGET_PROPERTIES(frtw
                      PROPERTIES
                      VERSION 1.3.5
                      SOVERSION 1)
		      
#Additional Platform specific libs
if(CMAKE_COMPILER_IS_GNUCC)
	set (PLATFORM_LIBS rt)
endif(CMAKE_COMPILER_IS_GNUCC)

# Link the executable to the NTTW library.
target_link_libraries (frtw nttw fftw3 ${PLATFORM_LIBS}) 

# Installation of the Library
# Installs the header files
install(FILES ${FRTW_HEADERS} DESTINATION "include")
# Installs the library
install(TARGETS frtw
        DESTINATION "lib"
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)