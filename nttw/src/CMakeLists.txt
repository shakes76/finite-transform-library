# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Add Sources
set(NTTW_SOURCES timing array image prime number)
# Add Headers into List
file(GLOB NTTW_HEADERS ${NTTW_SOURCE_DIR}/include/nttw/*.h) 

# 64-bit Build?
IF(CMAKE_SIZEOF_VOID_P EQUAL 8) #64-bit environment
  set (NTTW_64 ON) #True
  message("Detected 64-bit system. Using appropriate 64-bit integers.")
ELSEIF(CMAKE_SIZEOF_VOID_P EQUAL 8) #64-bit environment
  set (NTTW_64 OFF) #False
ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 8) #64-bit environment

# Create a library called "nttw" which includes the source files.
# The extension is already found. Any number of sources could be listed here.
include_directories (${NTTW_SOURCE_DIR}/include) 
add_library (nttw ${NTTW_SOURCES}) 

#MSVC
if(MSVC)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS) #remove "unsafe functions" warnings
endif(MSVC)

# Set the build version (VERSION) and the API version (SOVERSION)
SET_TARGET_PROPERTIES(nttw
                      PROPERTIES
                      VERSION 1.2.3
                      SOVERSION 1)

# Installation of the Library
# Installs the header files
install(FILES ${NTTW_HEADERS} DESTINATION "include/nttw")
# Installs the library
install(TARGETS nttw
        DESTINATION "lib"
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)