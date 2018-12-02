### NTTW configuration file that contains defines and library variables

NTTW_INCLUDEPATH +=  include
NTTW += .

NTTW_LIBSPATH += lib
unix:NTTW_LIBS += -lnttw
win32-msvc2005:NTTW_LIBS += nttw.lib
win32-g++:NTTW_LIBS += -lnttw

CONFIG -= qt
unix:DEFINES -= QT_SHARED
DEFINES -= QT_DLL
win32:DEFINES += _CRT_SECURE_NO_WARNINGS

DEFINES += NTTW_DLL

include(../cflags.pri) # compiler flags
include(../cxxflags.pri) # compiler flags

#32/64-bit Section
DEFINES += NTTW_64 #Uncomment to build 64-bit version
NTTW_BUILD_TYPE = 32-bit
contains(DEFINES,NTTW_64) {
    NTTW_BUILD_TYPE = 64-bit
}
message(Build Type: $$NTTW_BUILD_TYPE)
