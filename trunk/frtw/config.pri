### Ghost configuration file that contains defines and library variables

#FRTW Library
FRTW_INCLUDEPATH +=  include
FRTW += .

FRTW_LIBSPATH += lib
unix:FRTW_LIBS += -lfrtw
win32-msvc2005:FRTW_LIBS += frtw.lib
win32-g++:FRTW_LIBS += -lfrtw

#NTTW Library (EXTERNAL)
NTTW_LIBSPATH += ../nttw2/lib
unix:NTTW_LIBS += -lnttw
win32-msvc2005:NTTW_LIBS += nttw.lib
win32-g++:NTTW_LIBS += -lnttw

CONFIG -= qt
unix:DEFINES -= QT_SHARED
DEFINES -= QT_DLL
win32:DEFINES += _CRT_SECURE_NO_WARNINGS

#FFTW (EXTERNAL)
unix:FFTW_LIBS = -lfftw3
win32-msvc2005:FFTW_LIBS = fftw3.lib
win32-g++:FFTW_LIBS = -lfftw3

LIBS += $$FFTW_LIBS $$NTTW_LIBS

DEFINES += NTTW_DLL

include(cflags.pri) # compiler flags
include(cxxflags.pri) # compiler flags