######################################################################
# Automatically generated by qmake (2.01a) Fri Jul 17 10:46:33 2009
######################################################################
include(../../config.pri)

TEMPLATE = app
TARGET = 
DEPENDPATH += ../../$$FRTW_INCLUDEPATH
INCLUDEPATH += ../../$$FRTW_INCLUDEPATH
DESTDIR = ../../bin
OBJECTS_DIR += ../../obj

CONFIG += warn_on release console

unix:LIBS += -L../../$$FRTW_LIBSPATH $$FRTW_LIBS
win32-msvc2005:LIBS += /LIBPATH:"../../$$FRTW_LIBSPATH" $$FRTW_LIBS
win32-g++:LIBS += -L../../$$FRTW_LIBSPATH $$FRTW_LIBS

# Input
SOURCES += machine_properties.cpp
