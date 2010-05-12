### C Ghosts Project File

include(config.pri) # compiler flags

TEMPLATE = lib
TARGET = 
DESTDIR = lib
OBJECTS_DIR += obj

unix {
    INSTALLBASE    = /usr/local
}

win32 {
    INSTALLBASE    = C:/FRTW
}

DEPENDPATH += $$FRTW $$FRTW_INCLUDEPATH
INCLUDEPATH += $$FRTW_INCLUDEPATH
target.path    = $$INSTALLBASE/lib
headers.path   = $$INSTALLBASE/include
doc.path       = $$INSTALLBASE/doc

CONFIG += dll warn_on release

DEFINES += NTTW_MAKEDLL

message("FRTW C Library Build")
HEADERS += include/array_complex.h \
        include/radon.h \
	include/fourier.h
SOURCES += src/array_complex.c \
        src/radon.c \
	src/fourier.c

headers.files  = $$HEADERS
#doc.files      = $${DGV}/doc/html
INSTALLS = target headers