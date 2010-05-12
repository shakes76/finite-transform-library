### NTTW Project File

include(config.pri) # compiler flags

TEMPLATE = lib
TARGET = 
DESTDIR = lib
OBJECTS_DIR += obj

unix {
    INSTALLBASE    = /usr/local
}

win32 {
    INSTALLBASE    = C:/NTTW
}

DEPENDPATH += $$NTTW $$NTTW_INCLUDEPATH
INCLUDEPATH += $$NTTW_INCLUDEPATH
target.path    = $$INSTALLBASE/lib
headers.path   = $$INSTALLBASE/include/nttw
doc.path       = $$INSTALLBASE/doc

CONFIG += dll warn_on release

DEFINES += NTTW_MAKEDLL

message("NTTW Library Build")
HEADERS += include/global.h \
    include/timing.h \
	include/array.h \
	include/number32.h \
	include/image.h \
	include/prime.h 
SOURCES += src/array.c \
    src/timing.c \
	src/number32.c \
	src/image.c \
	src/prime.c

headers.files  = $$HEADERS
#doc.files      = $${DGV}/doc/html
INSTALLS = target headers