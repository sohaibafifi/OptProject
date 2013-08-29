#-------------------------------------------------
#
# Project created by QtCreator 2013-08-29T14:35:10
#
#-------------------------------------------------

QT       -= core gui

TARGET = OptProject
TEMPLATE = lib
CONFIG += staticlib

SOURCES += \
    src/Utilities/Timer.cpp \
    src/Utilities/Algorithms.cpp \
    src/Utilities/mtrand/mtrand.cpp \
    src/Utilities/getoptpp/getopt_pp.cpp

HEADERS += src/OptProject.hpp \
    src/Models/IProblem.hpp \
    src/Solvers/ISolver.hpp \
    src/Models/ISolution.hpp \
    src/Utilities/Algorithms.hpp \
    src/Utilities/Timer.hpp \
    src/Utilities/mtrand/mtrand.h \
    src/Utilities/getoptpp/getopt_pp.h \
    src/libs.hpp
unix:!symbian {
    maemo5 {
        target.path = /opt/usr/lib
    } else {
        target.path = /usr/lib
    }
    INSTALLS += target
}

OTHER_FILES += \
    premake4.lua
