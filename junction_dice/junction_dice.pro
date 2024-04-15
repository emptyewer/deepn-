QT       += core gui sql

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += sdk_no_version_check
CONFIG += c++20

TARGET = JunctionDice++

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

# DESTDIR = ../pythonqt/lib

#include (../build/common.prf)
#include (../build/PythonQt.prf)
#include (../build/PythonQt_QtAll.prf)

SOURCES += \
    jdworker.cpp \
    main.cpp \
    mainwindow.cpp \
    signals.cpp

HEADERS += \
    datastructs.h \
    jdworker.h \
    mainwindow.h \
    signals.h \
    zstr/strict_fstream.hpp \
    zstr/zstr.hpp

FORMS += \
    mainwindow.ui

# Default rules for deployment.
#qnx: target.path = /tmp/$${TARGET}/bin
#else: unix:!android: target.path = /opt/$${TARGET}/bin
#!isEmpty(target.path): INSTALLS += target

# Default rules for deployment
macx {
    ICON = ../icons/junction_dice.icns
    QMAKE_MACOSX_DEPLOYMENT_TARGET = 14.0
}

win32 {
    ICON = ../icons/junction_dice.ico
}

unix {
}

macx: LIBS += -L$$PWD/../zlib/build/ -lz.1.2.13
INCLUDEPATH += $$PWD/../zlib
DEPENDPATH += $$PWD/../zlib


