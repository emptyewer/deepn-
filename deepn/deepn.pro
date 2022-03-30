QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG+=sdk_no_version_check
CONFIG += c++17 conan_basic_setup

TARGET = DEEPN

include(../config.pri)
# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

include(../conanbuildinfo.pri)

SOURCES += \
    main.cpp \
    mainwindow.cpp \
    dqtextedit.cpp

HEADERS += \
    mainwindow.h \
    dqtextedit.h

FORMS += \
    mainwindow.ui

# Default rules for deployment
macx {
    ICON = ../icons/deepn.icns
    QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.15
    QMAKE_LFLAGS += -Bstatic
}

win32 {
    ICON = ../icons/deepn.ico
    QMAKE_LFLAGS += -static
}

unix {
}
