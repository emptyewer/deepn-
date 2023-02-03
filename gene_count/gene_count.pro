QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += sdk_no_version_check
CONFIG += c++17 conan_basic_setup

TARGET = GeneCount

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

DESTDIR = ../pythonqt/lib

include (../conanbuildinfo.pri)
include (../build/common.prf)
include (../build/PythonQt.prf)
include (../build/PythonQt_QtAll.prf)

SOURCES += \
    main.cpp \
    mainwindow.cpp \


HEADERS += \
    mainwindow.h \

FORMS += \
    mainwindow.ui

# Default rules for deployment
macx {
    ICON = ../icons/gene_count.icns
    QMAKE_MACOSX_DEPLOYMENT_TARGET = 12.0
}

win32 {
    ICON = ../icons/gene_count.ico
}

unix {
}


#copydata.commands = $(COPY_DIR) $$PWD/python/*.py $$OUT_PWD/GeneCount.app/Contents/MacOS/
#first.depends = $(first) copydata
#export(first.depends)
#export(copydata.commands)
#QMAKE_EXTRA_TARGETS += first copydata

