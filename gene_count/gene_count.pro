QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG+=sdk_no_version_check
CONFIG += c++17 conan_basic_setup

TARGET = GeneCount

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0
include(../conanbuildinfo.pri)


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
    QMAKE_LFLAGS += -Bstatic
    LIBS += -L/usr/local/opt/python@3.9/Frameworks/Python.framework/Versions/3.9/lib/python3.9/config-3.9-darwin -lpython3.9
    INCLUDEPATH += /usr/local/opt/python@3.9/Frameworks/Python.framework/Versions/3.9/include/python3.9
    DEPENDPATH += /usr/local/opt/python@3.9/Frameworks/Python.framework/Versions/3.9/include/python3.9
    PRE_TARGETDEPS += /usr/local/opt/python@3.9/Frameworks/Python.framework/Versions/3.9/lib/python3.9/config-3.9-darwin/libpython3.9.a
}

win32 {
    ICON = ../icons/gene_count.ico
    QMAKE_LFLAGS += -static
}

unix {
}

DISTFILES += \
    python/another.py \
    python/gene_count.py

