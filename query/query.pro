QT       += core gui charts

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += sdk_no_version_check
CONFIG += c++20


TARGET = MultiQuery++

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    dockwidget.cpp \
    main.cpp \
    mainwindow.cpp

HEADERS += \
    dockwidget.h \
    mainwindow.h

FORMS += \
    dockwidget.ui \
    mainwindow.ui

# Default rules for deployment
macx {
    ICON = ../icons/multi_query.icns
    QMAKE_MACOSX_DEPLOYMENT_TARGET = 14.0
}

win32 {
    ICON = ../icons/multi_query.ico
}


# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
