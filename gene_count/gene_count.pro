QT       += core gui sql

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += sdk_no_version_check
CONFIG += c++20

TARGET = GeneCount++

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

# DESTDIR = ../pythonqt/lib


# include (../build/common.prf)
# include (../build/PythonQt.prf)
# include (../build/PythonQt_QtAll.prf)
include (../xlsx/QSimpleXlsxWriter/QSimpleXlsxWriter.pri)

SOURCES += \
    gcworker.cpp \
    main.cpp \
    mainwindow.cpp \
    maphits.cpp \
    signals.cpp


HEADERS += \
    datastructs.h \
    gcworker.h \
    mainwindow.h \
    maphits.h \
    signals.h

FORMS += \
    mainwindow.ui


# Default rules for deployment
macx {
    ICON = ../icons/gene_count.icns
    QMAKE_MACOSX_DEPLOYMENT_TARGET = 14.0
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

#macx: LIBS += -L$$PWD/../samtools/htslib/ -lhts
#INCLUDEPATH += $$PWD/../samtools/htslib
#DEPENDPATH += $$PWD/../samtools/htslib
#macx: PRE_TARGETDEPS += $$PWD/../samtools/htslib/libhts.a
