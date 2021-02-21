source /etc/environment

qmake deepn.pro
make -j8
mkdir -p build-deepn-Desktop_Qt_5_15_2_linux_64bit-Release
mv main.o mainwindow.o moc_* ui_* Makefile DEEPN++ build-deepn-Desktop_Qt_5_15_2_linux_64bit-Release/

/usr/lib/mxe/usr/x86_64-w64-mingw32.static/qt5/bin/qmake deepn.pro
make -j8
rm -rf Makefile* deepn++_plugin_import.cpp ui_mainwindow.h debug
mv release build-deepn-Desktop_Qt_5_15_2_win64_64bit-Release
