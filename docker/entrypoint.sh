source /etc/environment
cd /builds/emptyewer/deepn-plus/

qmake deepn.pro
make -j8
mkdir -p deepn++_linux_64bit
mv main.o mainwindow.o moc_* ui_* Makefile DEEPN++ deepn++_linux_64bit/

/usr/lib/mxe/usr/x86_64-w64-mingw32.static/qt5/bin/qmake deepn.pro
make -j8
rm -rf Makefile* deepn++_plugin_import.cpp ui_mainwindow.h debug
mv release deepn++_windows_64bit

7z a deepn++_linux_64bit.7z deepn++_linux_64bit/*
7z a deepn++_windows_64bit.7z deepn++_windows_64bit/*
