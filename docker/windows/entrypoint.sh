source /etc/environment

ROOT_DIR=/builds/emptyewer/deepn-plus
# Compile and build DEEPN Windows 
rm -rf $ROOT_DIR/build_windows_release && mkdir -p $ROOT_DIR/build_windows_release
/mxe/usr/bin/x86_64-w64-mingw32.static-qmake-qt5 -o $ROOT_DIR/build_windows_release/Makefile $ROOT_DIR/deepn++.pro -spec win32-g++ CONFIG+=x86_64 CONFIG-=qtquickcompiler
make -j4 -C $ROOT_DIR/build_windows_release qmake_all
for d in $ROOT_DIR/build_windows_release/* ; do
    if [ -d "$d" ]; then
        make -j4 -C $d
    fi
done
mkdir -p $ROOT_DIR/build_windows_release/builds
for f in $ROOT_DIR/build_windows_release/*/*/*.exe ; do
    if [ -x "$f" ]; then
        cp $f $ROOT_DIR/build_windows_release/builds/
    fi
done
7z a -m0=ppmd -mx=9 $ROOT_DIR/deepn++_windows_64bit.7z $ROOT_DIR/build_windows_release/builds/*