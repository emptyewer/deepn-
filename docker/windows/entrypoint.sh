source /etc/environment

ROOT_DIR=/builds/emptyewer/deepn-plus
# Compile and build DEEPN Windows
rm -rf $ROOT_DIR/build_windows_release && mkdir -p $ROOT_DIR/build_windows_release
/usr/lib/mxe/usr/bin/x86_64-w64-mingw32.static-cmake -DCMAKE_BUILD_TYPE=Release -S $ROOT_DIR -B $ROOT_DIR/build_windows_release
make -C $ROOT_DIR/build_windows_release
7z a -m0=ppmd -mx=9 $ROOT_DIR/deepn++_windows_64bit.7z $ROOT_DIR/build_windows_release{/**,}/*++.exe