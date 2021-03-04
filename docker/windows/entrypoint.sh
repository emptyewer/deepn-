source /etc/environment

# Compile and build DEEPN Windows
rm -rf $CI_PROJECT_DIR/build_windows_release && mkdir -p $CI_PROJECT_DIR/build_windows_release
/usr/lib/mxe/usr/bin/x86_64-w64-mingw32.static-cmake -DCMAKE_BUILD_TYPE=Release -S $CI_PROJECT_DIR -B $CI_PROJECT_DIR/build_windows_release
make -C $CI_PROJECT_DIR/build_windows_release
7z a -m0=ppmd -mx=9 $CI_PROJECT_DIR/deepn++_windows_64bit.7z $CI_PROJECT_DIR/build_windows_release{/**,}/*++.exe