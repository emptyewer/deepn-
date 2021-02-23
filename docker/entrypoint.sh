source /etc/environment

# ROOT_DIR="/builds/emptyewer/deepn-plus/libs"
# # Compile Install Linux dependencies
# BUILD_DIR=build_linux
# for d in bpp-core bpp-seq bpp-seq-omics
# do
#     echo $ROOT_DIR"/"$d
#     rm -rf $ROOT_DIR"/"$d"/"$BUILD_DIR
#     mkdir -p $BUILD_DIR
#     cmake -DCMAKE_BUILD_TYPE=Release -S $ROOT_DIR"/"$d -B $ROOT_DIR"/"$d"/"$BUILD_DIR -DCMAKE_INSTALL_PREFIX=$BPP_DIR -DBUILD_STATIC=TRUE
#     make -j2 -C $ROOT_DIR"/"$d"/"$BUILD_DIR
# done

# # Compile Install Windows dependencies
# BUILD_DIR=build_windows
# for d in bpp-core bpp-seq bpp-seq-omics
# do
#     echo $ROOT_DIR"/"$d
#     rm -rf $ROOT_DIR"/"$d"/"$BUILD_DIR
#     mkdir -p $BUILD_DIR
#     /usr/lib/mxe/usr/bin/x86_64-w64-mingw32.static-cmake -DCMAKE_BUILD_TYPE=Release -S $ROOT_DIR"/"$d -B $ROOT_DIR"/"$d"/"$BUILD_DIR -DCMAKE_INSTALL_PREFIX=$BPP_DIR -DBUILD_STATIC=TRUE
#     make -j2 -C $ROOT_DIR"/"$d"/"$BUILD_DIR
# done

# Compile and build DEEPN Linux
ROOT_DIR=/builds/emptyewer/deepn-plus
cmake -DCMAKE_BUILD_TYPE=Release -S $ROOT_DIR -B $ROOT_DIR/build_linux_release
make -j2 -C $ROOT_DIR/build_linux_release
7z a $ROOT_DIR/deepn++_linux_64bit.7z $ROOT_DIR/build_linux_release/DEEPN++

/usr/lib/mxe/usr/bin/x86_64-w64-mingw32.static-cmake -DCMAKE_BUILD_TYPE=Release -S $ROOT_DIR -B $ROOT_DIR/build_windows_release
make -C $ROOT_DIR/build_windows_release
7z a $ROOT_DIR/deepn++_windows_64bit.7z $ROOT_DIR/build_windows_release/DEEPN++.exe
