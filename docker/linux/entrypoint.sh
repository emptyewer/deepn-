source /etc/environment

# Compile and build DEEPN Linux 
ROOT_DIR=/builds/emptyewer/deepn-plus
rm -rf $ROOT_DIR/build_linux_release && mkdir -p $ROOT_DIR/build_linux_release
cmake -DCMAKE_BUILD_TYPE=Release -S $ROOT_DIR -B $ROOT_DIR/build_linux_release
make -C $ROOT_DIR/build_linux_release
7z a -m0=ppmd -mx=9 $ROOT_DIR/deepn++_linux_64bit.7z $ROOT_DIR/build_linux_release{/**,}/*++