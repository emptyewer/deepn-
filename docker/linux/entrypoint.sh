source /etc/environment

# Compile and build DEEPN Linux
rm -rf $CI_PROJECT_DIR/build_linux_release && mkdir -p $CI_PROJECT_DIR/build_linux_release
cmake -DCMAKE_BUILD_TYPE=Release -S $CI_PROJECT_DIR -B $CI_PROJECT_DIR/build_linux_release
make -C $CI_PROJECT_DIR/build_linux_release
7z a -m0=ppmd -mx=9 $CI_PROJECT_DIR/deepn++_linux_64bit.7z $CI_PROJECT_DIR/build_linux_release{/**,}/*++