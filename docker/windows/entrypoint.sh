source /etc/environment

COMPILED_OUTPUT=${PROJECT_FILE%%.*}_windows_64bit.7z
# Compile and build DEEPN Windows
echo "Source ROOT directory: " $ROOT_DIR
echo "Project filename: " $PROJECT_FILE
echo "Compiled executable filename: " $COMPILED_OUTPUT
BUILD_DIR=$ROOT_DIR/build_windows_release
rm -rf $ROOT_DIR/build_windows_release && mkdir -p $BUILD_DIR
/mxe/usr/bin/x86_64-w64-mingw32.static-qmake-qt5 -o $BUILD_DIR/Makefile $ROOT_DIR/$PROJECT_FILE -spec win32-g++ CONFIG+=x86_64 CONFIG-=qtquickcompiler
make -j4 -C $BUILD_DIR qmake_all
for d in $BUILD_DIR/* ; do
    if [ -d "$d" ]; then
        make -j4 -C $d
    fi
done
mkdir -p $BUILD_DIR/builds
for f in $BUILD_DIR/*/*/*.exe ; do
    if [ -x "$f" ]; then
        cp $f $BUILD_DIR/builds/
    fi
done
7z a -m0=ppmd -mx=9 $ROOT_DIR/$COMPILED_OUTPUT $BUILD_DIR/builds/*