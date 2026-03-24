source /etc/environment

COMPILED_OUTPUT=${PROJECT_FILE%%.*}_linux_64bit.7z
# Compile and build DEEPN Linux 1
echo "Source ROOT directory: " $ROOT_DIR
echo "Project filename: " $PROJECT_FILE
echo "Compiled executable filename: " $COMPILED_OUTPUT
BUILD_DIR=$ROOT_DIR/build_linux_release
rm -rf $BUILD_DIR && mkdir -p $BUILD_DIR

conan install -b "*" -s build_type=Release -if $BUILD_DIR -g qmake $ROOT_DIR/conanfile.txt
/opt/qt5.15.2/bin/qmake -o $BUILD_DIR/Makefile $ROOT_DIR/$PROJECT_FILE -spec linux-g++-64 CONFIG+=x86_64 CONFIG-=qtquickcompiler
make -j4 -C $BUILD_DIR qmake_all
for d in $BUILD_DIR/* ; do
    if [ -d "$d" ]; then
        make -j4 -C $d
    fi
done
mkdir -p $BUILD_DIR/builds
for f in $BUILD_DIR/**/* ; do
    if [ -x "$f" ]; then
        cp $f $BUILD_DIR/builds/
    fi
done
7z a -m0=ppmd -mx=9 $ROOT_DIR/$COMPILED_OUTPUT $BUILD_DIR/builds/*