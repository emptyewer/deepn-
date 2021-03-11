source /etc/environment

# Compile and build DEEPN Linux 
ROOT_DIR=/builds/emptyewer/deepn-plus
rm -rf $ROOT_DIR/build_linux_release && mkdir -p $ROOT_DIR/build_linux_release
/opt/qt5.15.2/bin/qmake -o $ROOT_DIR/build_linux_release/Makefile $ROOT_DIR/deepn++.pro -spec linux-g++-64 CONFIG+=x86_64 CONFIG-=qtquickcompiler
make -j4 -C $ROOT_DIR/build_linux_release qmake_all
for d in $ROOT_DIR/build_linux_release/* ; do
    if [ -d "$d" ]; then
        make -j4 -C $d
    fi
done
mkdir -p $ROOT_DIR/build_linux_release/builds
for f in $ROOT_DIR/build_linux_release/**/* ; do
    if [ -x "$f" ]; then
        cp $f $ROOT_DIR/build_linux_release/builds/
    fi
done
7z a -m0=ppmd -mx=9 $ROOT_DIR/deepn++_linux_64bit.7z $ROOT_DIR/build_linux_release/builds/*