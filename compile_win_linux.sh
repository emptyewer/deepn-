# Download all the git


# Compile Install Linux dependencies
ROOT_DIR="/builds/emptyewer/deepn-plus/libs"
BUILD_DIR=build_linux
for d in bpp-core bpp-seq bpp-seq-omics
do
    echo $ROOT_DIR"/"$d
    rm -rf $ROOT_DIR"/"$d"/"$BUILD_DIR
    mkdir -p $BUILD_DIR
    cmake -S $ROOT_DIR"/"$d -B $ROOT_DIR"/"$d"/"$BUILD_DIR -DCMAKE_INSTALL_PREFIX=$BPP_DIR -DBUILD_STATIC=TRUE
    make -j8 -C $ROOT_DIR"/"$d"/"$BUILD_DIR
done

# docker run --rm -v $(pwd):/builds/emptyewer/deepn-plus/ venkykrishna/deepn:latest
