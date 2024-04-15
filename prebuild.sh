# zlib

cd ../zlib
if [ ! -d "build "]; then
    mkdir build
    cd build
    cmake ..
    make
fi

# ncbi
cd ../ncbi
if [ ! -f "Makefile" ]; then
    ./configure --with-static --with-mt --with-64 --with-static-exe --with-build-root=./build --without-boost --without-gui
    make -j12
fi

# samtools
cd ../samtools/htslib
if [ ! -f "configure" ]; then
    autoreconf -i
    ./configure

cd ..
if [ ! -f "configure" ]; then
    autoreconf -i
    ./configure
fi
make all all-htslib