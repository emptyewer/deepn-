# zlib
cd ../zlib
if [ ! -d "build " ]; then
    echo "Configuring zlib"
    mkdir build
    cd build
    /usr/local/bin/cmake ..
    cd ..
fi
echo "Building zlib"
make

# ncbi
cd ../ncbi
if [ ! -f "Makefile" ]; then
    echo "Configuring NCBI"
    ./configure --with-static --with-mt --with-64 --with-static-exe --with-build-root=./build --without-boost --without-gui
fi
echo "Building NCBI"
make -j12

# samtools
cd ../samtools/htslib
if [ ! -f "configure" ]; then
    echo "Configuring htslib"
    autoreconf -i
    ./configure
fi

cd ..
if [ ! -f "configure" ]; then
    echo "Configuring samtools"
    autoreconf -i
    ./configure
fi
echo "Building samtools"
make all all-htslib