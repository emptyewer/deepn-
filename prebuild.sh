SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# zlib
cd "$SCRIPT_DIR"
git submodule update --init zlib
cd zlib
if [ ! -d "build" ]; then
    echo "Configuring zlib"
    mkdir build
    cd build
    cmake ..
    cd ..
fi
echo "Building zlib"
cmake --build build

# ncbi
cd "$SCRIPT_DIR"
git submodule update --init ncbi
cd ncbi
if [ ! -f "build/build/Makefile" ]; then
    echo "Configuring NCBI"
    ./configure --with-static --with-mt --with-64 --with-static-exe --with-build-root=./build --without-boost --without-gui
fi
echo "Building NCBI"
make -C build/build -j12