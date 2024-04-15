# zlib
cd ../zlib
if [ -d "build " ]; then
    echo "Cleaning zlib"
    rm -rf build
fi
make

# ncbi
cd ../ncbi
if [ ! -d "build" ]; then
    echo "Cleaning NCBI"
    rm -rf build
fi

# samtools
echo "Cleaning samtools"
make clean clean-htslib