cd ../statgen
make -j12
cd ../samtools/htslib
if [ ! -f "configure" ]; then
    autoreconf -i
    ./configure
fi
cd ..
if [ ! -f "configure" ]; then
    autoreconf -i
    ./configure
fi
make all all-htslib
