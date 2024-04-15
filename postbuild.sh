
/bin/mkdir -p $(pwd)/deepn/DEEPN++.app/Contents/Data
rsync -avhu --ignore-errors ../data/deepn.json $(pwd)/deepn/DEEPN++.app/Contents/Data

## Uncomment for deployment ###

# for db in ../data/blast_db/*.unique.fasta; do
#     ../ncbi/build/bin/makeblastdb -in ${db} -dbtype nucl
#     ../ncbi/build/bin/makembindex -input ${db}
# done

files=(GeneCount++ JunctionDice++)
for file in "${files[@]}"; do  # loop through the array
    # samtools
    /bin/mkdir -p $(pwd)/pythonqt/lib/${file}.app/Contents/Tools
    /bin/mkdir -p $(pwd)/pythonqt/lib/${file}.app/Contents/Data
    if [[ ${file} == JunctionDice++ ]]; then
        rsync -avhu --ignore-errors ../ncbi/build/bin/blastn $(pwd)/pythonqt/lib/${file}.app/Contents/Tools
        rsync -avhu --ignore-errors ../blat/macos.x86_64/blat $(pwd)/pythonqt/lib/${file}.app/Contents/Tools
        rsync -avhu --ignore-errors ../data/blast_db/*.* $(pwd)/pythonqt/lib/${file}.app/Contents/Data
    else
        # rsync -avhu --ignore-errors ../samtools/samtools $(pwd)/pythonqt/lib/${file}.app/Contents/Tools
        # scripts & data
        /bin/mkdir -p $(pwd)/pythonqt/lib/${file}.app/Contents/Scripts/Data
        rsync -avhu --ignore-errors ../scripts/${file}/* $(pwd)/pythonqt/lib/${file}.app/Contents/Scripts
        rsync -avhu --ignore-errors ../data/deepn.sqlite $(pwd)/pythonqt/lib/${file}.app/Contents/Scripts/Data
        # HTSLIB
        rsync -avhu --ignore-errors ../samtools/htslib/*.dylib $(pwd)/pythonqt/lib/${file}.app/Contents/MacOS
        # change lib loader paths
        install_name_tool -change /usr/local/lib/libhts.3.dylib @executable_path/libhts.3.dylib $(pwd)/pythonqt/lib/${file}.app/Contents/MacOS/${file}
    fi
done

rsync -avhu --ignore-errors $(pwd)/pythonqt/lib/GeneCount++.app  $(pwd)/deepn/DEEPN++.app/Contents/Resources
rsync -avhu --ignore-errors $(pwd)/pythonqt/lib/JunctionDice++.app  $(pwd)/deepn/DEEPN++.app/Contents/Resources
rsync -avhu --ignore-errors $(pwd)/query/MultiQuery++.app  $(pwd)/deepn/DEEPN++.app/Contents/Resources
rsync -avhu --ignore-errors $(pwd)/read_depth/ReadDepth++.app  $(pwd)/deepn/DEEPN++.app/Contents/Resources