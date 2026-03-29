
/bin/mkdir -p $(pwd)/deepn/DEEPN++.app/Contents/Data
rsync -avhu --ignore-errors ../data/deepn.json $(pwd)/deepn/DEEPN++.app/Contents/Data

## Uncomment for deployment ###

# for db in ../data/blast_db/*.unique.fasta; do
#     ../ncbi/build/bin/makeblastdb -in ${db} -dbtype nucl
#     ../ncbi/build/bin/makembindex -input ${db}
# done

files=(GeneCount++ JunctionDice++)
for file in "${files[@]}"; do  # loop through the array
    /bin/mkdir -p $(pwd)/staging/${file}.app/Contents/Tools
    /bin/mkdir -p $(pwd)/staging/${file}.app/Contents/Data
    if [[ ${file} == JunctionDice++ ]]; then
        rsync -avhu --ignore-errors ../ncbi/build/bin/blastn $(pwd)/staging/${file}.app/Contents/Tools
        rsync -avhu --ignore-errors ../blat/macos.x86_64/blat $(pwd)/staging/${file}.app/Contents/Tools
        rsync -avhu --ignore-errors ../data/blast_db/*.* $(pwd)/staging/${file}.app/Contents/Data
    else
        # scripts & data
        /bin/mkdir -p $(pwd)/staging/${file}.app/Contents/Scripts/Data
        rsync -avhu --ignore-errors ../scripts/${file}/* $(pwd)/staging/${file}.app/Contents/Scripts
        rsync -avhu --ignore-errors ../data/deepn.sqlite $(pwd)/staging/${file}.app/Contents/Scripts/Data
    fi
done

rsync -avhu --ignore-errors $(pwd)/staging/GeneCount++.app  $(pwd)/deepn/DEEPN++.app/Contents/Resources
rsync -avhu --ignore-errors $(pwd)/staging/JunctionDice++.app  $(pwd)/deepn/DEEPN++.app/Contents/Resources
rsync -avhu --ignore-errors $(pwd)/query/MultiQuery++.app  $(pwd)/deepn/DEEPN++.app/Contents/Resources
rsync -avhu --ignore-errors $(pwd)/read_depth/ReadDepth++.app  $(pwd)/deepn/DEEPN++.app/Contents/Resources