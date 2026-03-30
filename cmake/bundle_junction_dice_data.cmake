# Bundle data and tools into JunctionDice++.app
# Called via: cmake -P bundle_junction_dice_data.cmake -DSRC=<project_root> -DDST=<bundle_Contents>

# Copy blastn (from ncbi prebuild, may not exist yet)
set(BLASTN "${SRC}/ncbi/build/bin/blastn")
if(EXISTS "${BLASTN}")
    file(COPY "${BLASTN}" DESTINATION "${DST}/Tools")
else()
    message(STATUS "blastn not found at ${BLASTN} — run prebuild.sh first")
endif()

# Copy blat
set(BLAT "${SRC}/blat/bin/blat")
if(EXISTS "${BLAT}")
    file(COPY "${BLAT}" DESTINATION "${DST}/Tools")
else()
    message(STATUS "blat not found at ${BLAT}")
endif()

# Copy blast databases (if populated)
file(GLOB BLAST_DB_FILES "${SRC}/data/blast_db/*")
if(BLAST_DB_FILES)
    file(COPY ${BLAST_DB_FILES} DESTINATION "${DST}/Data")
endif()

# Copy deepn.json config
file(COPY "${SRC}/data/deepn.json" DESTINATION "${DST}/Data")

# Unzip pre-built gene annotation SQLite databases
file(GLOB ANNOT_ZIPS "${SRC}/data/*.fasta.sqlite.zip")
foreach(Z ${ANNOT_ZIPS})
    get_filename_component(ZNAME "${Z}" NAME_WE)  # e.g. hg38GeneList2023.unique.fasta.sqlite
    set(SQLITE_OUT "${DST}/Data/${ZNAME}")
    if(NOT EXISTS "${SQLITE_OUT}")
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E tar xf "${Z}"
            WORKING_DIRECTORY "${DST}/Data"
        )
        # Remove macOS resource fork junk if present
        file(REMOVE_RECURSE "${DST}/Data/__MACOSX")
    endif()
endforeach()

# Copy genome FASTA JSON configs
file(GLOB FASTA_JSONS "${SRC}/data/*.unique.fasta.json")
foreach(F ${FASTA_JSONS})
    file(COPY "${F}" DESTINATION "${DST}/Data")
endforeach()
