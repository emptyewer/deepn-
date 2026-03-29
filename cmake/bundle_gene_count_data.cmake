# Bundle data files into GeneCount++.app
# Called via: cmake -P bundle_gene_count_data.cmake -DSRC=<data_dir> -DDST=<bundle_Contents>

# Unzip genome databases
file(GLOB SQLITE_ZIPS "${SRC}/*.sqlite.zip")
foreach(ZIP ${SQLITE_ZIPS})
    execute_process(
        COMMAND unzip -oq "${ZIP}" -d "${DST}/Data"
        RESULT_VARIABLE RES
    )
    if(NOT RES EQUAL 0)
        message(WARNING "Failed to unzip ${ZIP}")
    endif()
endforeach()

# Remove __MACOSX junk from unzip
file(REMOVE_RECURSE "${DST}/Data/__MACOSX")

# Unzip gene lists
file(MAKE_DIRECTORY "${DST}/Data/gene_list")
file(GLOB PRN_ZIPS "${SRC}/gene_list/*.prn.zip")
foreach(ZIP ${PRN_ZIPS})
    execute_process(
        COMMAND unzip -oq "${ZIP}" -d "${DST}/Data/gene_list"
        RESULT_VARIABLE RES
    )
    if(NOT RES EQUAL 0)
        message(WARNING "Failed to unzip ${ZIP}")
    endif()
endforeach()

file(REMOVE_RECURSE "${DST}/Data/gene_list/__MACOSX")

# Copy gene dictionaries (real pickle files, not LFS pointers)
file(COPY "${SRC}/gene_dictionary" DESTINATION "${DST}/Data")

# Copy JSON config files
file(GLOB JSON_FILES "${SRC}/*.json")
foreach(F ${JSON_FILES})
    file(COPY "${F}" DESTINATION "${DST}/Data")
endforeach()
