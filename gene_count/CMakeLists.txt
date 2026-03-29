set(ICON_FILE ${CMAKE_CURRENT_SOURCE_DIR}/../icons/gene_count.icns)
set_source_files_properties(${ICON_FILE} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")

qt_add_executable(gene_count
    MACOSX_BUNDLE
    main.cpp
    mainwindow.cpp
    mainwindow.h
    gcworker.cpp
    gcworker.h
    maphits.cpp
    maphits.h
    signals.cpp
    signals.h
    datastructs.h
    mainwindow.ui
    ${ICON_FILE}
)

set_target_properties(gene_count PROPERTIES
    OUTPUT_NAME "GeneCount++"
    MACOSX_BUNDLE_ICON_FILE gene_count.icns
)

if(WIN32)
    target_sources(gene_count PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../icons/gene_count_windows.rc)
endif()

target_link_libraries(gene_count PRIVATE
    Qt6::Core
    Qt6::Gui
    Qt6::Widgets
    Qt6::Sql
    QSimpleXlsxWriter
)

if(APPLE)
    add_custom_command(TARGET gene_count POST_BUILD
        # Python scripts
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${CMAKE_SOURCE_DIR}/scripts/GeneCount++"
            "$<TARGET_BUNDLE_DIR:gene_count>/Contents/Scripts"
        # Unzip and copy genome data (sqlite, gene lists, dictionaries, configs)
        COMMAND ${CMAKE_COMMAND}
            -DSRC=${CMAKE_SOURCE_DIR}/data
            -DDST=$<TARGET_BUNDLE_DIR:gene_count>/Contents
            -P ${CMAKE_SOURCE_DIR}/cmake/bundle_gene_count_data.cmake
        COMMENT "Bundling GeneCount++ scripts and data"
    )
endif()
