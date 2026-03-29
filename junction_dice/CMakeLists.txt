set(ICON_FILE ${CMAKE_CURRENT_SOURCE_DIR}/../icons/junction_dice.icns)
set_source_files_properties(${ICON_FILE} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")

qt_add_executable(junction_dice
    MACOSX_BUNDLE
    main.cpp
    mainwindow.cpp
    mainwindow.h
    jdworker.cpp
    jdworker.h
    signals.cpp
    signals.h
    datastructs.h
    mainwindow.ui
    ${ICON_FILE}
)

set_target_properties(junction_dice PROPERTIES
    OUTPUT_NAME "JunctionDice++"
    MACOSX_BUNDLE_ICON_FILE junction_dice.icns
)

if(WIN32)
    target_sources(junction_dice PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../icons/junction_make_windows.rc)
endif()

target_include_directories(junction_dice PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/zstr
)

target_link_libraries(junction_dice PRIVATE
    Qt6::Core
    Qt6::Gui
    Qt6::Widgets
    Qt6::Sql
    zlibstatic
)

add_dependencies(junction_dice ncbi_blast blat_build)

if(APPLE)
    add_custom_command(TARGET junction_dice POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
            "$<TARGET_BUNDLE_DIR:junction_dice>/Contents/Tools"
        COMMAND ${CMAKE_COMMAND} -E make_directory
            "$<TARGET_BUNDLE_DIR:junction_dice>/Contents/Data"
        # Copy tools (blastn, blat) and database configs
        COMMAND ${CMAKE_COMMAND}
            -DSRC=${CMAKE_SOURCE_DIR}
            -DDST=$<TARGET_BUNDLE_DIR:junction_dice>/Contents
            -P ${CMAKE_SOURCE_DIR}/cmake/bundle_junction_dice_data.cmake
        COMMENT "Bundling JunctionDice++ tools and data"
    )
endif()
