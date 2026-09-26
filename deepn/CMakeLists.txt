set(ICON_FILE ${CMAKE_CURRENT_SOURCE_DIR}/../icons/deepn.icns)
set_source_files_properties(${ICON_FILE} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")

qt_add_executable(deepn
    MACOSX_BUNDLE
    main.cpp
    mainwindow.cpp
    mainwindow.h
    darkpaint.h
    dqtextedit.cpp
    dqtextedit.h
    mainwindow.ui
    deepn.qrc
    ${ICON_FILE}
)

set_target_properties(deepn PROPERTIES
    OUTPUT_NAME "DEEPN++"
    MACOSX_BUNDLE_ICON_FILE deepn.icns
    MACOSX_BUNDLE_BUNDLE_NAME "DEEPN++"
    MACOSX_BUNDLE_GUI_IDENTIFIER "edu.uiowa.deepn"
    MACOSX_BUNDLE_BUNDLE_VERSION "${PROJECT_VERSION}"
    MACOSX_BUNDLE_SHORT_VERSION_STRING "${PROJECT_VERSION}"
)

if(WIN32)
    target_sources(deepn PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../icons/deepn_windows.rc)
endif()

target_include_directories(deepn PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(deepn PRIVATE
    Qt6::Core
    Qt6::Gui
    Qt6::Widgets
    Qt6::Sql
)

add_dependencies(deepn junction_dice gene_count query read_depth deseq2_app)

if(APPLE)
    get_target_property(_qt_qmake Qt6::qmake IMPORTED_LOCATION)
    get_filename_component(_qt_bin "${_qt_qmake}" DIRECTORY)
    find_program(MACDEPLOYQT_EXECUTABLE macdeployqt HINTS "${_qt_bin}")

    if(NOT MACDEPLOYQT_EXECUTABLE)
        message(WARNING
            "macdeployqt not found near ${_qt_bin}. "
            "Qt frameworks will not be embedded — the app will not be self-contained.")
    endif()

    # macdeployqt path is exported so child CMakeLists can use it
    set(MACDEPLOYQT_EXECUTABLE "${MACDEPLOYQT_EXECUTABLE}" CACHE INTERNAL "")

    # ── Step 1: assemble the parent bundle ───────────────────────────────────
    # Child apps already have macdeployqt run in their own CMakeLists.txt,
    # so they're self-contained with Qt frameworks when copied here.
    # rm -rf the destination before cp -a so that incremental rebuilds don't
    # create double-nesting (cp -a copies SOURCE *inside* DEST if DEST exists).
    # Use cp -a (not cmake -E copy_directory) to preserve framework symlinks.
    add_custom_command(TARGET deepn POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
            "$<TARGET_BUNDLE_DIR:deepn>/Contents/Data"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_SOURCE_DIR}/data/deepn.json"
            "$<TARGET_BUNDLE_DIR:deepn>/Contents/Data/"
        COMMAND ${CMAKE_COMMAND} -E rm -rf
            "$<TARGET_BUNDLE_DIR:deepn>/Contents/Resources/$<TARGET_FILE_BASE_NAME:gene_count>.app"
        COMMAND cp -a
            "$<TARGET_BUNDLE_DIR:gene_count>"
            "$<TARGET_BUNDLE_DIR:deepn>/Contents/Resources/$<TARGET_FILE_BASE_NAME:gene_count>.app"
        COMMAND ${CMAKE_COMMAND} -E rm -rf
            "$<TARGET_BUNDLE_DIR:deepn>/Contents/Resources/$<TARGET_FILE_BASE_NAME:junction_dice>.app"
        COMMAND cp -a
            "$<TARGET_BUNDLE_DIR:junction_dice>"
            "$<TARGET_BUNDLE_DIR:deepn>/Contents/Resources/$<TARGET_FILE_BASE_NAME:junction_dice>.app"
        COMMAND ${CMAKE_COMMAND} -E rm -rf
            "$<TARGET_BUNDLE_DIR:deepn>/Contents/Resources/$<TARGET_FILE_BASE_NAME:query>.app"
        COMMAND cp -a
            "$<TARGET_BUNDLE_DIR:query>"
            "$<TARGET_BUNDLE_DIR:deepn>/Contents/Resources/$<TARGET_FILE_BASE_NAME:query>.app"
        COMMAND ${CMAKE_COMMAND} -E rm -rf
            "$<TARGET_BUNDLE_DIR:deepn>/Contents/Resources/$<TARGET_FILE_BASE_NAME:read_depth>.app"
        COMMAND cp -a
            "$<TARGET_BUNDLE_DIR:read_depth>"
            "$<TARGET_BUNDLE_DIR:deepn>/Contents/Resources/$<TARGET_FILE_BASE_NAME:read_depth>.app"
        COMMAND ${CMAKE_COMMAND} -E rm -rf
            "$<TARGET_BUNDLE_DIR:deepn>/Contents/Resources/$<TARGET_FILE_BASE_NAME:deseq2_app>.app"
        COMMAND cp -a
            "$<TARGET_BUNDLE_DIR:deseq2_app>"
            "$<TARGET_BUNDLE_DIR:deepn>/Contents/Resources/$<TARGET_FILE_BASE_NAME:deseq2_app>.app"
        COMMENT "Bundling child apps and data into DEEPN++.app"
    )

    # ── Step 2: macdeployqt the parent for its own executable ────────────────
    if(MACDEPLOYQT_EXECUTABLE)
        add_custom_command(TARGET deepn POST_BUILD
            COMMAND ${MACDEPLOYQT_EXECUTABLE}
                "$<TARGET_BUNDLE_DIR:deepn>"
                -verbose=0 -codesign=-
            COMMENT "macdeployqt: embedding Qt frameworks into DEEPN++.app"
        )
    endif()

    # ── Step 3: code signing (Release builds only, when identity is provided) ─
    #
    # Set CODESIGN_IDENTITY at configure time to enable signing, e.g.:
    #   cmake -DCODESIGN_IDENTITY="Developer ID Application: Your Name (TEAMID)" ..
    #
    # Signing order enforced by cmake/codesign.cmake:
    #   dylibs → frameworks → nested .app bundles → DEEPN++.app
    # --options runtime + --timestamp satisfies Apple notarization requirements.
    set(CODESIGN_IDENTITY "" CACHE STRING
        "Code signing identity for release builds. \
Example: 'Developer ID Application: Your Name (TEAMID)'. \
Leave empty to skip signing.")

    if(CODESIGN_IDENTITY AND MACDEPLOYQT_EXECUTABLE)
        add_custom_command(TARGET deepn POST_BUILD
            COMMAND ${CMAKE_COMMAND}
                "-DAPP_PATH=$<TARGET_BUNDLE_DIR:deepn>"
                "-DIDENTITY=${CODESIGN_IDENTITY}"
                "-DENTITLEMENTS=${CMAKE_SOURCE_DIR}/cmake/entitlements.plist"
                -P "${CMAKE_SOURCE_DIR}/cmake/codesign.cmake"
            COMMENT "Code signing DEEPN++.app (identity: ${CODESIGN_IDENTITY})"
        )
    endif()
endif()
