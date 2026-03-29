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
)

if(WIN32)
    target_sources(deepn PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../icons/deepn_windows.rc)
endif()

target_include_directories(deepn PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(deepn PRIVATE
    Qt6::Core
    Qt6::Gui
    Qt6::Widgets
)

add_dependencies(deepn junction_dice gene_count query read_depth deseq2_app)

if(APPLE)
    add_custom_command(TARGET deepn POST_BUILD
        # Data for main app
        COMMAND ${CMAKE_COMMAND} -E make_directory
            "$<TARGET_BUNDLE_DIR:deepn>/Contents/Data"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_SOURCE_DIR}/data/deepn.json"
            "$<TARGET_BUNDLE_DIR:deepn>/Contents/Data/"
        # Bundle child apps into Resources
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            "$<TARGET_BUNDLE_DIR:gene_count>"
            "$<TARGET_BUNDLE_DIR:deepn>/Contents/Resources/$<TARGET_FILE_BASE_NAME:gene_count>.app"
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            "$<TARGET_BUNDLE_DIR:junction_dice>"
            "$<TARGET_BUNDLE_DIR:deepn>/Contents/Resources/$<TARGET_FILE_BASE_NAME:junction_dice>.app"
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            "$<TARGET_BUNDLE_DIR:query>"
            "$<TARGET_BUNDLE_DIR:deepn>/Contents/Resources/$<TARGET_FILE_BASE_NAME:query>.app"
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            "$<TARGET_BUNDLE_DIR:read_depth>"
            "$<TARGET_BUNDLE_DIR:deepn>/Contents/Resources/$<TARGET_FILE_BASE_NAME:read_depth>.app"
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            "$<TARGET_BUNDLE_DIR:deseq2_app>"
            "$<TARGET_BUNDLE_DIR:deepn>/Contents/Resources/$<TARGET_FILE_BASE_NAME:deseq2_app>.app"
        COMMENT "Bundling child apps and data into DEEPN++.app"
    )
endif()
