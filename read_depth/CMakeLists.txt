set(ICON_FILE ${CMAKE_CURRENT_SOURCE_DIR}/../icons/read_depth.icns)
set_source_files_properties(${ICON_FILE} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")

qt_add_executable(read_depth MACOSX_BUNDLE
    main.cpp
    mainwindow.cpp mainwindow.h
    depth_plot_widget.cpp depth_plot_widget.h
    boundary_detector.cpp boundary_detector.h
    mainwindow.ui
    ${ICON_FILE}
)

set_target_properties(read_depth PROPERTIES
    OUTPUT_NAME "ReadDepth++"
    MACOSX_BUNDLE_ICON_FILE read_depth.icns
    MACOSX_BUNDLE_BUNDLE_NAME "ReadDepth++"
    MACOSX_BUNDLE_GUI_IDENTIFIER "edu.uiowa.deepn.readdepth"
    MACOSX_BUNDLE_BUNDLE_VERSION "${PROJECT_VERSION}"
    MACOSX_BUNDLE_SHORT_VERSION_STRING "${PROJECT_VERSION}"
)

if(APPLE AND MACDEPLOYQT_EXECUTABLE)
    add_custom_command(TARGET read_depth POST_BUILD
        COMMAND ${MACDEPLOYQT_EXECUTABLE} "$<TARGET_BUNDLE_DIR:read_depth>" -verbose=0 -codesign=-
        COMMENT "macdeployqt: ReadDepth++.app"
    )
endif()

if(WIN32)
    target_sources(read_depth PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../icons/read_depth_windows.rc)
endif()

target_link_libraries(read_depth PRIVATE
    deepn_data
    deepn_widgets
    Qt6::Core Qt6::Gui Qt6::Widgets Qt6::Sql Qt6::Charts Qt6::Svg
)
