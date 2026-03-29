set(ICON_FILE ${CMAKE_CURRENT_SOURCE_DIR}/../icons/read_depth.icns)
set_source_files_properties(${ICON_FILE} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")

qt_add_executable(read_depth
    MACOSX_BUNDLE
    main.cpp
    mainwindow.cpp
    mainwindow.h
    mainwindow.ui
    ${ICON_FILE}
)

set_target_properties(read_depth PROPERTIES
    OUTPUT_NAME "ReadDepth++"
    MACOSX_BUNDLE_ICON_FILE read_depth.icns
)

if(WIN32)
    target_sources(read_depth PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../icons/read_depth_windows.rc)
endif()

target_link_libraries(read_depth PRIVATE
    Qt6::Core
    Qt6::Gui
    Qt6::Widgets
    Qt6::Charts
)
