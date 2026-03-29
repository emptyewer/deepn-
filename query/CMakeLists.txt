set(ICON_FILE ${CMAKE_CURRENT_SOURCE_DIR}/../icons/multi_query.icns)
set_source_files_properties(${ICON_FILE} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")

qt_add_executable(query
    MACOSX_BUNDLE
    main.cpp
    mainwindow.cpp
    mainwindow.h
    dockwidget.cpp
    dockwidget.h
    mainwindow.ui
    dockwidget.ui
    ${ICON_FILE}
)

set_target_properties(query PROPERTIES
    OUTPUT_NAME "MultiQuery++"
    MACOSX_BUNDLE_ICON_FILE multi_query.icns
)

if(WIN32)
    target_sources(query PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../icons/multi_query_windows.rc)
endif()

target_link_libraries(query PRIVATE
    Qt6::Core
    Qt6::Gui
    Qt6::Widgets
    Qt6::Charts
)
