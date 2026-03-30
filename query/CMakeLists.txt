set(ICON_FILE ${CMAKE_CURRENT_SOURCE_DIR}/../icons/multi_query.icns)
set_source_files_properties(${ICON_FILE} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")

qt_add_executable(query
    MACOSX_BUNDLE
    main.cpp
    mainwindow.cpp
    mainwindow.h
    junction_plot_widget.cpp
    junction_plot_widget.h
    junction_table_model.cpp
    junction_table_model.h
    gene_detail_panel.cpp
    gene_detail_panel.h
    comparison_manager.cpp
    comparison_manager.h
    mainwindow.ui
    ${ICON_FILE}
)

set_target_properties(query PROPERTIES
    OUTPUT_NAME "MultiQuery++"
    MACOSX_BUNDLE_ICON_FILE multi_query.icns
    MACOSX_BUNDLE_BUNDLE_NAME "MultiQuery++"
    MACOSX_BUNDLE_GUI_IDENTIFIER "edu.uiowa.deepn.multiquery"
)

if(WIN32)
    target_sources(query PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../icons/multi_query_windows.rc)
endif()

target_link_libraries(query PRIVATE
    deepn_data
    deepn_widgets
    Qt6::Core
    Qt6::Gui
    Qt6::Widgets
    Qt6::Sql
    Qt6::Charts
    Qt6::Svg
)
