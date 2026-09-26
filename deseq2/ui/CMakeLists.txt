qt_add_executable(deseq2_app
    MACOSX_BUNDLE
    src/main.cpp
    src/main_window.cpp
    src/gene_count_handler.cpp
    include/main_window.h
    include/gene_count_handler.h
    src/mainwindow.ui
)

set_target_properties(deseq2_app PROPERTIES
    OUTPUT_NAME "StatMaker++"
    MACOSX_BUNDLE_BUNDLE_NAME "StatMaker++"
    MACOSX_BUNDLE_GUI_IDENTIFIER "edu.uiowa.deepn.statmaker"
    MACOSX_BUNDLE_BUNDLE_VERSION "${PROJECT_VERSION}"
    MACOSX_BUNDLE_SHORT_VERSION_STRING "${PROJECT_VERSION}"
)

if(APPLE AND MACDEPLOYQT_EXECUTABLE)
    add_custom_command(TARGET deseq2_app POST_BUILD
        COMMAND ${MACDEPLOYQT_EXECUTABLE} "$<TARGET_BUNDLE_DIR:deseq2_app>" -verbose=0 -codesign=-
        COMMENT "macdeployqt: StatMaker++.app"
    )
endif()

target_include_directories(deseq2_app PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(deseq2_app PRIVATE
    Qt6::Core
    Qt6::Gui
    Qt6::Widgets
    Qt6::Sql
    Qt6::Charts
    Qt6::PrintSupport
    Qt6::Svg
    deseq2_statistics
)
