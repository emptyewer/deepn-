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
)

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
