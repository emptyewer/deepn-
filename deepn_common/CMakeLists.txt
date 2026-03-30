# deepn_common: shared data layer and widgets for MultiQuery++ and ReadDepth++

# --- Data layer library (no Qt Widgets dependency) ---
add_library(deepn_data STATIC
    include/data_structures.h
    include/gene_annotation_db.h
    include/sqlite_junction_loader.h
    include/depth_calculator.h
    include/export_engine.h
    include/batch_runner.h
    src/gene_annotation_db.cpp
    src/sqlite_junction_loader.cpp
    src/depth_calculator.cpp
    src/export_engine.cpp
    src/batch_runner.cpp
)

target_include_directories(deepn_data PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(deepn_data PUBLIC
    Qt6::Core
    Qt6::Gui
    Qt6::Widgets
    Qt6::Sql
    Qt6::Svg
    Qt6::PrintSupport
)

set_target_properties(deepn_data PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    AUTOMOC ON
)

# --- Widget library (Qt Widgets + Charts dependency) ---
add_library(deepn_widgets STATIC
    include/gene_selector_widget.h
    include/mrna_track_widget.h
    include/sync_controller.h
    src/gene_selector_widget.cpp
    src/mrna_track_widget.cpp
    src/sync_controller.cpp
)

target_include_directories(deepn_widgets PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(deepn_widgets PUBLIC
    deepn_data
    Qt6::Widgets
    Qt6::Charts
)

set_target_properties(deepn_widgets PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    AUTOMOC ON
)
