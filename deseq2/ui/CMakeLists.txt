cmake_minimum_required(VERSION 3.16)
project(DESeq2GUI VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard to C++17 as required by workspace rules
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Force static linking for Qt5 but allow dynamic system libraries
set(BUILD_SHARED_LIBS OFF)

# Set Qt5 static installation path
set(QT5_STATIC_PATH "/Users/venky/Softwares/qt5-static")
set(CMAKE_PREFIX_PATH ${QT5_STATIC_PATH})

# Add library search paths for static Qt5
link_directories(${QT5_STATIC_PATH}/lib)

# Add QXlsx subdirectory if not already defined (gene_count may have added it)
if(NOT TARGET QXlsx)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../gene_count/QXlsx/QXlsx ${CMAKE_CURRENT_BINARY_DIR}/QXlsx)
endif()

# Find Qt5 components (static)
find_package(Qt5 REQUIRED COMPONENTS Core Widgets UiTools PrintSupport)

# Find Eigen3
find_package(Eigen3 REQUIRED)

# Enable automoc for Qt
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

# QCustomPlot sources
set(QCUSTOMPLOT_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../../QCustomplot/qcustomplot.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../../QCustomplot/qcustomplot.h
)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../statistics/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../QCustomplot)

# Source files following workspace rules for small classes
set(UI_SOURCES
    src/main.cpp
    src/main_window.cpp
    src/gene_count_handler.cpp
    include/main_window.h
    include/gene_count_handler.h
)

# Statistics library sources
set(STATISTICS_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../statistics/src/deseq_dataset.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../statistics/src/deseq_stats.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../statistics/src/utils.cpp
)

# UI file for Qt Designer
set(UI_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/mainwindow.ui)

# Create executable with macOS bundle
if(APPLE)
    add_executable(DESeq2 MACOSX_BUNDLE ${UI_SOURCES} ${STATISTICS_SOURCES} ${QCUSTOMPLOT_SOURCES} ${UI_FILES})
    
    # Set bundle properties
    set_target_properties(DESeq2 PROPERTIES
        MACOSX_BUNDLE_INFO_PLIST "${CMAKE_CURRENT_SOURCE_DIR}/Info.plist"
        MACOSX_BUNDLE_GUI_IDENTIFIER "com.deseq2.app"
        MACOSX_BUNDLE_BUNDLE_VERSION "${PROJECT_VERSION}"
        MACOSX_BUNDLE_SHORT_VERSION_STRING "${PROJECT_VERSION}"
        MACOSX_BUNDLE_NAME "DESeq2"
        MACOSX_BUNDLE_INFO_STRING "DESeq2 Differential Expression Analysis"
        MACOSX_BUNDLE_COPYRIGHT "Copyright © 2024"
        MACOSX_BUNDLE_ICON_FILE "DESeq2.icns"
        POSITION_INDEPENDENT_CODE OFF
    )
    
    # Copy icon to bundle
    add_custom_command(TARGET DESeq2 POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_CURRENT_SOURCE_DIR}/DESeq2.icns"
        "$<TARGET_FILE_DIR:DESeq2>/../Resources/DESeq2.icns"
    )
else()
    add_executable(DESeq2 ${UI_SOURCES} ${STATISTICS_SOURCES} ${QCUSTOMPLOT_SOURCES} ${UI_FILES})
    set_target_properties(DESeq2 PROPERTIES
        POSITION_INDEPENDENT_CODE OFF
    )
endif()

# Link Qt5 libraries, Eigen3, QXlsx, and PrintSupport (for QCustomPlot)
target_link_libraries(DESeq2
    Qt5::Core
    Qt5::Widgets
    Qt5::UiTools
    Qt5::PrintSupport
    Eigen3::Eigen
    QXlsx::QXlsx
)

# Set compiler flags following workspace rules
target_compile_options(DESeq2 PRIVATE -Wall -Wextra -O3)

# Set output directory
set_target_properties(DESeq2 PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Platform-specific linking flags
if(APPLE)
    # macOS specific linking flags
    target_link_options(DESeq2 PRIVATE
        -Wl,-dead_strip
        -Wl,-framework,CoreFoundation
        -Wl,-framework,ApplicationServices
        -Wl,-framework,CoreGraphics
        -Wl,-framework,CoreText
        -Wl,-framework,ImageIO
        -Wl,-framework,IOKit
        -Wl,-framework,OpenGL
        -Wl,-framework,AGL
        -Wl,-framework,Carbon
        -Wl,-framework,Cocoa
        -Wl,-framework,AudioToolbox
        -Wl,-framework,AudioUnit
        -Wl,-framework,CoreAudio
        -Wl,-framework,CoreVideo
        -Wl,-framework,QuartzCore
        -Wl,-framework,Security
        -Wl,-framework,SystemConfiguration
    )
elseif(UNIX AND NOT APPLE)
    # Linux specific linking flags
    target_link_options(DESeq2 PRIVATE
        -Wl,--whole-archive
        -Wl,--no-whole-archive
        -lm
        -lz
    )
elseif(WIN32)
    # Windows specific linking flags
    target_link_options(DESeq2 PRIVATE
        /SUBSYSTEM:WINDOWS
        /ENTRY:mainCRTStartup
    )
endif()

# Print configuration summary
message(STATUS "DESeq2GUI Configuration Summary:")
message(STATUS "  Executable: DESeq2")
message(STATUS "  Qt5 Static Path: ${QT5_STATIC_PATH}")
message(STATUS "  Output Directory: ${CMAKE_BINARY_DIR}/bin")
message(STATUS "  Platform: ${CMAKE_SYSTEM_NAME}") 