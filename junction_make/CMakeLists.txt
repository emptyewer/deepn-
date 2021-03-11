cmake_minimum_required(VERSION 3.5)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(JunctionMake++ LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15")
set(CONAN_DISABLE_CHECK_COMPILER TRUE)

# add conan
include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

find_package(QT NAMES Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt${QT_VERSION_MAJOR} COMPONENTS Widgets REQUIRED)

set(PROJECT_SOURCES
        main.cpp
        mainwindow.cpp
        mainwindow.h
        mainwindow.ui
)

set(EXECUTABLE "${PROJECT_NAME}")

# OS Specific configurations
if (APPLE)
    set(EXECUTABLE "${PROJECT_NAME}")
    set(MACOSX_BUNDLE_ICON_FILE "junction_make.icns")
    set(JUNCTION_MAKE_ICON ${CMAKE_CURRENT_SOURCE_DIR}/../icons/junction_make.icns)
    set_source_files_properties(${JUNCTION_MAKE_ICON} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")
    add_executable(${EXECUTABLE} MACOSX_BUNDLE ${PROJECT_SOURCES} ${JUNCTION_MAKE_ICON})
elseif (WIN32)
    set(JUNCTION_MAKE_WINRC ../icons/junction_make_windows.rc)
    set_property(SOURCE ${PROJECT_SOURCE_DIR}/../icons/junction_make_windows.rc APPEND PROPERTY OBJECT_DEPENDS ${PROJECT_SOURCE_DIR}/../icons/junction_make.ico)
    add_executable(${EXECUTABLE} WIN32 ${PROJECT_SOURCES} ${JUNCTION_MAKE_WINRC})
else()
    add_executable(${EXECUTABLE} ${PROJECT_SOURCES})
endif()

target_link_libraries(${EXECUTABLE} PRIVATE Qt${QT_VERSION_MAJOR}::Widgets ${CONAN_LIBS_RELEASE})


