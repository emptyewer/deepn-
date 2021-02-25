cmake_minimum_required(VERSION 3.5)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(DEEPN++ LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15")
set(CONAN_DISABLE_CHECK_COMPILER TRUE)

# Sources
set(PROJECT_SOURCES
        main.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/mainwindow.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/mainwindow.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/mainwindow.ui
)

set(EXECUTABLE "${PROJECT_NAME}")

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

# Libraries
find_package(QT NAMES Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt${QT_VERSION_MAJOR} COMPONENTS Widgets REQUIRED)
find_program(7ZIP NAMES 7z 7zip REQUIRED)

# OS Specific configurations
if (APPLE)
    set(EXECUTABLE "${PROJECT_NAME}")
    set(MACOSX_BUNDLE_ICON_FILE "deepn.icns")
    set(DEEPN_ICON ${CMAKE_CURRENT_SOURCE_DIR}/icons/deepn.icns)
    set_source_files_properties(${DEEPN_ICON} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")
    add_executable(${EXECUTABLE} MACOSX_BUNDLE ${PROJECT_SOURCES} ${DEEPN_ICON})
    get_filename_component(QT_BIN_DIRECTORY ${QT_QMAKE_EXECUTABLE} DIRECTORY)
    set(DEPLOY_EXECUTABLE "${QT_BIN_DIRECTORY}/macdeployqt")
    set(EXECUTABLE_PATH "${${EXECUTABLE}_BINARY_DIR}/${EXECUTABLE}.app")
elseif (WIN32)
    set(DEEPN_WINRC deepn_windows.rc)
    set_property(SOURCE deepn_windows.rc APPEND PROPERTY OBJECT_DEPENDS ${PROJECT_SOURCE_DIR}/icons/deepn.ico)
    add_executable(${EXECUTABLE} WIN32 ${PROJECT_SOURCES} ${DEEPN_WINRC})
else()
    add_executable(${EXECUTABLE} ${PROJECT_SOURCES})
endif()

target_link_libraries(${EXECUTABLE} PUBLIC Qt${QT_VERSION_MAJOR}::Widgets ${CONAN_LIBS_RELEASE})

if (APPLE)
    add_custom_command(TARGET ${EXECUTABLE} POST_BUILD COMMAND ${DEPLOY_EXECUTABLE} ${EXECUTABLE_PATH} VERBATIM)
    add_custom_command(TARGET ${EXECUTABLE} POST_BUILD COMMAND ${7ZIP} a -m0=ppmd -mx=9 ${CMAKE_CURRENT_SOURCE_DIR}/deepn++_mac_64bit.7z ${EXECUTABLE_PATH} VERBATIM)
endif()
