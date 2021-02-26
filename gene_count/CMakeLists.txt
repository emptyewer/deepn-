cmake_minimum_required(VERSION 3.5)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(GeneCount++ LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15")
set(CONAN_DISABLE_CHECK_COMPILER TRUE)

find_package(QT NAMES Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt${QT_VERSION_MAJOR} COMPONENTS Widgets REQUIRED)

set(PROJECT_SOURCES
        main.cpp
        mainwindow.cpp
        mainwindow.h
        mainwindow.ui
)

set(EXECUTABLE "${PROJECT_NAME}")

include(${CMAKE_BINARY_DIR}/conan.cmake)
conan_cmake_configure(REQUIRES fmt/6.1.2 GENERATORS cmake_find_package)
conan_cmake_autodetect(settings)
conan_cmake_install(PATH_OR_REFERENCE .
                    BUILD missing
                    REMOTE conan-center
                    SETTINGS ${settings})

# OS Specific configurations
if (APPLE)
    set(EXECUTABLE "${PROJECT_NAME}")
    set(MACOSX_BUNDLE_ICON_FILE "gene_count.icns")
    set(GENE_COUNT_ICON ${CMAKE_CURRENT_SOURCE_DIR}/../icons/gene_count.icns)
    set_source_files_properties(${GENE_COUNT_ICON} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")
    add_executable(${EXECUTABLE} MACOSX_BUNDLE ${PROJECT_SOURCES} ${GENE_COUNT_ICON})
    get_filename_component(QT_BIN_DIRECTORY ${QT_QMAKE_EXECUTABLE} DIRECTORY)
    set(DEPLOY_EXECUTABLE "${QT_BIN_DIRECTORY}/macdeployqt")
    set(EXECUTABLE_PATH "${${EXECUTABLE}_BINARY_DIR}/${EXECUTABLE}.app")
elseif (WIN32)
    set(GENE_COUNT_WINRC ../icons/gene_count_windows.rc)
    set_property(SOURCE ${PROJECT_SOURCE_DIR}/../icons/gene_count_windows.rc APPEND PROPERTY OBJECT_DEPENDS ${PROJECT_SOURCE_DIR}/../icons/gene_count.ico)
    add_executable(${EXECUTABLE} WIN32 ${PROJECT_SOURCES} ${GENE_COUNT_WINRC})
else()
    add_executable(${EXECUTABLE} ${PROJECT_SOURCES})
endif()

target_link_libraries(${EXECUTABLE} PRIVATE Qt${QT_VERSION_MAJOR}::Widgets ${CONAN_LIBS_RELEASE})

if (APPLE)
    add_custom_command(TARGET ${EXECUTABLE} POST_BUILD COMMAND ${DEPLOY_EXECUTABLE} ${EXECUTABLE_PATH} VERBATIM)
endif()

