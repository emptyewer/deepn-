cmake_minimum_required(VERSION 3.5)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(DEEPN++ LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15")
set(CONAN_DISABLE_CHECK_COMPILER TRUE)

# macros
macro(print_all_variables)
    message(STATUS "print_all_variables------------------------------------------{")
    get_cmake_property(_variableNames VARIABLES)
    foreach (_variableName ${_variableNames})
        message(STATUS "${_variableName}=${${_variableName}}")
    endforeach()
    message(STATUS "print_all_variables------------------------------------------}")
endmacro()

# add conan
include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

# Sources
set(PROJECT_SOURCES
        main.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/deepn/mainwindow.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/deepn/mainwindow.h
        ${CMAKE_CURRENT_SOURCE_DIR}/deepn/dqtextedit.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/deepn/dqtextedit.h
        ${CMAKE_CURRENT_SOURCE_DIR}/deepn/mainwindow.ui
)
set(EXECUTABLE "${PROJECT_NAME}")

# Libraries
find_package(QT NAMES Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt${QT_VERSION_MAJOR} COMPONENTS Widgets REQUIRED)
find_program(7ZIP NAMES 7z 7zip REQUIRED)

# Subdirectories
add_subdirectory(gene_count)
add_subdirectory(junction_make)

# OS Specific configurations
if (APPLE)
    set(MACOSX_BUNDLE_ICON_FILE "deepn.icns")
    set(DEEPN_ICON ${CMAKE_CURRENT_SOURCE_DIR}/icons/deepn.icns)
    set_source_files_properties(${DEEPN_ICON} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")
    add_executable(${EXECUTABLE} MACOSX_BUNDLE ${PROJECT_SOURCES} ${DEEPN_ICON})
elseif (WIN32)
    set(DEEPN_WINRC icons/deepn_windows.rc)
    set_property(SOURCE ${PROJECT_SOURCE_DIR}/icons/deepn_windows.rc APPEND PROPERTY OBJECT_DEPENDS ${PROJECT_SOURCE_DIR}/icons/deepn.ico)
    add_executable(${EXECUTABLE} WIN32 ${PROJECT_SOURCES} ${DEEPN_WINRC})
else()
    add_executable(${EXECUTABLE} ${PROJECT_SOURCES})
endif()

target_link_libraries(${EXECUTABLE} PRIVATE Qt${QT_VERSION_MAJOR}::Widgets ${CONAN_LIBS_RELEASE})
print_all_variables()
if (APPLE)
    add_custom_command(TARGET ${EXECUTABLE} PRE_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin/${EXECUTABLE}.app/Contents/Apps)
    add_custom_command(TARGET ${EXECUTABLE} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory ${GeneCount++_BINARY_DIR}/bin/GeneCount++.app ${CMAKE_BINARY_DIR}/bin/${EXECUTABLE}.app/Contents/Apps/GeneCount++.app VERBATIM)
    add_custom_command(TARGET ${EXECUTABLE} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory ${JunctionMake++_BINARY_DIR}/bin/JunctionMake++.app ${CMAKE_BINARY_DIR}/bin/${EXECUTABLE}.app/Contents/Apps/JunctionMake++.app VERBATIM)
    add_custom_command(TARGET ${EXECUTABLE} POST_BUILD COMMAND ${7ZIP} a -m0=ppmd -mx=9 ${CMAKE_CURRENT_SOURCE_DIR}/deepn++_mac_64bit.7z ${EXECUTABLE_PATH} VERBATIM)
endif()
