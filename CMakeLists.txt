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


if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
  message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
  file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/v0.16.1/conan.cmake"
                "${CMAKE_BINARY_DIR}/conan.cmake"
                TLS_VERIFY ON)
endif()

include(${CMAKE_BINARY_DIR}/conan.cmake)
conan_cmake_configure(REQUIRES fmt/6.1.2 boost/1.75.0 spdlog/1.8.2 nlohmann_json/3.9.1 yaml-cpp/0.6.3 seqan/2.4.0 GENERATORS cmake_find_package)
conan_cmake_autodetect(settings)
conan_cmake_install(PATH_OR_REFERENCE .
                    BUILD missing
                    REMOTE conan-center
                    SETTINGS ${settings})

# Libraries
find_package(QT NAMES Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt${QT_VERSION_MAJOR} COMPONENTS Widgets REQUIRED)
find_program(7ZIP NAMES 7z 7zip REQUIRED)

add_subdirectory(gene_count)
add_subdirectory(junction_make)

# OS Specific configurations
if (APPLE)
    set(EXECUTABLE "${PROJECT_NAME}")
    set(MACOSX_BUNDLE_ICON_FILE "deepn.icns")
    set(DEEPN_ICON ${CMAKE_CURRENT_SOURCE_DIR}/icons/deepn.icns)
    set_source_files_properties(${DEEPN_ICON} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")
    add_executable(${EXECUTABLE} MACOSX_BUNDLE ${PROJECT_SOURCES} ${DEEPN_ICON})
    set(EXECUTABLE_PATH "${${EXECUTABLE}_BINARY_DIR}/${EXECUTABLE}.app")
elseif (WIN32)
    set(DEEPN_WINRC icons/deepn_windows.rc)
    set_property(SOURCE ${PROJECT_SOURCE_DIR}/icons/deepn_windows.rc APPEND PROPERTY OBJECT_DEPENDS ${PROJECT_SOURCE_DIR}/icons/deepn.ico)
    add_executable(${EXECUTABLE} WIN32 ${PROJECT_SOURCES} ${DEEPN_WINRC})
else()
    add_executable(${EXECUTABLE} ${PROJECT_SOURCES})
endif()

target_link_libraries(${EXECUTABLE} PRIVATE Qt${QT_VERSION_MAJOR}::Widgets ${CONAN_LIBS_RELEASE})

if (APPLE)
    add_custom_command(TARGET ${EXECUTABLE} PRE_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${EXECUTABLE_PATH}/Contents/Apps)
    add_custom_command(TARGET ${EXECUTABLE} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_BINARY_DIR}/gene_count/GeneCount++.app ${EXECUTABLE_PATH}/Contents/Apps/GeneCount++.app VERBATIM)
    add_custom_command(TARGET ${EXECUTABLE} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_BINARY_DIR}/junction_make/JunctionMake++.app ${EXECUTABLE_PATH}/Contents/Apps/JunctionMake++.app VERBATIM)
    add_custom_command(TARGET ${EXECUTABLE} POST_BUILD COMMAND ${7ZIP} a -m0=ppmd -mx=9 ${CMAKE_CURRENT_SOURCE_DIR}/deepn++_mac_64bit.7z ${EXECUTABLE_PATH} VERBATIM)
endif()
