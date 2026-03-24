cmake_minimum_required(VERSION 3.5...3.20)
project(StatMaker++)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOUIC_SEARCH_PATHS .)

set(CMAKE_PREFIX_PATH
        "$ENV{HOME}/Softwares/qt5-static/lib/cmake"
        ${CMAKE_PREFIX_PATH}
)

find_package(Qt5 COMPONENTS Core Widgets Concurrent REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
        main.cpp
        mainwindow.cpp
        statsmake.cpp
        customlistview.cpp
        signals.cpp
        datastructs.cpp
        AdvancedNBDispersion.cpp
        AdvancedNBDispersion.h
)

set(HEADERS
        mainwindow.h
        statsmake.h
        customlistview.h
        signals.h
        datastructs.h
)

set(FORMS
        mainwindow.ui
)

if (APPLE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET 14.0)
    set(MACOSX_BUNDLE_ICON_FILE "deepn.icns")
    set(APP_ICON_MACOS "${CMAKE_CURRENT_SOURCE_DIR}/../icons/deepn.icns")
    set_source_files_properties(${APP_ICON_MACOS} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")
    add_executable(StatMaker++ MACOSX_BUNDLE ${SOURCES} ${HEADERS} ${FORMS} ${APP_ICON_MACOS})
    #    add_custom_command(TARGET StatMaker++ POST_BUILD
    #            COMMAND macdeployqt $<TARGET_BUNDLE_DIR:StatMaker++>
    #            COMMAND /bin/cp ${CMAKE_CURRENT_SOURCE_DIR}/../build/JunctionDice-Info.plist $<TARGET_BUNDLE_DIR:JunctionDice++>/Contents/Info.plist
    #    )
elseif (WIN32)
    set(APP_ICON_WIN "${CMAKE_CURRENT_SOURCE_DIR}/../icons/deepn.ico")
    add_executable(StatMaker++ ${SOURCES} ${HEADERS} ${FORMS} ${APP_ICON_WIN})
endif ()


if (APPLE)
    set_target_properties(StatMaker++ PROPERTIES
            MACOSX_BUNDLE TRUE
            OUTPUT_NAME StatMaker++
            RUNTIME_OUTPUT_DIRECTORY ../components
    )
elseif (WIN32)
    set_target_properties(StatMaker++ PROPERTIES
            WIN32_EXECUTABLE TRUE
            OUTPUT_NAME StatMaker++
            RUNTIME_OUTPUT_DIRECTORY ../components
    )
endif ()


target_link_libraries(StatMaker++ Qt5::Core Qt5::Widgets Qt5::Concurrent
        "-framework Security")
