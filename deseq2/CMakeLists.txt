cmake_minimum_required(VERSION 3.16)

# Set CMake policies
cmake_policy(SET CMP0167 NEW)

project(DESeq2 VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard to C++17 as required by workspace rules
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Force static linking for Qt5 but allow dynamic system libraries
set(BUILD_SHARED_LIBS OFF)

# Set Qt5 static installation path
set(QT5_STATIC_PATH "/Users/venky/Softwares/qt5-static")
list(APPEND CMAKE_PREFIX_PATH ${QT5_STATIC_PATH})

# Add library search paths for static Qt5
link_directories(${QT5_STATIC_PATH}/lib)

# Set static linking for Qt5
set(QT_USE_STATIC ON)
set(QT_USE_STATIC_RUNTIME ON)

# Find required packages for statistics library
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem)

# Find Qt5 components for UI (static)
find_package(Qt5 REQUIRED COMPONENTS Core Widgets UiTools)

# Enable automoc for Qt
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

# Set global compiler flags following workspace rules
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3")

# Platform-specific settings
if(APPLE)
    # macOS specific settings - clang++ doesn't support -static-libgcc
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15")
elseif(UNIX AND NOT APPLE)
    # Linux specific settings
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
elseif(WIN32)
    # Windows specific settings
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
endif()

# Add subdirectories in dependency order
add_subdirectory(statistics)
add_subdirectory(ui)

# Print configuration summary
message(STATUS "DESeq2 Configuration Summary:")
message(STATUS "  Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Qt5 Static Path: ${QT5_STATIC_PATH}")
message(STATUS "  Eigen3 Found: ${EIGEN3_FOUND}")
message(STATUS "  Boost Found: ${Boost_FOUND}")
message(STATUS "  Qt5 Components: Core=${Qt5Core_FOUND}, Widgets=${Qt5Widgets_FOUND}, UiTools=${Qt5UiTools_FOUND}")

# Create a custom target for cleaning all build artifacts
add_custom_target(clean-all
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_SOURCE_DIR}/build
    COMMENT "Cleaning all build artifacts"
)

 