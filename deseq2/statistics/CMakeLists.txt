cmake_minimum_required(VERSION 3.16)

# Set CMake policies
cmake_policy(SET CMP0167 NEW)

project(DESeq2Statistics VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard to C++17 as required by workspace rules
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Force static linking for this library
set(BUILD_SHARED_LIBS OFF)

# Find required packages
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${EIGEN3_INCLUDE_DIR})

# Source files following workspace rules for small classes
set(STATISTICS_SOURCES
    src/deseq_dataset.cpp
    src/deseq_stats.cpp
    src/utils.cpp
)

# Create static library
add_library(deseq2_statistics STATIC ${STATISTICS_SOURCES})

# Set library properties
set_target_properties(deseq2_statistics PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/include/deseq_dataset.h;${CMAKE_CURRENT_SOURCE_DIR}/include/deseq_stats.h;${CMAKE_CURRENT_SOURCE_DIR}/include/utils.h"
    POSITION_INDEPENDENT_CODE ON
)

# Link libraries (static)
target_link_libraries(deseq2_statistics
    Eigen3::Eigen
    Boost::filesystem
)

# Set compiler flags following workspace rules
target_compile_options(deseq2_statistics PRIVATE -Wall -Wextra -O3)

# Create test executables
add_executable(test_deseq2 tests/test_deseq2.cpp)
target_link_libraries(test_deseq2 deseq2_statistics)

add_executable(example_usage tests/example_usage.cpp)
target_link_libraries(example_usage deseq2_statistics)

add_executable(test_ppm_conversion tests/test_ppm_conversion.cpp)
target_link_libraries(test_ppm_conversion deseq2_statistics)

add_executable(example_ppm_deseq2_analysis tests/example_ppm_deseq2_analysis.cpp)
target_link_libraries(example_ppm_deseq2_analysis deseq2_statistics)

# Set compiler flags for all executables
set(EXECUTABLES test_deseq2 example_usage test_ppm_conversion example_ppm_deseq2_analysis)
foreach(EXEC ${EXECUTABLES})
    target_compile_options(${EXEC} PRIVATE -Wall -Wextra -O3)
endforeach()

# Print configuration summary
message(STATUS "DESeq2Statistics Configuration Summary:")
message(STATUS "  Library: deseq2_statistics (static)")
message(STATUS "  Test Executables: ${EXECUTABLES}")
message(STATUS "  Eigen3 Version: ${EIGEN3_VERSION}")
message(STATUS "  Boost Version: ${Boost_VERSION}") 