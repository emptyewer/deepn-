cmake_minimum_required(VERSION 3.21)
project(deepn++ VERSION 1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

if(APPLE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "14.0" CACHE STRING "Minimum macOS deployment target")
endif()

# --- Qt6 ---
find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets Sql Charts PrintSupport Svg)

# --- Eigen3 (header-only, for DESeq2 statistics) ---
find_package(Eigen3 REQUIRED)

# --- zlib (static, from submodule) ---
set(ZLIB_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
add_subdirectory(zlib EXCLUDE_FROM_ALL)

# --- QSimpleXlsxWriter (from xlsx submodule) ---
set(XLSX_SRC ${CMAKE_CURRENT_SOURCE_DIR}/xlsx/simplexlsx-code)

add_library(QSimpleXlsxWriter STATIC
    ${XLSX_SRC}/PathManager.cpp
    ${XLSX_SRC}/Xlsx/Chart.cpp
    ${XLSX_SRC}/Xlsx/Chartsheet.cpp
    ${XLSX_SRC}/Xlsx/Drawing.cpp
    ${XLSX_SRC}/Xlsx/SimpleXlsxDef.cpp
    ${XLSX_SRC}/Xlsx/Workbook.cpp
    ${XLSX_SRC}/Xlsx/Worksheet.cpp
    ${XLSX_SRC}/Xlsx/XlsxHeaders.cpp
    ${XLSX_SRC}/XLSXColors/clsRGBColorRecord.cpp
    ${XLSX_SRC}/XLSXColors/XLSXColorLib.cpp
    ${XLSX_SRC}/Zip/zip.cpp
)

target_include_directories(QSimpleXlsxWriter PUBLIC
    ${XLSX_SRC}
    ${XLSX_SRC}/Xlsx
    ${XLSX_SRC}/XLSXColors
    ${XLSX_SRC}/Zip
)

target_compile_features(QSimpleXlsxWriter PRIVATE cxx_std_11)

# --- NCBI BLAST (blastn from submodule, autotools build) ---
include(ExternalProject)
include(ProcessorCount)
ProcessorCount(NPROC)
if(NPROC EQUAL 0)
    set(NPROC 4)
endif()

set(NCBI_BUILD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ncbi/build/build)

ExternalProject_Add(ncbi_blast
    SOURCE_DIR        ${CMAKE_CURRENT_SOURCE_DIR}/ncbi
    CONFIGURE_COMMAND bash -c "test -f '${NCBI_BUILD_DIR}/Makefile' || '${CMAKE_CURRENT_SOURCE_DIR}/ncbi/configure' --with-static --with-mt --with-64 --with-static-exe --with-build-root='${CMAKE_CURRENT_SOURCE_DIR}/ncbi/build' --without-boost --without-gui"
    BUILD_COMMAND     bash -c "make -C '${NCBI_BUILD_DIR}' all_r -k -j${NPROC} || true"
    BUILD_IN_SOURCE   TRUE
    INSTALL_COMMAND   ""
    BUILD_BYPRODUCTS  ${CMAKE_CURRENT_SOURCE_DIR}/ncbi/build/bin/blastn
)

# --- blat (from submodule, plain Makefile) ---
ExternalProject_Add(blat_build
    SOURCE_DIR        ${CMAKE_CURRENT_SOURCE_DIR}/blat
    CONFIGURE_COMMAND ""
    BUILD_COMMAND     ${CMAKE_CURRENT_SOURCE_DIR}/cmake/build_blat.sh
        ${CMAKE_CURRENT_SOURCE_DIR}/blat ${NPROC}
    BUILD_IN_SOURCE   TRUE
    INSTALL_COMMAND   ""
    BUILD_BYPRODUCTS  ${CMAKE_CURRENT_SOURCE_DIR}/blat/bin/blat
)

# --- DESeq2 statistics library (header-only Eigen3, no Qt) ---
set(DESEQ2_STATS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/deseq2/statistics)
add_library(deseq2_statistics STATIC
    ${DESEQ2_STATS_DIR}/src/deseq_dataset.cpp
    ${DESEQ2_STATS_DIR}/src/deseq_stats.cpp
    ${DESEQ2_STATS_DIR}/src/utils.cpp
    ${DESEQ2_STATS_DIR}/src/y2h_scores.cpp
)
target_include_directories(deseq2_statistics PUBLIC
    ${DESEQ2_STATS_DIR}/include
)
target_link_libraries(deseq2_statistics PUBLIC Eigen3::Eigen)
set_target_properties(deseq2_statistics PROPERTIES POSITION_INDEPENDENT_CODE ON)

# DESeq2 validation test
add_executable(test_deseq2_validation ${DESEQ2_STATS_DIR}/tests/test_validation.cpp)
target_link_libraries(test_deseq2_validation deseq2_statistics)

# Y2H-SCORES validation test
add_executable(test_y2h_scores ${DESEQ2_STATS_DIR}/tests/test_y2h_scores.cpp)
target_link_libraries(test_y2h_scores deseq2_statistics)

# --- Subprojects ---
add_subdirectory(junction_dice)
add_subdirectory(gene_count)
add_subdirectory(query)
add_subdirectory(read_depth)
add_subdirectory(deseq2/ui)
add_subdirectory(deepn)
