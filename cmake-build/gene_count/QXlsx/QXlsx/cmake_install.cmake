# Install script for directory: /Users/venky/Projects/deepn-plus/gene_count/QXlsx/QXlsx

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "devel" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/Users/venky/Projects/deepn-plus/cmake-build/gene_count/QXlsx/QXlsx/libQXlsxQt5.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libQXlsxQt5.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libQXlsxQt5.a")
    execute_process(COMMAND "/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libQXlsxQt5.a")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "devel" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/QXlsxQt5" TYPE FILE FILES
    "/Users/venky/Projects/deepn-plus/gene_count/QXlsx/QXlsx/header/xlsxabstractooxmlfile.h"
    "/Users/venky/Projects/deepn-plus/gene_count/QXlsx/QXlsx/header/xlsxabstractsheet.h"
    "/Users/venky/Projects/deepn-plus/gene_count/QXlsx/QXlsx/header/xlsxabstractsheet_p.h"
    "/Users/venky/Projects/deepn-plus/gene_count/QXlsx/QXlsx/header/xlsxcellformula.h"
    "/Users/venky/Projects/deepn-plus/gene_count/QXlsx/QXlsx/header/xlsxcell.h"
    "/Users/venky/Projects/deepn-plus/gene_count/QXlsx/QXlsx/header/xlsxcelllocation.h"
    "/Users/venky/Projects/deepn-plus/gene_count/QXlsx/QXlsx/header/xlsxcellrange.h"
    "/Users/venky/Projects/deepn-plus/gene_count/QXlsx/QXlsx/header/xlsxcellreference.h"
    "/Users/venky/Projects/deepn-plus/gene_count/QXlsx/QXlsx/header/xlsxchart.h"
    "/Users/venky/Projects/deepn-plus/gene_count/QXlsx/QXlsx/header/xlsxchartsheet.h"
    "/Users/venky/Projects/deepn-plus/gene_count/QXlsx/QXlsx/header/xlsxconditionalformatting.h"
    "/Users/venky/Projects/deepn-plus/gene_count/QXlsx/QXlsx/header/xlsxdatavalidation.h"
    "/Users/venky/Projects/deepn-plus/gene_count/QXlsx/QXlsx/header/xlsxdatetype.h"
    "/Users/venky/Projects/deepn-plus/gene_count/QXlsx/QXlsx/header/xlsxdocument.h"
    "/Users/venky/Projects/deepn-plus/gene_count/QXlsx/QXlsx/header/xlsxformat.h"
    "/Users/venky/Projects/deepn-plus/gene_count/QXlsx/QXlsx/header/xlsxglobal.h"
    "/Users/venky/Projects/deepn-plus/gene_count/QXlsx/QXlsx/header/xlsxrichstring.h"
    "/Users/venky/Projects/deepn-plus/gene_count/QXlsx/QXlsx/header/xlsxworkbook.h"
    "/Users/venky/Projects/deepn-plus/gene_count/QXlsx/QXlsx/header/xlsxworksheet.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/venky/Projects/deepn-plus/cmake-build/gene_count/QXlsx/QXlsx/CMakeFiles/QXlsx.dir/install-cxx-module-bmi-noconfig.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "devel" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/QXlsxQt5/QXlsxQt5Targets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/QXlsxQt5/QXlsxQt5Targets.cmake"
         "/Users/venky/Projects/deepn-plus/cmake-build/gene_count/QXlsx/QXlsx/CMakeFiles/Export/9160ef171b5927dbe66bf41de9e1c9c5/QXlsxQt5Targets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/QXlsxQt5/QXlsxQt5Targets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/QXlsxQt5/QXlsxQt5Targets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/QXlsxQt5" TYPE FILE FILES "/Users/venky/Projects/deepn-plus/cmake-build/gene_count/QXlsx/QXlsx/CMakeFiles/Export/9160ef171b5927dbe66bf41de9e1c9c5/QXlsxQt5Targets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^()$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/QXlsxQt5" TYPE FILE FILES "/Users/venky/Projects/deepn-plus/cmake-build/gene_count/QXlsx/QXlsx/CMakeFiles/Export/9160ef171b5927dbe66bf41de9e1c9c5/QXlsxQt5Targets-noconfig.cmake")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/QXlsxQt5" TYPE FILE FILES
    "/Users/venky/Projects/deepn-plus/cmake-build/gene_count/QXlsx/QXlsx/QXlsxQt5Config.cmake"
    "/Users/venky/Projects/deepn-plus/cmake-build/gene_count/QXlsx/QXlsx/QXlsxQt5ConfigVersion.cmake"
    )
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/venky/Projects/deepn-plus/cmake-build/gene_count/QXlsx/QXlsx/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
