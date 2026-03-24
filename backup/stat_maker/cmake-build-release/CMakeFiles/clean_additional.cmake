# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "CMakeFiles/deepn_cpp_qt_autogen.dir/AutogenUsed.txt"
  "CMakeFiles/deepn_cpp_qt_autogen.dir/ParseCache.txt"
  "deepn_cpp_qt_autogen"
  )
endif()
