# fix_child_rpaths.cmake
# Fixes @rpath in child apps nested inside DEEPN++.app/Contents/Resources/
# so they find Qt frameworks in the parent bundle's Contents/Frameworks/.
#
# Usage: cmake -DAPP_BUNDLE=/path/to/DEEPN++.app -P fix_child_rpaths.cmake

if(NOT APP_BUNDLE)
    message(FATAL_ERROR "APP_BUNDLE not set")
endif()

set(RESOURCES_DIR "${APP_BUNDLE}/Contents/Resources")

# Find all child app executables
file(GLOB CHILD_APPS "${RESOURCES_DIR}/*.app")

foreach(CHILD_APP ${CHILD_APPS})
    get_filename_component(APP_NAME "${CHILD_APP}" NAME_WE)
    set(EXECUTABLE "${CHILD_APP}/Contents/MacOS/${APP_NAME}")

    if(NOT EXISTS "${EXECUTABLE}")
        message(WARNING "Executable not found: ${EXECUTABLE}")
        continue()
    endif()

    message(STATUS "Fixing rpaths: ${APP_NAME}")

    # Get current rpaths
    execute_process(
        COMMAND otool -l "${EXECUTABLE}"
        OUTPUT_VARIABLE OTOOL_OUTPUT
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # Remove all existing LC_RPATH entries (they contain build machine paths)
    string(REGEX MATCHALL "path [^\n]+" RPATH_MATCHES "${OTOOL_OUTPUT}")
    foreach(RPATH_MATCH ${RPATH_MATCHES})
        string(REGEX REPLACE "^path ([^ ]+).*" "\\1" RPATH_VALUE "${RPATH_MATCH}")
        # Skip if it's already a relative path we want
        if(RPATH_VALUE MATCHES "^@")
            continue()
        endif()
        execute_process(
            COMMAND install_name_tool -delete_rpath "${RPATH_VALUE}" "${EXECUTABLE}"
            ERROR_QUIET
        )
    endforeach()

    # Add the correct relative rpath to the parent bundle's Frameworks
    # Child exe: DEEPN++.app/Contents/Resources/Foo.app/Contents/MacOS/Foo
    # Frameworks: DEEPN++.app/Contents/Frameworks/
    # Relative:   @loader_path/../../../../Frameworks
    execute_process(
        COMMAND install_name_tool
            -add_rpath "@loader_path/../../../../Frameworks"
            "${EXECUTABLE}"
        ERROR_QUIET
    )

    # Also fix any child-specific dylibs/frameworks if present
    file(GLOB_RECURSE CHILD_DYLIBS "${CHILD_APP}/Contents/*.dylib")
    foreach(DYLIB ${CHILD_DYLIBS})
        execute_process(
            COMMAND otool -l "${DYLIB}"
            OUTPUT_VARIABLE DYLIB_OTOOL
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        string(REGEX MATCHALL "path [^\n]+" DYLIB_RPATHS "${DYLIB_OTOOL}")
        foreach(RPATH_MATCH ${DYLIB_RPATHS})
            string(REGEX REPLACE "^path ([^ ]+).*" "\\1" RPATH_VALUE "${RPATH_MATCH}")
            if(RPATH_VALUE MATCHES "^@")
                continue()
            endif()
            execute_process(
                COMMAND install_name_tool -delete_rpath "${RPATH_VALUE}" "${DYLIB}"
                ERROR_QUIET
            )
        endforeach()
    endforeach()
endforeach()

message(STATUS "Child app rpaths fixed")
