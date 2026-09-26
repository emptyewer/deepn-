# fix_frameworks.cmake
# Fixes Qt framework bundle structure for codesigning.
# macdeployqt in Qt 6 sometimes creates Versions/Current as a directory
# instead of a symlink to Versions/A, causing "bundle format is ambiguous".
#
# Usage: cmake -DAPP_PATH=/path/to/Foo.app -P fix_frameworks.cmake

if(NOT APP_PATH)
    message(FATAL_ERROR "APP_PATH is required")
endif()

file(GLOB_RECURSE _fw_candidates
    LIST_DIRECTORIES true
    "${APP_PATH}/*.framework"
)

set(_frameworks "")
foreach(_d IN LISTS _fw_candidates)
    if(IS_DIRECTORY "${_d}" AND _d MATCHES "\\.framework$")
        list(APPEND _frameworks "${_d}")
    endif()
endforeach()

foreach(_fw IN LISTS _frameworks)
    set(_versions_dir "${_fw}/Versions")
    set(_current "${_versions_dir}/Current")
    set(_a_dir "${_versions_dir}/A")

    # Skip if no Versions/A exists
    if(NOT IS_DIRECTORY "${_a_dir}")
        continue()
    endif()

    # Fix Versions/Current: should be symlink to A
    if(IS_DIRECTORY "${_current}" AND NOT IS_SYMLINK "${_current}")
        file(REMOVE_RECURSE "${_current}")
        execute_process(COMMAND ln -sf A "${_current}")
    endif()

    # Fix top-level entries: should be symlinks into Versions/Current
    get_filename_component(_fw_name "${_fw}" NAME_WE)

    # Fix the main binary symlink
    set(_top_binary "${_fw}/${_fw_name}")
    if(EXISTS "${_a_dir}/${_fw_name}" AND NOT IS_SYMLINK "${_top_binary}")
        file(REMOVE "${_top_binary}")
        execute_process(COMMAND ln -sf "Versions/Current/${_fw_name}" "${_top_binary}")
    endif()

    # Fix Resources symlink
    set(_top_resources "${_fw}/Resources")
    if(IS_DIRECTORY "${_a_dir}/Resources" AND NOT IS_SYMLINK "${_top_resources}")
        file(REMOVE_RECURSE "${_top_resources}")
        execute_process(COMMAND ln -sf "Versions/Current/Resources" "${_top_resources}")
    endif()

    # Fix Headers symlink (if present)
    set(_top_headers "${_fw}/Headers")
    if(IS_DIRECTORY "${_a_dir}/Headers" AND NOT IS_SYMLINK "${_top_headers}")
        file(REMOVE_RECURSE "${_top_headers}")
        execute_process(COMMAND ln -sf "Versions/Current/Headers" "${_top_headers}")
    endif()
endforeach()

message(STATUS "Framework symlinks fixed in ${APP_PATH}")
