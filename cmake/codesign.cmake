# cmake/codesign.cmake
# Deep-signs a macOS .app bundle for notarization.
# Signs components in the correct order: dylibs → frameworks → nested .apps → main bundle.
#
# Usage (called via add_custom_command POST_BUILD):
#   cmake -DAPP_PATH=<path/to/Foo.app>
#         -DIDENTITY=<signing-identity>
#         -DENTITLEMENTS=<path/to/entitlements.plist>
#         -P codesign.cmake

cmake_minimum_required(VERSION 3.21)

if(NOT APP_PATH)
    message(FATAL_ERROR "codesign.cmake: APP_PATH is required")
endif()
if(NOT IDENTITY)
    message(FATAL_ERROR "codesign.cmake: IDENTITY is required")
endif()
if(NOT ENTITLEMENTS OR NOT EXISTS "${ENTITLEMENTS}")
    message(WARNING "codesign.cmake: ENTITLEMENTS not found at '${ENTITLEMENTS}', signing without custom entitlements")
    set(ENTITLEMENTS "")
endif()

# Build the base codesign flags list
set(_base_flags --force --timestamp --options runtime --verbose=0 --sign "${IDENTITY}")
if(ENTITLEMENTS)
    list(APPEND _base_flags --entitlements "${ENTITLEMENTS}")
endif()

function(_sign_item PATH)
    execute_process(
        COMMAND codesign ${_base_flags} "${PATH}"
        RESULT_VARIABLE _res
        ERROR_VARIABLE  _err
        OUTPUT_QUIET
    )
    if(NOT _res EQUAL 0)
        message(WARNING "codesign failed for ${PATH}:\n${_err}")
    else()
        message(STATUS "  signed: ${PATH}")
    endif()
endfunction()

message(STATUS "codesign: signing '${APP_PATH}' with identity '${IDENTITY}'")

# ── Step 1: Sign all loose dylibs ────────────────────────────────────────────
file(GLOB_RECURSE _dylibs "${APP_PATH}/Contents/*.dylib")
foreach(_f IN LISTS _dylibs)
    _sign_item("${_f}")
endforeach()

# ── Step 2: Sign .framework bundles (whole bundle, deepest-first) ─────────────
# Collect all .framework paths; sort so that nested frameworks come before outer ones
file(GLOB_RECURSE _fw_candidates
    LIST_DIRECTORIES true
    "${APP_PATH}/Contents/*.framework"
)
set(_frameworks "")
foreach(_d IN LISTS _fw_candidates)
    if(IS_DIRECTORY "${_d}" AND _d MATCHES "\\.framework$")
        list(APPEND _frameworks "${_d}")
    endif()
endforeach()
# Reverse-sort by path depth so deepest (most nested) frameworks sign first
list(SORT _frameworks ORDER DESCENDING)
foreach(_fw IN LISTS _frameworks)
    _sign_item("${_fw}")
endforeach()

# ── Step 3: Sign nested .app bundles in Resources/ ───────────────────────────
file(GLOB _nested_apps "${APP_PATH}/Contents/Resources/*.app")
foreach(_nested IN LISTS _nested_apps)
    if(IS_DIRECTORY "${_nested}")
        _sign_item("${_nested}")
    endif()
endforeach()

# ── Step 4: Sign the main .app bundle itself ──────────────────────────────────
_sign_item("${APP_PATH}")

message(STATUS "codesign: done")
