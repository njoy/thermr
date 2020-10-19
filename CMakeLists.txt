########################################################################
# Preamble
########################################################################

cmake_minimum_required( VERSION 3.14 )
project( thermr LANGUAGES CXX )


########################################################################
# Project-wide setup
########################################################################

set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED YES )

option( thermr_unit_tests
    "Compile the thermr unit tests and integrate with ctest" ON
    )
option( strict_compile
    "Treat all warnings as errors." ON
    )

# Compile flags
set( common_flags "-Wall" "-Wextra" "-Wpedantic" )
set( strict_flags "-Werror" )
set( release_flags "-O3" )
set( debug_flags "-O0" "-g" )


########################################################################
# Dependencies
########################################################################

set( REPOSITORIES "release"
    CACHE STRING
    "Options for where to fetch repositories: develop, release, local"
    )

if( REPOSITORIES STREQUAL "develop" )
    include( cmake/develop_dependencies.cmake )

elseif( REPOSITORIES STREQUAL "release" )
    include( cmake/release_dependencies.cmake )

elseif( REPOSITORIES STREQUAL "local" )
    include( cmake/local_dependencies.cmake )

endif()


########################################################################
# Project targets
########################################################################

add_library( thermr 
    src/thermr.cpp
    )
target_include_directories( thermr PUBLIC src/ )
target_link_libraries( thermr
    PUBLIC ENDFtk
    PUBLIC catch-adapter
    PUBLIC range-v3-adapter
    )
target_compile_options( thermr PRIVATE
    ${common_flags}
    $<$<BOOL:${strict_compile}>:${strict_flags}>
    $<$<CONFIG:DEBUG>:${debug_flags}>
    $<$<CONFIG:RELEASE>:${release_flags}>
    )


#######################################################################
# Top-level Only
#######################################################################

if( CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR )

    # unit testing
    if( thermr_unit_tests )
        include( cmake/unit_testing.cmake )
    endif()

endif()
