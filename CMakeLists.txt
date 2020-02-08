
cmake_minimum_required( VERSION 3.2 ) 
set( CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Supported configuration types" FORCE )
        
if( NOT ROOT_DIRECTORY )
    set( ROOT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    if ( NOT fetched_subprojects )
        if ( NOT PYTHON_EXECUTABLE )
            find_package( PythonInterp )
            if ( NOT PYTHONINTERP_FOUND )
                message( FATAL_ERROR "Python interpeter installation was not found." )
            endif()
        endif()
        execute_process( COMMAND ${PYTHON_EXECUTABLE} "./metaconfigure/fetch_subprojects.py"
                         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
                         RESULT_VARIABLE fetch_failure )
        if ( NOT fetch_failure )
            set( fetched_subprojects TRUE CACHE BOOL "fetch script ran" )
        else()
            message( FATAL_ERROR "Failed to fetch dependencies" )
        endif()
    endif()
endif()

project( thermr LANGUAGES CXX )
get_directory_property( is_subproject PARENT_DIRECTORY )
include( CMakeDependentOption REQUIRED )

set( thermr_AppleClang_minimum_version 8.0 )
set( thermr_GNU_minimum_version 6.2 )
set( thermr_Clang_minimum_version 3.8 )

if( thermr_${CMAKE_CXX_COMPILER_ID}_minimum_version )
    if( CMAKE_CXX_COMPILER_VERSION AND
        CMAKE_CXX_COMPILER_VERSION VERSION_LESS
        ${thermr_${CMAKE_CXX_COMPILER_ID}_minimum_version} )
        message( FATAL_ERROR "${CMAKE_CXX_COMPILER_ID} version must be greater than ${thermr_${CMAKE_CXX_COMPILER_ID}_minimum_version}" )
    endif()
endif()

# general properties
option( thermr_strict "Compile time warnings are converted to errors" ON )

# binary instrumentation
option( coverage "Enable binary instrumentation to collect test coverage information in the DEBUG configuration" )
option( profile_generate "Enable binary instrumentation to generation execution profiles in the RELEASE configuration which may be used to guide later optimization" )

# additional optimizations
option( link_time_optimization "Enable link time optimization in the RELEASE configuration" )
option( profile_use "In the RELEASE configuration, leverage previously generated exeution profile to inform optimization decisions" )
option( nonportable_optimization "Enable optimizations which compromise portability of resulting binary in the RELEASE configuration" )

# libraries and linking
option( static "Statically link component and environment libraries" OFF )
if ( static AND ( "${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin" ) )
    message( FATAL_ERROR "Static binaries not supported on OSX" )
endif()

CMAKE_DEPENDENT_OPTION( static_libraries "Statically link component libraries" OFF "NOT static" ON )
CMAKE_DEPENDENT_OPTION( static_thermr "Statically link the thermr component library" OFF "NOT static;NOT static_libraries" ON ) 

option( unit_tests "Compile the thermr unit tests and integrate with ctest" ON ) 

if ( profile_generate AND profile_use )
    message( FATAL_ERROR "Cannot both generate and use execution profile in the same configuration" )
endif()

if ( profile_generate )
    file( MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/profiling" )
endif()

set( thermr_AppleClang_Darwin_common_flags "-stdlib=libc++" "-Wall" "-Wextra" "-Wpedantic" "-std=c++14" )
set( thermr_AppleClang_Darwin_DEBUG_flags "-O0" "-g" "-gdwarf-3" )
set( thermr_AppleClang_Darwin_RELEASE_flags "-O3" "-DNDEBUG" )
set( thermr_AppleClang_Darwin_strict_flags "-Werror" )
set( thermr_AppleClang_Darwin_coverage_flags "--coverage" )
set( thermr_AppleClang_Darwin_subproject_flags  )
set( thermr_AppleClang_Darwin_base_project_flags  )
set( thermr_AppleClang_Darwin_profile_generate_flags "-fprofile-generate='${CMAKE_BINARY_DIR}/profiling'" )
set( thermr_AppleClang_Darwin_link_time_optimization_flags "-flto" )
set( thermr_AppleClang_Darwin_profile_use_flags "-fprofile-use='${CMAKE_BINARY_DIR}/profiling'" )
set( thermr_AppleClang_Darwin_nonportable_optimization_flags "-march=native" )
set( thermr_AppleClang_Darwin_static_flags "-static" )
set( thermr_GNU_Linux_common_flags "-Wall" "-Wextra" "-Wpedantic" "-std=c++14" )
set( thermr_GNU_Linux_DEBUG_flags "-O0" "-g" "-gdwarf-3" "-fsignaling-nans" )
set( thermr_GNU_Linux_RELEASE_flags "-O3" "-DNDEBUG" )
set( thermr_GNU_Linux_strict_flags "-Werror" )
set( thermr_GNU_Linux_coverage_flags "--coverage" )
set( thermr_GNU_Linux_subproject_flags  )
set( thermr_GNU_Linux_base_project_flags  )
set( thermr_GNU_Linux_profile_generate_flags "-fprofile-generate='${CMAKE_BINARY_DIR}/profiling'" )
set( thermr_GNU_Linux_link_time_optimization_flags "-flto" )
set( thermr_GNU_Linux_profile_use_flags "-fprofile-use='${CMAKE_BINARY_DIR}/profiling'" )
set( thermr_GNU_Linux_nonportable_optimization_flags "-march=native" )
set( thermr_GNU_Linux_static_flags "-static" )
set( thermr_GNU_Windows_common_flags "-Wall" "-Wextra" "-Wpedantic" "-std=c++14" )
set( thermr_GNU_Windows_DEBUG_flags "-O0" "-g" "-gdwarf-3" "-fsignaling-nans" )
set( thermr_GNU_Windows_RELEASE_flags "-O3" "-DNDEBUG" )
set( thermr_GNU_Windows_strict_flags "-Werror" )
set( thermr_GNU_Windows_coverage_flags "--coverage" )
set( thermr_GNU_Windows_subproject_flags  )
set( thermr_GNU_Windows_base_project_flags  )
set( thermr_GNU_Windows_profile_generate_flags "-fprofile-generate='${CMAKE_BINARY_DIR}/profiling'" )
set( thermr_GNU_Windows_link_time_optimization_flags "-flto" )
set( thermr_GNU_Windows_profile_use_flags "-fprofile-use='${CMAKE_BINARY_DIR}/profiling'" )
set( thermr_GNU_Windows_nonportable_optimization_flags "-march=native" )
set( thermr_GNU_Windows_static_flags "-static" )
set( thermr_GNU_Darwin_common_flags "-Wall" "-Wextra" "-Wpedantic" "-std=c++14" )
set( thermr_GNU_Darwin_DEBUG_flags "-O0" "-g" "-gdwarf-3" "-fsignaling-nans" )
set( thermr_GNU_Darwin_RELEASE_flags "-O3" "-DNDEBUG" )
set( thermr_GNU_Darwin_strict_flags "-Werror" )
set( thermr_GNU_Darwin_coverage_flags "--coverage" )
set( thermr_GNU_Darwin_subproject_flags  )
set( thermr_GNU_Darwin_base_project_flags  )
set( thermr_GNU_Darwin_profile_generate_flags "-fprofile-generate='${CMAKE_BINARY_DIR}/profiling'" )
set( thermr_GNU_Darwin_link_time_optimization_flags "-flto" )
set( thermr_GNU_Darwin_profile_use_flags "-fprofile-use='${CMAKE_BINARY_DIR}/profiling'" )
set( thermr_GNU_Darwin_nonportable_optimization_flags "-march=native" )
set( thermr_GNU_Darwin_static_flags "-static" )
set( thermr_Clang_Linux_common_flags "-stdlib=libstdc++" "-Wall" "-Wextra" "-Wpedantic" "-std=c++14" )
set( thermr_Clang_Linux_DEBUG_flags "-O0" "-g" "-gdwarf-3" )
set( thermr_Clang_Linux_RELEASE_flags "-O3" "-DNDEBUG" )
set( thermr_Clang_Linux_strict_flags "-Werror" )
set( thermr_Clang_Linux_coverage_flags "--coverage" )
set( thermr_Clang_Linux_subproject_flags  )
set( thermr_Clang_Linux_base_project_flags  )
set( thermr_Clang_Linux_profile_generate_flags "-fprofile-generate='${CMAKE_BINARY_DIR}/profiling'" )
set( thermr_Clang_Linux_link_time_optimization_flags "-flto" )
set( thermr_Clang_Linux_profile_use_flags "-fprofile-use='${CMAKE_BINARY_DIR}/profiling'" )
set( thermr_Clang_Linux_nonportable_optimization_flags "-march=native" )
set( thermr_Clang_Linux_static_flags "-static" )
set( thermr_Clang_Windows_common_flags "-stdlib=libc++" "-Wall" "-Wextra" "-Wpedantic" "-std=c++14" )
set( thermr_Clang_Windows_DEBUG_flags "-O0" "-g" "-gdwarf-3" )
set( thermr_Clang_Windows_RELEASE_flags "-O3" "-DNDEBUG" )
set( thermr_Clang_Windows_strict_flags "-Werror" )
set( thermr_Clang_Windows_coverage_flags "--coverage" )
set( thermr_Clang_Windows_subproject_flags  )
set( thermr_Clang_Windows_base_project_flags  )
set( thermr_Clang_Windows_profile_generate_flags "-fprofile-generate='${CMAKE_BINARY_DIR}/profiling'" )
set( thermr_Clang_Windows_link_time_optimization_flags "-flto" )
set( thermr_Clang_Windows_profile_use_flags "-fprofile-use='${CMAKE_BINARY_DIR}/profiling'" )
set( thermr_Clang_Windows_nonportable_optimization_flags "-march=native" )
set( thermr_Clang_Windows_static_flags "-static" )
set( thermr_Clang_Darwin_common_flags "-stdlib=libc++" "-Wall" "-Wextra" "-Wpedantic" "-std=c++14" )
set( thermr_Clang_Darwin_DEBUG_flags "-O0" "-g" "-gdwarf-3" )
set( thermr_Clang_Darwin_RELEASE_flags "-O3" "-DNDEBUG" )
set( thermr_Clang_Darwin_strict_flags "-Werror" )
set( thermr_Clang_Darwin_coverage_flags "--coverage" )
set( thermr_Clang_Darwin_subproject_flags  )
set( thermr_Clang_Darwin_base_project_flags  )
set( thermr_Clang_Darwin_profile_generate_flags "-fprofile-generate='${CMAKE_BINARY_DIR}/profiling'" )
set( thermr_Clang_Darwin_link_time_optimization_flags "-flto" )
set( thermr_Clang_Darwin_profile_use_flags "-fprofile-use='${CMAKE_BINARY_DIR}/profiling'" )
set( thermr_Clang_Darwin_nonportable_optimization_flags "-march=native" )
set( thermr_Clang_Darwin_static_flags "-static" )

if ( static_thermr )
    set( thermr_library_linkage STATIC )
else ()
    set( thermr_library_linkage SHARED )
endif () 

set( CMAKE_SKIP_BUILD_RPATH FALSE )
set( CMAKE_BUILD_WITH_INSTALL_RPATH FALSE )
if ( CMAKE_SYSTEM_NAME STREQUAL "Darwin" )
    set( rpath_prefix "@loader_path" )
else()
    set( rpath_prefix "\\$ORIGIN" )
endif()
list( INSERT 0 CMAKE_INSTALL_RPATH "${rpath_prefix}/../lib" )
set( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )
get_directory_property( is_subproject PARENT_DIRECTORY )

if( NOT TARGET catch-adapter )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/catch-adapter )
endif()

if( NOT TARGET disco )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/disco )
endif()

if( NOT TARGET hana-adapter )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/hana-adapter )
endif()

if( NOT TARGET fmt-adapter )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/fmt-adapter )
endif()

if( NOT TARGET spdlog-adapter )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/spdlog-adapter )
endif()

if( NOT TARGET Log )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/Log )
endif()

if( NOT TARGET header-utilities )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/header-utilities )
endif()

if( NOT TARGET range-v3-adapter )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/range-v3-adapter )
endif()

if( NOT TARGET variant-adapter )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/variant-adapter )
endif()

if( NOT TARGET ENDFtk )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/ENDFtk )
endif()


if ( NOT GIT_EXECUTABLE )
    find_package( Git )
    if ( NOT GIT_FOUND )
        message( FATAL_ERROR "git installation was not found." )
    endif()
endif()

execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
) 

message( STATUS "" )
message( STATUS "-----------------------------------------------------------" )
message( STATUS "" )
message( STATUS "thermr" )
message( STATUS "Git current branch: ${GIT_BRANCH}" )
message( STATUS "Git commit hash: ${GIT_HASH}" )
message( STATUS "" )
message( STATUS "-----------------------------------------------------------" )

add_library( thermr ${thermr_library_linkage} 
             "${CMAKE_CURRENT_SOURCE_DIR}/src/approximations/freeGas/freeGas.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/calcem/calcem.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/calcem/calcem_util/e_ep_mu.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/calcem/calcem_util/e_mu_ep.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/calcem/calcem_util/sig.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/calcem/calcem_util/e_ep_mu_util/ep_mu.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/calcem/calcem_util/e_ep_mu_util/sigl.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/calcem/calcem_util/e_mu_ep_util/adaptiveReconstruction.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/calcem/calcem_util/e_mu_ep_util/mainLoop.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/calcem/calcem_util/e_mu_ep_util/sigu.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/calcem/calcem_util/sig_util/terp1.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/calcem/calcem_util/sig_util/terpq.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/coh/coh.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/coh/coh_util/sigcoh.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/coh/coh_util/upstk.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/coh/coh_util/sigcoh_util/form.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/coh/coh_util/sigcoh_util/legndr.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/coh/coh_util/sigcoh_util/terp.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/general_util/sigfig.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/generalTools/testing.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/iel/incoherentElastic.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/temp_eff/temp_eff.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/temp_eff/temp_eff_util/gatef2.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/temp_eff/temp_eff_util/gateff.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/thermr.cpp" )

target_include_directories( thermr PUBLIC src/ )

set( PREFIX thermr_${CMAKE_CXX_COMPILER_ID}_${CMAKE_SYSTEM_NAME} )

target_compile_options( thermr PRIVATE 
${${PREFIX}_common_flags}
$<$<BOOL:${thermr_strict}>:${${PREFIX}_strict_flags}>
$<$<BOOL:${static}>:${${PREFIX}_static_flags}>
$<$<BOOL:${is_subproject}>:${${PREFIX}_subproject_flags}>
$<$<NOT:$<BOOL:${is_subproject}>>:${${PREFIX}_base_project_flags}>
$<$<CONFIG:DEBUG>:
${${PREFIX}_DEBUG_flags} 
$<$<BOOL:${coverage}>:${${PREFIX}_coverage_flags}>>
$<$<CONFIG:RELEASE>:
${${PREFIX}_RELEASE_flags} 
$<$<BOOL:${profile_generate}>:${${PREFIX}_profile_generate_flags}> 
$<$<BOOL:${profile_use}>:${${PREFIX}_profile_use_flags}> 
$<$<BOOL:${link_time_optimization}>:${${PREFIX}_link_time_optimization_flags}> 
$<$<BOOL:${nonportable_optimization}>:${${PREFIX}_nonportable_optimization_flags}>>
${CXX_appended_flags} ${thermr_appended_flags} )

target_link_libraries( thermr PUBLIC "$<$<AND:$<CONFIG:RELEASE>,$<BOOL:${link_time_optimization}>>:${${PREFIX}_RELEASE_flags};${${PREFIX}_link_time_optimization_flags}$<$<BOOL:${profile_generate}>:${${PREFIX}_profile_generate_flags};>$<$<BOOL:${profile_use}>:${${PREFIX}_profile_use_flags};>$<$<BOOL:${nonportable_optimization}>:${${PREFIX}_nonportable_optimization_flags};>>$<$<CONFIG:DEBUG>:$<$<BOOL:${coverage}>:${${PREFIX}_coverage_flags};>>$<$<BOOL:CXX_appended_flags>:${CXX_appended_flags};>$<$<BOOL:thermr_appended_flags>:${thermr_appended_flags};>" )

target_link_libraries( thermr PUBLIC catch-adapter PUBLIC ENDFtk PUBLIC range-v3-adapter )

if( NOT is_subproject )
    enable_testing() 
    if ( unit_tests )
        add_subdirectory( src/approximations/freeGas/test )
        add_subdirectory( src/calcem/test )
        add_subdirectory( src/calcem/calcem_util/test )
        add_subdirectory( src/calcem/calcem_util/e_ep_mu_util/test )
        add_subdirectory( src/calcem/calcem_util/e_mu_ep_util/test )
        add_subdirectory( src/calcem/calcem_util/sig_util/test )
        add_subdirectory( src/coh/test )
        add_subdirectory( src/coh/coh_util/test )
        add_subdirectory( src/coh/coh_util/sigcoh_util/test )
        add_subdirectory( src/general_util/test )
        add_subdirectory( src/iel/test )
        add_subdirectory( src/temp_eff/test )
        add_subdirectory( src/temp_eff/temp_eff_util/test )
        add_subdirectory( src/test )
    endif() 
endif()

install( DIRECTORY src// DESTINATION include
         FILE_PERMISSIONS OWNER_READ OWNER_WRITE 
                          GROUP_READ 
                          WORLD_READ
         DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE 
                               GROUP_READ 
                               WORLD_READ
         FILES_MATCHING REGEX ".*\.h[+][+]$|.*\.hxx$|.*\.hpp$|.*\.hh$|.*\.h$"  )
