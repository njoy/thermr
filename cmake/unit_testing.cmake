#######################################################################
# Setup
#######################################################################

message( STATUS "Adding thermr unit testing" )
enable_testing()


#######################################################################
# Unit testing directories
#######################################################################

add_subdirectory( src/approximations/freeGas/test )
add_subdirectory( src/inelastic/e_ep_mu_util/test )
add_subdirectory( src/inelastic/e_mu_ep_util/test )
add_subdirectory( src/inelastic/sig_util/test )
add_subdirectory( src/inelastic/test )
add_subdirectory( src/coherentElastic/coherentElastic_util/sigcoh_util/test )
add_subdirectory( src/coherentElastic/coherentElastic_util/test )
add_subdirectory( src/coherentElastic/test )
add_subdirectory( src/general_util/test )
add_subdirectory( src/incoherentElastic/test )
add_subdirectory( src/test )
