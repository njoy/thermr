
add_executable( inelastic.e_mu_ep_util.test
                mainLoop.test.cpp
                sigu.test.cpp )
target_compile_options( inelastic.e_mu_ep_util.test PRIVATE ${${PREFIX}_common_flags}
$<$<BOOL:${strict}>:${${PREFIX}_strict_flags}>$<$<CONFIG:DEBUG>:
${${PREFIX}_DEBUG_flags}
$<$<BOOL:${coverage}>:${${PREFIX}_coverage_flags}>>
$<$<CONFIG:RELEASE>:
${${PREFIX}_RELEASE_flags}
$<$<BOOL:${link_time_optimization}>:${${PREFIX}_link_time_optimization_flags}>
$<$<BOOL:${nonportable_optimization}>:${${PREFIX}_nonportable_optimization_flags}>>

${CXX_appended_flags} ${thermr_appended_flags} )
target_link_libraries( inelastic.e_mu_ep_util.test PUBLIC thermr ) 
add_test( NAME inelastic.e_mu_ep_util COMMAND inelastic.e_mu_ep_util.test )
