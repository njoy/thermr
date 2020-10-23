
add_executable( coherentElastic.coherentElastic_util.sigcoh_util.test
                form.test.cpp
                legndr.test.cpp
                terp.test.cpp )
target_compile_options( coherentElastic.coherentElastic_util.sigcoh_util.test PRIVATE ${${PREFIX}_common_flags}
$<$<BOOL:${strict}>:${${PREFIX}_strict_flags}>$<$<CONFIG:DEBUG>:
${${PREFIX}_DEBUG_flags}
$<$<BOOL:${coverage}>:${${PREFIX}_coverage_flags}>>
$<$<CONFIG:RELEASE>:
${${PREFIX}_RELEASE_flags}
$<$<BOOL:${link_time_optimization}>:${${PREFIX}_link_time_optimization_flags}>
$<$<BOOL:${nonportable_optimization}>:${${PREFIX}_nonportable_optimization_flags}>>

${CXX_appended_flags} ${thermr_appended_flags} )
target_link_libraries( coherentElastic.coherentElastic_util.sigcoh_util.test PUBLIC thermr ) 
add_test( NAME coherentElastic.coherentElastic_util.sigcoh_util COMMAND coherentElastic.coherentElastic_util.sigcoh_util.test )
