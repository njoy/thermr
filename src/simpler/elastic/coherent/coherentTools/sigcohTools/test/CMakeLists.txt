
add_executable( simpler.elastic.coherent.coherentTools.sigcohTools.test
                form.test.cpp
                legndr.test.cpp
                terp.test.cpp )
target_compile_options( simpler.elastic.coherent.coherentTools.sigcohTools.test PRIVATE ${${PREFIX}_common_flags}
$<$<BOOL:${strict}>:${${PREFIX}_strict_flags}>$<$<CONFIG:DEBUG>:
${${PREFIX}_DEBUG_flags}
$<$<BOOL:${coverage}>:${${PREFIX}_coverage_flags}>>
$<$<CONFIG:RELEASE>:
${${PREFIX}_RELEASE_flags}
$<$<BOOL:${link_time_optimization}>:${${PREFIX}_link_time_optimization_flags}>
$<$<BOOL:${nonportable_optimization}>:${${PREFIX}_nonportable_optimization_flags}>>

${CXX_appended_flags} ${thermr_appended_flags} )
target_link_libraries( simpler.elastic.coherent.coherentTools.sigcohTools.test PUBLIC thermr ) 
add_test( NAME simpler.elastic.coherent.coherentTools.sigcohTools COMMAND simpler.elastic.coherent.coherentTools.sigcohTools.test )