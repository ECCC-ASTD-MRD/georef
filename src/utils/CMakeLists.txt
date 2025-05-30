add_link_options(-lpthread)

IF(DEFINED ENV{INTELCOMP_HOME})
    add_link_options(-Wl,-rpath $ENV{INTELCOMP_HOME}/lib/intel64_lin) 
    add_link_options(-lintlc -lifcore -lifport)
endif()

add_executable(georef_interp georef_interp.c)

target_link_libraries(georef_interp georef-ompi-static rmn::rmn-ompi m)
set_target_properties(georef_interp PROPERTIES VERSION ${PROJECT_VERSION})

install(TARGETS georef_interp RUNTIME DESTINATION bin)
