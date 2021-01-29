message(STATUS "Generating test binaries")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)
add_link_options(-Wl,-rpath=${CMAKE_INSTALL_PREFIX}/lib)

IF(DEFINED ENV{INTELCOMP_HOME})
#CMAKE_CXX_COMPILER_ID MATCHES "Intel" )
   add_link_options(-Wl,-rpath $ENV{INTELCOMP_HOME}/lib/intel64_lin) 
   add_link_options(-lintlc -lifcore -lifport -fopenmp)
endif()

#add_executable(interpolate interpolate.c)
#target_link_libraries(interpolate georef -lm ${RMN_LIBRARIES} ${VGRID_LIBRARIES})
#set_target_properties(interpolate PROPERTIES VERSION ${PROJECT_VERSION})

#install(TARGETS interpolate
#   RUNTIME DESTINATION bin)
