message(STATUS "Generating libgeoref binaries")

#add_link_options(-Wl,-rpath=${CMAKE_INSTALL_PREFIX}/lib)
add_link_options(-lpthread)

IF(DEFINED ENV{INTELCOMP_HOME})
#CMAKE_CXX_COMPILER_ID MATCHES "Intel" )
   add_link_options(-Wl,-rpath $ENV{INTELCOMP_HOME}/lib/intel64_lin) 
   add_link_options(-lintlc -lifcore -lifport)
endif()

add_executable(EZTiler EZTiler.c)
add_executable(ReGrid ReGrid.c)
add_executable(Interpolate Interpolate.c)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)

target_link_libraries(EZTiler georef -lm ${RMN_LIBRARIES} ${VGRID_LIBRARIES})
set_target_properties(EZTiler PROPERTIES VERSION ${PROJECT_VERSION})

target_link_libraries(ReGrid georef -lm ${RMN_LIBRARIES} ${VGRID_LIBRARIES})
set_target_properties(ReGrid PROPERTIES VERSION ${PROJECT_VERSION})

target_link_libraries(Interpolate georef -lm ${RMN_LIBRARIES} ${VGRID_LIBRARIES})
set_target_properties(Interpolate PROPERTIES VERSION ${PROJECT_VERSION})

install(TARGETS EZTiler ReGrid Interpolate
   RUNTIME DESTINATION bin)