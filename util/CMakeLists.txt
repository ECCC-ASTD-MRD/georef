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

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../lib)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../lib)

target_link_libraries(EZTiler georef -lm ${RMN_LIBRARIES} ${VGRID_LIBRARIES})
set_target_properties(EZTiler PROPERTIES VERSION ${PROJECT_VERSION})

target_link_libraries(ReGrid georef -lm ${RMN_LIBRARIES} ${VGRID_LIBRARIES})
set_target_properties(ReGrid PROPERTIES VERSION ${PROJECT_VERSION})

target_link_libraries(Interpolate georef -lm ${RMN_LIBRARIES} ${VGRID_LIBRARIES})
set_target_properties(Interpolate PROPERTIES VERSION ${PROJECT_VERSION})

install(TARGETS EZTiler ReGrid Interpolate
   RUNTIME DESTINATION bin)

set(LIBGEOREF_DATA /local/drive1/afsr005/GRIDS)
#add_test(GridInterp_L2L   Interpolate -i ${LIBGEOREF_DATA}/L.fstd -g ${LIBGEOREF_DATA}/L.fstd -o /tmp/libgeoref.fstd -n GRID -e L2L -v DEBUG)
#add_test(GridInterp_PS2L  Interpolate -i ${LIBGEOREF_DATA}/PS.fstd -g ${LIBGEOREF_DATA}/L.fstd -o /tmp/libgeoref.fstd -n GRID -e PS2L -v DEBUG)
#add_test(GridInterp_PSN2L Interpolate -i ${LIBGEOREF_DATA}/PSN.fstd -g ${LIBGEOREF_DATA}/L.fstd -o /tmp/libgeoref.fstd -n GRID -e PSN2L -v DEBUG)
#add_test(GridInterp_PSS2L Interpolate -i ${LIBGEOREF_DATA}/PSS.fstd -g ${LIBGEOREF_DATA}/L.fstd -o /tmp/libgeoref.fstd -n GRID -e PSS2L -v DEBUG)
#add_test(GridInterp_ZE2L  Interpolate -i ${LIBGEOREF_DATA}/ZE.fstd -g ${LIBGEOREF_DATA}/L.fstd -o /tmp/libgeoref.fstd -n GRID -e ZE2L -v DEBUG)
#add_test(GridInterp_ZL2L  Interpolate -i ${LIBGEOREF_DATA}/ZL.fstd -g ${LIBGEOREF_DATA}/L.fstd -o /tmp/libgeoref.fstd -n GRID -e ZL2L -v DEBUG)
#add_test(GridInterp_ZW2L  Interpolate -i ${LIBGEOREF_DATA}/ZW.fstd -g ${LIBGEOREF_DATA}/L.fstd -o /tmp/libgeoref.fstd -n GRID -e ZW2L -v DEBUG)
#add_test(GridInterp_O2L   Interpolate -i ${LIBGEOREF_DATA}/O.fstd -g ${LIBGEOREF_DATA}/L.fstd -o /tmp/libgeoref.fstd -n GRID -e O2L -v DEBUG)
#add_test(GridInterp_M2L   Interpolate -i ${LIBGEOREF_DATA}/M.fstd -g ${LIBGEOREF_DATA}/L.fstd -o /tmp/libgeoref.fstd -n GRID -e M2L -v DEBUG)

#add_test(GridInterp_L2PS  Interpolate -g ${LIBGEOREF_DATA}/PS.fstd -i ${LIBGEOREF_DATA}/L.fstd -o /tmp/libgeoref.fstd -n GRID -e L2PS -v DEBUG)
#add_test(GridInterp_L2PSN Interpolate -g ${LIBGEOREF_DATA}/PSN.fstd -i ${LIBGEOREF_DATA}/L.fstd -o /tmp/libgeoref.fstd -n GRID -e L2PSN -v DEBUG)
#add_test(GridInterp_L2PSS Interpolate -g ${LIBGEOREF_DATA}/PSS.fstd -i ${LIBGEOREF_DATA}/L.fstd -o /tmp/libgeoref.fstd -n GRID -e L2PSS -v DEBUG)
#add_test(GridInterp_L2ZE  Interpolate -g ${LIBGEOREF_DATA}/ZE.fstd -i ${LIBGEOREF_DATA}/L.fstd -o /tmp/libgeoref.fstd -n GRID -e L2ZE -v DEBUG)
#add_test(GridInterp_L2ZL  Interpolate -g ${LIBGEOREF_DATA}/ZL.fstd -i ${LIBGEOREF_DATA}/L.fstd -o /tmp/libgeoref.fstd -n GRID -e L2ZL -v DEBUG)
#add_test(GridInterp_L2ZW  Interpolate -g ${LIBGEOREF_DATA}/ZW.fstd -i ${LIBGEOREF_DATA}/L.fstd -o /tmp/libgeoref.fstd -n GRID -e L2ZW -v DEBUG)
#add_test(GridInterp_L2O   Interpolate -g ${LIBGEOREF_DATA}/O.fstd -i ${LIBGEOREF_DATA}/L.fstd -o /tmp/libgeoref.fstd -n GRID -e L2O -v DEBUG)
#add_test(GridInterp_L2M   Interpolate -g ${LIBGEOREF_DATA}/M.fstd -i ${LIBGEOREF_DATA}/L.fstd -o /tmp/libgeoref.fstd -n GRID -e L2M -v DEBUG)

add_test(GridInterp_ZE2PS  Interpolate -g ${LIBGEOREF_DATA}/PS.fstd -i ${LIBGEOREF_DATA}/HRDPS.fstd -o /tmp/libgeoref.fstd -n TT -e ZE2PS -v DEBUG)
add_test(GridInterp_ZE2PSN Interpolate -g ${LIBGEOREF_DATA}/PSN.fstd -i ${LIBGEOREF_DATA}/HRDPS.fstd -o /tmp/libgeoref.fstd -n TT -e ZE2PSN -v DEBUG)
add_test(GridInterp_ZE2PSS Interpolate -g ${LIBGEOREF_DATA}/PSS.fstd -i ${LIBGEOREF_DATA}/HRDPS.fstd -o /tmp/libgeoref.fstd -n TT -e ZE2PSS -v DEBUG)
add_test(GridInterp_ZE2ZE  Interpolate -g ${LIBGEOREF_DATA}/ZE.fstd -i ${LIBGEOREF_DATA}/HRDPS.fstd -o /tmp/libgeoref.fstd -n TT -e ZE2ZE -v DEBUG)
add_test(GridInterp_ZE2ZL  Interpolate -g ${LIBGEOREF_DATA}/ZL.fstd -i ${LIBGEOREF_DATA}/HRDPS.fstd -o /tmp/libgeoref.fstd -n TT -e ZE2ZL -v DEBUG)
add_test(GridInterp_ZE2ZW  Interpolate -g ${LIBGEOREF_DATA}/ZW.fstd -i ${LIBGEOREF_DATA}/HRDPS.fstd -o /tmp/libgeoref.fstd -n TT -e ZE2ZW -v DEBUG)
add_test(GridInterp_ZE2O   Interpolate -g ${LIBGEOREF_DATA}/O.fstd -i ${LIBGEOREF_DATA}/HRDPS.fstd -o /tmp/libgeoref.fstd -n TT -e ZE2O -v DEBUG)
add_test(GridInterp_ZE2M   Interpolate -g ${LIBGEOREF_DATA}/M.fstd -i ${LIBGEOREF_DATA}/HRDPS.fstd -o /tmp/libgeoref.fstd -n TT -e ZE2M -v DEBUG)
