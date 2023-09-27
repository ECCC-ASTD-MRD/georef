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

target_link_libraries(EZTiler georef -lm rmn::rmn-ompi)
set_target_properties(EZTiler PROPERTIES VERSION ${PROJECT_VERSION})

target_link_libraries(ReGrid georef -lm rmn::rmn-ompi)
set_target_properties(ReGrid PROPERTIES VERSION ${PROJECT_VERSION})

target_link_libraries(Interpolate georef -lm rmn::rmn-ompi)
set_target_properties(Interpolate PROPERTIES VERSION ${PROJECT_VERSION})

install(TARGETS EZTiler ReGrid Interpolate
   RUNTIME DESTINATION bin)

#set(LIBGEOREF_DATA /local/drive1/afsr005/GRIDS)
#add_test(GridInterp_L2L   Interpolate -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/L.fstd -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/L.fstd -o /tmp/libgeoref.fstd -n GRID -e L2L -v DEBUG)
#add_test(GridInterp_PS2L  Interpolate -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/PS.fstd -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/L.fstd -o /tmp/libgeoref.fstd -n GRID -e PS2L -v DEBUG)
#add_test(GridInterp_PSN2L Interpolate -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/PSN.fstd -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/L.fstd -o /tmp/libgeoref.fstd -n GRID -e PSN2L -v DEBUG)
#add_test(GridInterp_PSS2L Interpolate -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/PSS.fstd -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/L.fstd -o /tmp/libgeoref.fstd -n GRID -e PSS2L -v DEBUG)
#add_test(GridInterp_ZE2L  Interpolate -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/ZE.fstd -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/L.fstd -o /tmp/libgeoref.fstd -n GRID -e ZE2L -v DEBUG)
#add_test(GridInterp_ZL2L  Interpolate -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/ZL.fstd -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/L.fstd -o /tmp/libgeoref.fstd -n GRID -e ZL2L -v DEBUG)
#add_test(GridInterp_ZW2L  Interpolate -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/ZW.fstd -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/L.fstd -o /tmp/libgeoref.fstd -n GRID -e ZW2L -v DEBUG)
#add_test(GridInterp_O2L   Interpolate -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/O.fstd -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/L.fstd -o /tmp/libgeoref.fstd -n GRID -e O2L -v DEBUG)
#add_test(GridInterp_M2L   Interpolate -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/M.fstd -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/L.fstd -o /tmp/libgeoref.fstd -n GRID -e M2L -v DEBUG)

#add_test(GridInterp_L2PS  Interpolate -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/PS.fstd -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/L.fstd -o /tmp/libgeoref.fstd -n GRID -e L2PS -v DEBUG)
#add_test(GridInterp_L2PSN Interpolate -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/PSN.fstd -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/L.fstd -o /tmp/libgeoref.fstd -n GRID -e L2PSN -v DEBUG)
#add_test(GridInterp_L2PSS Interpolate -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/PSS.fstd -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/L.fstd -o /tmp/libgeoref.fstd -n GRID -e L2PSS -v DEBUG)
#add_test(GridInterp_L2ZE  Interpolate -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/ZE.fstd -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/L.fstd -o /tmp/libgeoref.fstd -n GRID -e L2ZE -v DEBUG)
#add_test(GridInterp_L2ZL  Interpolate -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/ZL.fstd -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/L.fstd -o /tmp/libgeoref.fstd -n GRID -e L2ZL -v DEBUG)
#add_test(GridInterp_L2ZW  Interpolate -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/ZW.fstd -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/L.fstd -o /tmp/libgeoref.fstd -n GRID -e L2ZW -v DEBUG)
#add_test(GridInterp_L2O   Interpolate -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/O.fstd -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/L.fstd -o /tmp/libgeoref.fstd -n GRID -e L2O -v DEBUG)
#add_test(GridInterp_L2M   Interpolate -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/M.fstd -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/L.fstd -o /tmp/libgeoref.fstd -n GRID -e L2M -v DEBUG)

add_test(GridInterp_ZE2PS  Interpolate -m LINEAR -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/PS.fstd -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/HRDPS.fstd -o /tmp/libgeoref.fstd -n TT -e ZE2PS -v DEBUG)
add_test(GridInterp_ZE2PSN Interpolate -m LINEAR -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/PSN.fstd -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/HRDPS.fstd -o /tmp/libgeoref.fstd -n TT -e ZE2PSN -v DEBUG)
add_test(GridInterp_ZE2PSS Interpolate -m LINEAR -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/PSS.fstd -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/HRDPS.fstd -o /tmp/libgeoref.fstd -n TT -e ZE2PSS -v DEBUG)
add_test(GridInterp_ZE2ZE  Interpolate -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/ZE.fstd -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/HRDPS.fstd -o /tmp/libgeoref.fstd -n TT -e ZE2ZE -v DEBUG)
add_test(GridInterp_ZE2ZL  Interpolate -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/ZL.fstd -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/HRDPS.fstd -o /tmp/libgeoref.fstd -n TT -e ZE2ZL -v DEBUG)
add_test(GridInterp_ZE2ZW  Interpolate -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/ZW.fstd -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/HRDPS.fstd -o /tmp/libgeoref.fstd -n TT -e ZE2ZW -v DEBUG)
add_test(GridInterp_ZE2O   Interpolate -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/O.fstd -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/HRDPS.fstd -o /tmp/libgeoref.fstd -n TT -e ZE2O -v DEBUG)
add_test(GridInterp_ZE2M   Interpolate -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/M.fstd -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/HRDPS.fstd -o /tmp/libgeoref.fstd -n TT -e ZE2M -v DEBUG)

add_test(GridInterp_O2ZE   Interpolate -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/ZE.fstd -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/O.fstd -o /tmp/libgeoref.fstd -n GRID -e O2ZE -v DEBUG)
