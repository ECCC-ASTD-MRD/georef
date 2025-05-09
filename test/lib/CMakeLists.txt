if (EC_INIT_DONE LESS 2)
    include(CTest)
    add_custom_target(check COMMAND CTEST_OUTPUT_ON_FAILURE=true ${CMAKE_CTEST_COMMAND})

    if(BUILD_TESTING)
        if ("$ENV{ECCI_DATA_DIR}" STREQUAL "")
            message(FATAL_ERROR "The ECCI_DATA_DIR environment variable isn't defined!\nTests won't be able to find the required data files.\nPlease define that variable or set BUILD_TESTING=NO.")
        endif()

        add_link_options(-lpthread)

        IF(DEFINED ENV{INTELCOMP_HOME})
            add_link_options(-Wl,-rpath $ENV{INTELCOMP_HOME}/lib/intel64_lin) 
            add_link_options(-lintlc -lifcore -lifport)
        endif()

        add_executable(TestF90 EXCLUDE_FROM_ALL test.F90)
        add_dependencies(check TestF90)
        target_link_libraries(TestF90 georef-shared rmn::rmn m)

        add_test(TestF90 TestF90 $ENV{ECCI_DATA_DIR}/georef/grids/L.fstd $ENV{ECCI_DATA_DIR}/georef/grids/PS.fstd)

        add_executable(Interpolate EXCLUDE_FROM_ALL Interpolate.c)
        add_dependencies(check Interpolate)

        target_link_libraries(Interpolate georef-ompi-shared rmn::rmn-ompi m)
        set_target_properties(Interpolate PROPERTIES VERSION ${PROJECT_VERSION})

        file(REMOVE libgeoref.fstd)

        add_test(GridInterp_L2L   Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/L.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/L.fstd -o ./georef.fstd -n GRID -e L2L -v DEBUG)
        add_test(GridInterp_PS2L  Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/L.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/PS.fstd -o ./georef.fstd -n GRID -e PS2L -v DEBUG)
        add_test(GridInterp_PSN2L Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/L.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/PSN.fstd -o ./georef.fstd -n GRID -e PSN2L -v DEBUG)
        add_test(GridInterp_PSS2L Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/L.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/PSS.fstd -o ./georef.fstd -n GRID -e PSS2L -v DEBUG)
        add_test(GridInterp_ZE2L  Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/L.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/ZE.fstd -o ./georef.fstd -n GRID -e ZE2L -v DEBUG)
        add_test(GridInterp_ZL2L  Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/L.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/ZL.fstd -o ./georef.fstd -n GRID -e ZL2L -v DEBUG)
        add_test(GridInterp_O2L   Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/L.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/O.fstd -o ./georef.fstd -n GRID -e O2L -v DEBUG)
        add_test(GridInterp_M2L   Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/L.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/M.fstd -o ./georef.fstd -n GRID -e M2L -v DEBUG)

        add_test(GridInterp_L2PS  Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/PS.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/L.fstd -o ./georef.fstd -n GRID -e L2PS -v DEBUG)
        add_test(GridInterp_L2PSN Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/PSN.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/L.fstd -o ./georef.fstd -n GRID -e L2PSN -v DEBUG)
        add_test(GridInterp_L2PSS Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/PSS.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/L.fstd -o ./georef.fstd -n GRID -e L2PSS -v DEBUG)
        add_test(GridInterp_L2ZE  Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/ZE.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/L.fstd -o ./georef.fstd -n GRID -e L2ZE -v DEBUG)
        add_test(GridInterp_L2ZL  Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/ZL.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/L.fstd -o ./georef.fstd -n GRID -e L2ZL -v DEBUG)
        add_test(GridInterp_L2O   Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/O.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/L.fstd -o ./georef.fstd -n GRID -e L2O -v DEBUG)
        add_test(GridInterp_L2M   Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/M.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/L.fstd -o ./georef.fstd -n GRID -e L2M -v DEBUG)

        add_test(GridInterp_LZE2PS  Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/PS.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/HRDPS.fstd -o ./georef.fstd -n TT -e LZE2PS -v DEBUG)
        add_test(GridInterp_LZE2PSN Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/PSN.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/HRDPS.fstd -o ./georef.fstd -n TT -e LZE2PSN -v DEBUG)
        add_test(GridInterp_LZE2PSS Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/PSS.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/HRDPS.fstd -o ./georef.fstd -n TT -e LZE2PSS -v DEBUG)
        add_test(GridInterp_LZE2ZE  Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/ZE.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/HRDPS.fstd -o ./georef.fstd -n TT -e LZE2ZE -v DEBUG)
        add_test(GridInterp_LZE2ZL  Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/ZL.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/HRDPS.fstd -o ./georef.fstd -n TT -e LZE2ZL -v DEBUG)
        add_test(GridInterp_LZE2O   Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/O.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/HRDPS.fstd -o ./georef.fstd -n TT -e LZE2O -v DEBUG)
        add_test(GridInterp_LZE2M   Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/M.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/HRDPS.fstd -o ./georef.fstd -n TT -e LZE2M -v DEBUG)

        add_test(GridInterp_O2ZE   Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/ZE.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/O.fstd -o ./georef.fstd -n GRID -e O2ZE -v DEBUG)
        add_test(GridInterp_LZE2U  Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/U.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/HRDPS.fstd -o ./georef.fstd -n TT -e LZE2U -v DEBUG)
        add_test(GridInterp_O2U    Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/U.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/O.fstd -o ./georef.fstd -n GRID -e O2U -v DEBUG)
        add_test(GridInterp_U2ZE   Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/ZE.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/U.fstd -o ./georef.fstd -n GRID -e U2ZE -v DEBUG)
        add_test(GridInterp_U2O    Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/O.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/U.fstd -o ./georef.fstd -n GRID -e U2O -v DEBUG)

        if(GDAL_FOUND)
            add_test(GridInterp_ZW2L  Interpolate -i $ENV{ECCI_DATA_DIR}/georef/grids/ZW.fstd -g $ENV{ECCI_DATA_DIR}/georef/grids/L.fstd -o ./georef.fstd -n GRID -e ZW2L -v DEBUG)
            add_test(GridInterp_L2ZW  Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/ZW.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/L.fstd -o ./georef.fstd -n GRID -e L2ZW -v DEBUG)
            add_test(GridInterp_LZE2ZW  Interpolate -g $ENV{ECCI_DATA_DIR}/georef/grids/ZW.fstd -i $ENV{ECCI_DATA_DIR}/georef/grids/HRDPS.fstd -o ./georef.fstd -n TT -e LZE2ZW -v DEBUG)
        endif()

    endif()
endif()
