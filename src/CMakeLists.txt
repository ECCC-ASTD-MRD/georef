message(STATUS "Generating libgeoref librairie")

file(GLOB PROJECT_INCLUDE_FILES include/*.h include/georef/*.h)
file(GLOB PROJECT_C_FILES *.c)
file(GLOB PROJECT_F_FILES *.F90 F90/*.F90)

set(shared_targets georef-shared georef-ompi-shared)
set(static_targets georef-static georef-ompi-static)
set(targets ${static_targets} ${shared_targets})

add_library(georef-static STATIC ${PROJECT_INCLUDE_FILES} ${PROJECT_C_FILES} ${PROJECT_F_FILES})
add_library(georef-shared SHARED $<TARGET_OBJECTS:georef-static>)
add_library(georef-ompi-static STATIC $<TARGET_OBJECTS:georef-static>)
add_library(georef-ompi-shared SHARED $<TARGET_OBJECTS:georef-static>)

add_library(georef::georef ALIAS georef-static)

foreach(LINK_TYPE static shared)
    foreach(OMPI "" -ompi)
        set_target_properties(georef${OMPI}-${LINK_TYPE} PROPERTIES
            VERSION ${PROJECT_VERSION}
            PUBLIC_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/include/GeoRef.h
            INTERFACE_COMPILE_DEFINITIONS HAVE_GEOREF="${PROJECT_VERSION}"
            Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include
            OUTPUT_NAME georef${OMPI}
            POSITION_INDEPENDENT_CODE ON)

        target_include_directories(georef${OMPI}-${LINK_TYPE} PUBLIC 
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/georef>
                $<INSTALL_INTERFACE:include>
                $<INSTALL_INTERFACE:include/georef>)
        #        target_include_directories(georef${OMPI}-${LINK_TYPE} INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>)
                
        target_compile_definitions(georef${OMPI}-${LINK_TYPE} PRIVATE HAVE_GPC)
        add_dependencies(georef${OMPI}-${LINK_TYPE} ${PROJECT_NAME}_build_info)
        add_library(georef::georef${OMPI}-${LINK_TYPE} ALIAS georef${OMPI}-${LINK_TYPE})
            
        if(OpenMP_FOUND)
            target_link_libraries(georef${OMPI}-${LINK_TYPE} PUBLIC OpenMP::OpenMP_C OpenMP::OpenMP_Fortran)
        endif()

        if(rmn_FOUND)
            target_link_libraries(georef${OMPI}-${LINK_TYPE} PUBLIC rmn::rmn${OMPI}-${LINK_TYPE})
        endif()

        if(vgrid_FOUND)
            target_link_libraries(georef${OMPI}-${LINK_TYPE} PUBLIC vgrid::vgrid)
        endif()

        if(GDAL_FOUND)
            target_link_libraries(georef${OMPI}-${LINK_TYPE} PUBLIC GDAL::GDAL)
        endif()

    endforeach()
    install(TARGETS ${${LINK_TYPE}_targets} EXPORT ${PROJECT_NAME}-${LINK_TYPE}-targets)
endforeach()

install(DIRECTORY include/georef DESTINATION include)
install(FILES include/GeoRef.inc TYPE INCLUDE)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/georef_mod.mod TYPE INCLUDE)

# Symbolic link for archive numbered version
add_custom_command(TARGET ${PROJECT_NAME}-static POST_BUILD
    COMMAND ln -sfv lib${PROJECT_NAME}.a lib${PROJECT_NAME}_${PROJECT_VERSION}.a
    BYPRODUCTS lib${PROJECT_NAME}_${PROJECT_VERSION}.a)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib${PROJECT_NAME}_${PROJECT_VERSION}.a DESTINATION lib)
