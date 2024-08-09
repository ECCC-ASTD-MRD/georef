message(STATUS "Generating libgeoref librairie")

file(GLOB PROJECT_INCLUDE_FILES *.h)
file(GLOB PROJECT_C_FILES *.c)
file(GLOB PROJECT_F_FILES f90/*.f90 )

set(targets georef-static georef-shared)
set(shared_targets georef-shared)
set(static_targets georef-static)

add_library(georef-static STATIC ${PROJECT_INCLUDE_FILES} ${PROJECT_C_FILES} ${PROJECT_F_FILES})
target_include_directories(georef-static PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_include_directories(georef-static INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>)

add_library(georef-shared SHARED $<TARGET_OBJECTS:georef-static>)

add_dependencies(georef-static ${PROJECT_NAME}_build_info)
set_target_properties(georef-static PROPERTIES 
   VERSION ${PROJECT_VERSION}
   PUBLIC_HEADER "${PROJECT_INCLUDE_FILES}"
   INTERFACE_COMPILE_DEFINITIONS HAVE_GEOREF="${PROJECT_VERSION}"
   Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include
   OUTPUT_NAME georef
   POSITION_INDEPENDENT_CODE  ON
   COMPILE_FLAGS "${MPI_COMPILE_FLAGS}"
   LINK_FLAGS "${MPI_LINK_FLAGS}")

set_target_properties(georef-shared PROPERTIES 
      VERSION ${PROJECT_VERSION}
      PUBLIC_HEADER "${PROJECT_INCLUDE_FILES}"
      INTERFACE_COMPILE_DEFINITIONS HAVE_GEOREF="${PROJECT_VERSION}"
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include
      OUTPUT_NAME georef
      POSITION_INDEPENDENT_CODE  ON
      COMPILE_FLAGS "${MPI_COMPILE_FLAGS}"
      LINK_FLAGS "${MPI_LINK_FLAGS}")
   
add_library(georef::georef-static ALIAS georef-static)
add_library(georef::georef-shared ALIAS georef-shared)
add_library(georef::georef ALIAS georef-static)

if (rmn_FOUND)
   target_link_libraries(georef-static rmn::rmn-ompi)
   target_link_libraries(georef-shared rmn::rmn-ompi)
endif()

if (vgrid_FOUND)
   target_link_libraries(georef-static vgrid::vgrid)
   target_link_libraries(georef-shared vgrid::vgrid)
endif()

#target_link_libraries(georef ${LIBXML2_LIBRARIES})
if(GDAL_FOUND)
   target_link_libraries(georef-static GDAL::GDAL)
   target_link_libraries(georef-shared GDAL::GDAL)
endif()

install(TARGETS ${shared_targets}
   EXPORT ${PROJECT_NAME}-shared-targets
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib
   INCLUDES DESTINATION include
   PUBLIC_HEADER DESTINATION include)
   
install(TARGETS ${static_targets}
   EXPORT ${PROJECT_NAME}-static-targets
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib
   INCLUDES DESTINATION include
   PUBLIC_HEADER DESTINATION include)  

# Symbolic link for archive numbered version
add_custom_command(TARGET ${PROJECT_NAME}-static POST_BUILD
   COMMAND ln -sfv lib${PROJECT_NAME}.a lib${PROJECT_NAME}_${PROJECT_VERSION}.a
   BYPRODUCTS lib${PROJECT_NAME}_${PROJECT_VERSION}.a)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib${PROJECT_NAME}_${PROJECT_VERSION}.a DESTINATION lib)
