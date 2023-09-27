message(STATUS "Generating libgeoref librairie")

file(GLOB PROJECT_INCLUDE_FILES *.h)
file(GLOB PROJECT_C_FILES *.c)
file(GLOB PROJECT_F_FILES f90/*.f90 )

add_library(georef ${PROJECT_INCLUDE_FILES} ${PROJECT_C_FILES} ${PROJECT_F_FILES})
target_include_directories(georef PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_include_directories(georef INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>)

set_target_properties(georef PROPERTIES 
   VERSION ${PROJECT_VERSION}
   PUBLIC_HEADER "${PROJECT_INCLUDE_FILES}"
   INTERFACE_COMPILE_DEFINITIONS HAVE_GEOREF="${PROJECT_VERSION}"
   POSITION_INDEPENDENT_CODE  ON
   COMPILE_FLAGS "${MPI_COMPILE_FLAGS}"
   LINK_FLAGS "${MPI_LINK_FLAGS}")

add_library(georef::georef ALIAS georef)

add_dependencies(georef ${PROJECT_NAME}_build_info)

if (rmn_FOUND)
   target_link_libraries(georef rmn::rmn-ompi)
endif()

#if (vgrid_FOUND)
#   target_link_libraries(georef vgrid::vgrid)
#endif()

#target_link_libraries(georef ${LIBXML2_LIBRARIES})
if(GDAL_FOUND)
   target_link_libraries(georef GDAL::GDAL)
endif()

list(APPEND shared_targets georef)
#list(APPEND static_targets rmn-ompi)

install(TARGETS ${shared_targets}
   EXPORT georef-shared-targets
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib
   INCLUDES DESTINATION include
   PUBLIC_HEADER DESTINATION include)
   