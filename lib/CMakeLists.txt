message(STATUS "Generating libgeoref librairie")

file(GLOB PROJECT_INCLUDE_FILES *.h)
file(GLOB PROJECT_C_FILES *.c)
file(GLOB PROJECT_F_FILES f90/*.f90 )

add_library(georef ${PROJECT_INCLUDE_FILES} ${PROJECT_C_FILES} ${PROJECT_F_FILES})

set_target_properties(georef PROPERTIES 
   VERSION ${PROJECT_VERSION}
   PUBLIC_HEADER "${PROJECT_INCLUDE_FILES}"
   POSITION_INDEPENDENT_CODE  ON
   COMPILE_FLAGS "${MPI_COMPILE_FLAGS}"
   LINK_FLAGS "${MPI_LINK_FLAGS}")

add_dependencies(georef ${PROJECT_NAME}_build_info)

if (rmn_FOUND)
   target_link_libraries(georef rmn)
endif()

if (vgrid_FOUND)
   target_link_libraries(georef vgrid)
endif()

#target_link_libraries(georef ${LIBXML2_LIBRARIES})
if(GDAL_FOUND)
   target_link_libraries(georef ${GDAL_LIBRARIES})
endif()

install(TARGETS georef 
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib
   PUBLIC_HEADER DESTINATION include)
   