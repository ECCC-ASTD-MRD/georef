message(STATUS "Generating libgeoref librairie")

file(GLOB PROJECT_INCLUDE_FILES *.h)
file(GLOB PROJECT_C_FILES *.c)
file(GLOB PROJECT_F_FILES f90/*.f90 )

add_library(georef SHARED ${PROJECT_INCLUDE_FILES} ${PROJECT_C_FILES} ${PROJECT_F_FILES})

set_target_properties(georef$ENV{OMPI} PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(georef$ENV{OMPI} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(georef$ENV{OMPI} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
set_target_properties(georef$ENV{OMPI} PROPERTIES PUBLIC_HEADER "${PROJECT_INCLUDE_FILES}")

add_dependencies(georef$ENV{OMPI} build_info)

target_link_libraries(georef$ENV{OMPI} ${RMN_LIBRARIES} ${VGRID_LIBRARIES})
#target_link_libraries(georef$ENV{OMPI} ${LIBXML2_LIBRARIES})
if(GDAL_FOUND)
   target_link_libraries(georef$ENV{OMPI} ${GDAL_LIBRARIES})
endif()

install(TARGETS georef$ENV{OMPI} 
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib
   PUBLIC_HEADER DESTINATION include)
