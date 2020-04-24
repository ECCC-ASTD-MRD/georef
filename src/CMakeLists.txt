message(STATUS "Generating libgeoref librairie")


file(GLOB PROJECT_INCLUDE_FILES *.h)
file(GLOB PROJECT_C_FILES *.c)
file(GLOB PROJECT_F_FILES *.f90 *.f)
add_library(georef SHARED ${PROJECT_INCLUDE_FILES} ${PROJECT_C_FILES} ${PROJECT_F_FILES})

set_target_properties(georef PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(georef PROPERTIES PUBLIC_HEADER "${PROJECT_INCLUDE_FILES}")
target_link_libraries(georef ${RMN_LIBRARIES})

add_dependencies(georef timestamp)

install(TARGETS georef
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib
   PUBLIC_HEADER DESTINATION include)
