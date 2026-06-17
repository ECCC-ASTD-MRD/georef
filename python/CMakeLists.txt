
# For the version of the 'georef' python package that is in the build directory
# it is more convenient for it to contain links to the python files in the source
# otherwise we would have to re-run CMake everytime we change them.  At install
# time of course, we copy the python files.
file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/georef ${CMAKE_CURRENT_BINARY_DIR}/georef SYMBOLIC)

install(DIRECTORY georef DESTINATION lib/python FILES_MATCHING PATTERN *.py)
