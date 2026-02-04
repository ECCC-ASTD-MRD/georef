
# For the version of the 'rmn' python package that is in the build directory
# it is more convenient for it to contain links to the python files in the source
# otherwise we would have to re-run CMake everytime we change them.  At install
# time of course, we copy the python files.
#
set(python_files
    __init__.py
    _georef_c_bindings.py
    constants.py
    enums.py
    georef.py
    shared_lib.py
    structs.py
)
message(STATUS "python_files=${python_files}")
foreach(f IN LISTS python_files)
    message(STATUS "f=${f}")
    file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/${f} ${CMAKE_CURRENT_BINARY_DIR}/${f} SYMBOLIC)
endforeach()

install(FILES ${python_files} DESTINATION lib/python/georef)

configure_file(dev-setup-python.in.sh ${CMAKE_BINARY_DIR}/dev-setup-python.sh @ONLY)
add_test(NAME run-python-interpolate COMMAND bash -c "source ${CMAKE_BINARY_DIR}/dev-setup-python.sh; ${CMAKE_SOURCE_DIR}/test/python/run-python-test.sh")

