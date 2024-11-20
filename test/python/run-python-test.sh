this_dir=$(cd $(dirname $0) && pwd)
# Librmn with the python package must be loaded
# This python package must be loaded.  This can either be the
# one from the build directory or the one form the install directory and this
# is done by prepending ${build}/python or ${install}/lib/python to PYTHONPATH
python3.8 ${this_dir}/test/python/interpolate.py \
    --input ${ECCI_DATA_DIR}/libgeoref/grids/L.fstd \
    --grid ${ECCI_DATA_DIR}/libgeoref/grids/L.fstd \
    --output output_file.fstd \
    -n GRID \
    -e L2L
