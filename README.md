Georeference / Grids management and interpolation library


## Environment variables

* GEOREF_DESCRIPTOR_64  : Use 64 bit precision when writing of grid descriptors, (default: 32 bit)
* GEOREF_INDEX_SIZE_HINT: Hint for size of interpolation weight index (defaul: 1024)

## Build dependencies

- CMake 3.12+
- librmn

Note: **cmake_rpn** is included as a submodule.  Please clone with the
**--recursive** flag or run **git submodule update --init --recursive** in the
git repo after having cloned.

## At CMC

Source the right file depending on the architecture you need from the env directory.
This will load the specified compiler and define the ECCI_DATA_DIR variable for the test datasets

- Example for PPP3 and skylake specific architecture:

```bash
. $ECCI_ENV/latest/ubuntu-18.04-skylake-64/intel-19.0.3.199.sh
```

- Example for XC50 on intel-19.0.5

```bash
. $ECCI_ENV/latest/sles-15-skylake-64/intel-19.0.5.281.sh
```

- Example for CMC network and gnu 7.5:

```bash
. $ECCI_ENV/latest/ubuntu-18.04-amd-64/gnu-7.5.0.sh
```

## Build and install

```bash
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=[your install path] -DWITH_OMPI=[TRUE|FALSE] -Drmn_ROOT=[rmnlib location]
make -j 4
make test
make install
```