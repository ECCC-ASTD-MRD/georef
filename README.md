# Georeference / Grid management and interpolation library

## Description

This package provides C/Fortran/Python functions to project latitude/longitude coordinates to x/y grid space coordinates and vice-versa, get values at latitude/longitude or x/y coordinate values and interpolation from one grid to another.

Georef builds on the ezscint package by adding:
* More interpolation methods (conservative, average, geometric, rasterization)
* More projections through PROJ4 ('W'), Needs to build with optional GDAL library
* Manage tripole grids ('O'), triangular meshes ('M') and cube sphere ('C')
* All coordinates transformations are done at double precision
* Can read and save grid records in 64 bits (export GEOREF_DESCRIPTOR_64=TRUE)
* Can save and read interpolation indexes and weights
* Many levels of internal caches to speed-up transformations
* Functions are re-entrant (thread safety)

## Environment variables

* GEOREF_PRESERVE       : How many of the first georef created will be kept preserved in cache (default: 10)
* GEOREF_MAXSET         : Maximum number of interpolation geoset to keep in cache (default: 64)
* GEOREF_DESCRIPTOR_64  : Use 64 bit precision when writing of grid descriptors, (default: 32 bit)
* GEOREF_INDEX_SIZE_HINT: Hint for size of interpolation weight index (default: 1024)

## Data structures

* TGeoOptions : Defines various intepolation parameters. There is a global instance to TGeoOption *GeoRef_Options that should only be used in non threaded environment. Use a per thread instance of this type as parametre to functions needing them for threaded environment.
* TGeoRef     : Geo-Reference definition
* TGeoSet     : Links 2 Geo-Reference for interpolation. Contains cached values and index/weights to be reused on each intepolation between the 2 geo-reference

## Example usage:

```C
   // Define structures
   TGeoSet    *gset=NULL;
   TGeoRef    *refin=NULL,*refout=NULL;
   fst_record  grid=default_fst_record,record=default_fst_record;
   fst_file   *fin,*fout;

   ...

   // Create output georeference from a read record
   refout=GeoRef_Create(grid.ni,grid.nj,grid.grtyp,grid.ig1,grid.ig2,grid.ig3,grid.ig4,(fst_file*)grid.file);
   if(refout == NULL){
       App_Log(APP_ERROR, "Problem creating TGeoRef object\n");
       return(FALSE);
   }

   ...

   // Define interpolation options
   GeoRef_Options.Interp=IR_LINEAR;    // Linear interpolation
   GeoRef_Options.NoData=-999.0;       // No data value
   GeoRef_Options.Extrap=ER_VALUE;     // Use nodata value when extrapolating

   ...

   // Loop on records
   fst_query* query = fst24_new_query(fin,&crit,NULL);
   while(fst24_read_next(query,&record)>0 && n++<10) {

      // Create record georeference 
      refin=GeoRef_CreateFromRecord(&record);

      // Proceed with interpolation
      if (!GeoRef_Interp(refout,refin,&GeoRef_Options,grid.data,record.data)) {
         App_Log(APP_ERROR,"Interpolation problem");
         return(FALSE);
      }   
    
      // Write results
	  fst24_record_copy_metadata(&grid,&record,FST_META_TIME|FST_META_INFO);
      if (fst24_write(fout, &grid, FST_NO) < 0) {
         App_Log(APP_ERROR, "Unable to write record\n");
         return FALSE;
      }
   }

   ...

   // Get index/weights and write to file 
   gset=GeoRef_SetGet(refout,refin,NULL);
   if (GeoRef_SetHasIndex(gset)) {
      App_Log(APP_DEBUG,"Saving index containing %i items\n",gset->IndexSize);
      
      if (!GeoRef_SetWrite(gset,fout)){
         return(0);
      }
   }

   ...    
```

```fortran

   type(fst_file)   :: file1,file2,fileout
   type(fst_record) :: record1,record2
   type(georef)     :: gref1,gref2

   logical :: success
   integer :: len, status
   character(len=4096) :: argument
   real(kind = real32), dimension(:), pointer :: data_array1, data_array2

   ! Define interpolation options
   georef_options%interp=IR_LINEAR;    ! Linear interpolation
   georef_options%nodata=-999.0;       ! No data value
   georef_options%extrap=ER_VALUE;     ! Use nodata value when extrapolating

   ! Read record from first file
   call get_command_argument(1,argument,len,status)
   success = file1%open(trim(argument))
   if (.not. success) then
      call App_log(APP_ERROR, 'Unable to open FST file')
      call exit(-1)
   end if

   ! Read record from second file
   call get_command_argument(2,argument,len,status)
   success = file2%open(trim(argument))
   if (.not. success) then
      call App_log(APP_ERROR, 'Unable to open FST file')
      call exit(-1)
   end if

   ! Create geo-references
   success=file1%read(record1,NOMVAR="GRID")
   success=gref1%fromrecord(record1)

   success=file2%read(record2,NOMVAR="GRID")
   success=gref2%fromrecord(record2)

   ! Interpolation from gref2 to gref1
   call record1%get_data_array(data_array1)
   call record2%get_data_array(data_array2)

   len=gref1%interp(gref2,data_array1,data_array2)
   success=fileout%write(record1)

   success = file1%close()
   success = file2%close()
   success = fileout%close()

```

# Compilation

## Build dependencies

- CMake 3.12+
- cmake\_rpn
- librmn
- vgrid
- GDAL (optional)

## At CMC

load the right [code-tools](https://gitlab.science.gc.ca/RPN-SI/code-tools) environment depending on the architecture you need from the env directory.
This will load the specified compiler, cmake_rpn and the specific compiler parameters per platform

- Example for PPP6/SC6 and icelake specific architecture:

```bash
. r.load.dot mrd/rpn/code-tools/latest/env/rhel-8-icelake-64@intelonapi-2025.1.0
```

- Example for generic architecture on ppp6/SC6

```bash
   r.load.dot mrd/rpn/code-tools/latest/env/rhel-8-amd64-64@intelonapi-2025.1.0
```

- Example for GNU on any architecture:

```bash
   r.load.dot mrd/rpn/code-tools/latest/env/gnu
```

## Outside CMC

We need to intall each dependency and make them available to the project.

### cmake\_rpn

This project and all other RPN-SI projects are built with CMake and use
[`cmake_rpn`](https://github.com/ECCC-ASTD-MRD/cmake_rpn). There are two ways to get it.

- Install it once with
  ```
  git clone https://github.com/ECCC-ASTD-MRD/cmake_rpn
  ```
  and set the `EC_CMAKE_MODULE_PATH` to be `<location-of-cmake_rpn>/modules`.
  All RPN-SI projects look at this variable.

- Use the included submodule
  ```bash
  # From the root of this projec
  git submodule update --init --remote cmake_rpn
  ```
  *NOTE* We will not be updating the `cmake_rpn` submodule as it has become
  stable enough that using the submodule is not worth it anymore.  However we
  do update it with new features and using `--remote` will get the latest
  version rather than the one registered to the current commit of `georef`.

### Compiled projects

The RPN-SI projects which are dependencies must be installed somewhere according
to their instructions and made findable by CMake via the
[`CMAKE_PREFIX_PATH` environment variable](https://cmake.org/cmake/help/latest/envvar/CMAKE_PREFIX_PATH.html#envvar:CMAKE_PREFIX_PATH)

They can be found at
- [librmn](https://github.com/ECCC-ASTD-MRD/librmn)
- [vgrid](https://github.com/ECCC-ASTD-MRD/vgrid)

```bash
export CMAKE_PREFIX_PATH=<your librmn install path>:${CMAKE_PREFIX_PATH}
export CMAKE_PREFIX_PATH=<your vgrid install path>:${CMAKE_PREFIX_PATH}
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

Note: If using [`CMAKE_PREFIX_PATH` environment variable](https://cmake.org/cmake/help/latest/envvar/CMAKE_PREFIX_PATH.html#envvar:CMAKE_PREFIX_PATH)
the `-Drmn_ROOT` option is not necessary in the CMake command.

Once installed, this package can be made available to other CMake projects the
same way its dependencies were made available to it via the
[`CMAKE_PREFIX_PATH` environment variable](https://cmake.org/cmake/help/latest/envvar/CMAKE_PREFIX_PATH.html#envvar:CMAKE_PREFIX_PATH)

```bash
export CMAKE_PREFIX_PATH=<your install path>:${CMAKE_PREFIX_PATH}
```
