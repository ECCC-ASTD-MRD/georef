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
- librmn
- GDAL (optional)

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
