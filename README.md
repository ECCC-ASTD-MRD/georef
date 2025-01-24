# Description

Georeference / Grids management and interpolation library

Georef builds on the ezscint package by adding:
* More interpolation methods (conservative, average, geometric, rasterization)
* More projections through PROJ4
* Manage tripole grids ('O'), triangular meshes ('M') and cube sphere ('C')
* All coordinates transformations are done at double precision
* Can save and read interpolation indexes and weights
* Thread safety

## Environment variables

* GEOREF_PRESERVE       : How many of the first georef created will be kept preserved in cache (default: 10)
* GEOREF_DESCRIPTOR_64  : Use 64 bit precision when writing of grid descriptors, (default: 32 bit)
* GEOREF_INDEX_SIZE_HINT: Hint for size of interpolation weight index (default: 1024)

## Structures

* TGeoOptions : Defines various intepolation parameters
* TGeoRef     : Geo-Reference definition
* TGeoSet     : Links 2 Geo-Reference for interpolation. Contains cached values and index/weights to be reused on each intepolation betwee the 2 geo-reference

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
   GeoRef_Options.Interp=IR_LINEAR;
   GeoRef_Options.NoData=0.0;
   GeoRef_Options.Extrap=ER_VALUE;
   GeoRef_Options.ExtrapValue=-999.0;

   ...

   // Loop on records
   fst_query* query = fst24_new_query(fin,&crit,NULL);
   while(fst24_read_next(query,&record)>0 && n++<10) {

      // Create record georeference 
      refin=GeoRef_Create(record.ni,record.nj,record.grtyp,record.ig1,record.ig2,record.ig3,record.ig4,(fst_file*)record.file);

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