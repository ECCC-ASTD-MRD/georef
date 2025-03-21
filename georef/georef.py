    def write_fst(self, file: 'FSTFile') -> bool:
        """Write geographic set data to an FST file.
        
        This method wraps the libgeoref function GeoRef_SetWriteFST() found in src/GeoRef.c.
        
        Args:
            file (FSTFile): FST file to write to
            
        Returns:
            bool: True if successful, False otherwise
            
        Raises:
            GeoRefError: If writing fails or GeoSet is uninitialized
            
        Note:
            The underlying C function returns:
            - 0 for failure (NULL references or write error)
            - 1 for successful write operation
        """
        if not self._ptr:
            raise GeoRefError("Cannot write uninitialized GeoSet")
            
        val = _geoset_writefst(self._ptr, file._ptr)
        return val == 1 