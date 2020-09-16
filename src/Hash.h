


#define CHUNK           128
#define LOG2_CHUNK      7
#define MAX_LOG_CHUNK   12



// Previously c_gdkey2rowcol
static inline void* Hash_Get(THash *Hash,int Key) {
   return(Hash[key >> LOG2_CHUNK][key % CHUNK])
}

// Previously c_gdrowcol2key
static inline int Hash_Put(THash *Hash,void *Data, int row, int col) {
  *key = row * CHUNK + col;
}
