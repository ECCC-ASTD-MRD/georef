import numpy as np
from rmn import fst24_file, fst_record, FstDataType


# Fonction Angular Distance
def angular_dist_0(lon, lat):
    f1 = np.cos(lon) * np.sin(lat)
    f2 = np.sin(lon)
    return np.abs(np.arctan2(np.sqrt(f1**2 + f2**2), np.cos(lon) * np.cos(lat)))

# Generation du fichier Q.fst
def write_test_fst(lons, lats, ni):

    data = angular_dist_0(np.radians(lons), np.radians(lats)).astype(np.float64)
    
    f_fst = fst24_file("Q_python.fst", "R/W") 
    rec = fst_record()

    rec.data_type = FstDataType.FST_TYPE_REAL
    rec.data_bits = 64
    rec.pack_bits = 32
    rec.dateo = 0      
    rec.deet  = 0      
    rec.npas  = 0      
    rec.datev = 0    
    rec.data = data 
    rec.ni = ni
    rec.nj = ni * 6
    rec.nk = 1
    rec.nomvar = "DIST"
    rec.etiket = "TEST_FIELD"
    rec.grtyp = "Q"
    rec.typvar = "X"  
    rec.ip1 = 0
    rec.ip2 = 0
    rec.ip3 = 0
    rec.ig1 = 0
    rec.ig2 = 0
    rec.ig3 = 0
    rec.ig4 = 0

    try:
        FST_SKIP = 0 
        result = f_fst.write(rec, FST_SKIP)

        if result is not None and result <= 0:
            print(f"Erreur d'écriture : code {result}")
        else:
            print("Fichier généré")
            
    except Exception as e:
        print(f"Erreur lors de l'appel à write : {e}")

if __name__ == "__main__":
    write_test_fst(np.array([0, 45]), np.array([0, 10]), 180)