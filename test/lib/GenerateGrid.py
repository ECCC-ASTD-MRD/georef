import georef
import numpy as np
import os
from rmn import fst24_file, fst_record, FstDataType


# Fonction Angular Distance
def angular_dist_0(lon, lat):
    f1 = np.cos(lon) * np.sin(lat)
    f2 = np.sin(lon)
    return np.abs(np.arctan2(np.sqrt(f1**2 + f2**2), np.cos(lon) * np.cos(lat)))

# Generation des fichiers .fst
def generate_grid(lons, lats, filename, ni, nj, grtyp="Q", ig1=0, ig2=0, ig3=0, ig4=0, nomvar="DIST", etiket="TEST"):
    """
    Génère un fichier FST : grille + champ
    """        
    
    if os.path.exists(filename):
        os.remove(filename)

    # Ouverture du fichier
    with fst24_file(filename, "RSF+R/W")  as f_fst:

        # Grille
        geo = georef.GeoRef(ni, nj, grtyp, ig1, ig2, ig3, ig4, f_fst)
        lats, lons = geo.getll()
        
        print("GRID générée")

        # Calcul des données
        data = angular_dist_0(np.radians(lons), np.radians(lats)).astype(np.float32)

        #if np.isnan(data).any():
         #   raise ValueError(f"Données")

        geo.write_fst(f_fst, ig1, ig2, ig3, ig4, "my_grid")

        # Champ
        rec = fst_record()
        rec.data_type = FstDataType.FST_TYPE_REAL
        rec.data_bits = 32
        rec.pack_bits = 32
        rec.dateo = 0      
        rec.deet  = 0      
        rec.npas  = 0      
        rec.datev = 0      
        rec.data = data 
        rec.ni = ni 
        rec.nj = nj
        rec.nk = 1
        rec.nomvar = nomvar
        rec.etiket = etiket
        rec.grtyp  = grtyp
        rec.typvar = "X"
        rec.ip1 = ig1
        rec.ip2 = ig2
        rec.ip3 = ig3
        rec.ig1 = 0
        rec.ig2 = 0 
        rec.ig3 = 0
        rec.ig4 = 0 

        # Écriture du champ
        FST_REWRITE = 2
        f_fst.write(rec, FST_REWRITE)


