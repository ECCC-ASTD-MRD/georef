import georef
import numpy as np
import os
from rmn import fst24_file, fst_record, FstDataType
from georef.cubed_sphere import decode_angle, decode_ig4


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
    SUPPORTED_GRIDS = {'A', 'B', 'E', 'G', 'H', 'L', 'N', 'Q', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '!', '#'}     

    if grtyp.upper() not in SUPPORTED_GRIDS:
        available = ", ".join(sorted(SUPPORTED_GRIDS))
        raise ValueError(f"Type de grille '{grtyp}' inconnu. Types supportés : {available}")
    
    if os.path.exists(filename):
        os.remove(filename)

    # Ouverture du fichier
    with fst24_file(filename, "RSF+R/W")  as f_fst:

        # Grille
        if grtyp == "Q":
            num_elem, num_solpts = decode_ig4(ig4)
            geo = georef.CubedSphereRef(decode_angle(ig1), decode_angle(ig2), decode_angle(ig3), num_elem, num_solpts, f_fst)
        else:
            geo = georef.GeoRef(ni, nj, grtyp, ig1, ig2, ig3, ig4, f_fst)
        lats, lons = geo.getll()
        
        print(f"GRID générée {grtyp}")

        # Calcul des données
        data = angular_dist_0(np.radians(lons), np.radians(lats)).astype(np.float32)

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
        rec.ni = geo.shape[0] 
        rec.nj = geo.shape[1] 
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


