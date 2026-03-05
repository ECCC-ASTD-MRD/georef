import numpy as np
import georef
from rmn import fst24_file, fst_record, FstDataType


# Fonction Angular Distance
def angular_dist_0(lon, lat):
    f1 = np.cos(lon) * np.sin(lat)
    f2 = np.sin(lon)
    return np.abs(np.arctan2(np.sqrt(f1**2 + f2**2), np.cos(lon) * np.cos(lat)))

def encode_cs_ig4(num_elements, num_solpts):
    return ((num_elements & 0x1ffff) << 7) | (num_solpts & 0x3f)

# Generation du fichier Q.fst
def write_complete_fst(lons, lats, filename, ni, num_elem, num_solpts, nomvar="DIST", etiket="TEST"):
    """
    Génère un fichier FST : grille + champ
    """
    nj = ni * 6
    grtyp = "Q"
    
    ig1 = 0x420000  
    ig2 = 0xa4fa00   
    ig3 = 0x660000   
    ig4 = encode_cs_ig4(num_elem, num_solpts)

    # Calcul des données
    data = angular_dist_0(np.radians(lons), np.radians(lats)).astype(np.float64)

    # Ouverture du fichier
    try:
        f_fst = fst24_file(filename, "R/W") 
    except Exception as e:
        print(f"Erreur d'ouverture de {filename} : {e}")
        return

    try:
        # Grille
        geo = georef.GeoRef(ni, nj, grtyp, ig1, ig2, ig3, ig4, f_fst)
        geo.write_fst(f_fst, ig1, ig2, ig3, ig4, "my_grid")
        print("GRID généré")

        # Champ
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
        FST_SKIP = 0 
        result = f_fst.write(rec, FST_SKIP)
        
        if result is not None and result <= 0:
            print(f"Erreur d'écriture du champ {nomvar}")
        else:
            print(f"Record {nomvar} généré")

    except Exception as e:
        print(f"Erreur durant la génération du fichier : {e}")
    
    finally:
        f_fst.close()
        print(f"Fichier {filename} généré")

if __name__ == "__main__":
    write_complete_fst(np.array([0, 45]), np.array([0, 10]),"Q_python_grid_field.fst", 180, 36, 5)