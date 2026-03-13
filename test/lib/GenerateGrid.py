import georef
from georef.cubed_sphere import encodeig4
import numpy as np
import os
from rmn import fst24_file, fst_record, FstDataType


# Fonction Angular Distance
def angular_dist_0(lon, lat):
    f1 = np.cos(lon) * np.sin(lat)
    f2 = np.sin(lon)
    return np.abs(np.arctan2(np.sqrt(f1**2 + f2**2), np.cos(lon) * np.cos(lat)))

# Generation des fichiers .fst
def generate_grid(lons, lats, filename, ni, nj=None, num_elem=0, num_solpts=0, grtyp="Q", ig1=0, ig2=0, ig3=0, ig4=0, nomvar="DIST", etiket="TEST"):
    """
    Génère un fichier FST : grille + champ
    """        
    
    if os.path.exists(filename):
        os.remove(filename)
    
    if grtyp == "Q":
        actual_nj = ni * 6
        ig4 = encodeig4(num_elem, num_solpts) 
    else:
        actual_nj = nj if nj is not None else ni // 2

    # Ouverture du fichier
    try:
        f_fst = fst24_file(filename, "RSF+R/W") 

        # Grille
        geo = georef.GeoRef(ni, actual_nj, grtyp, ig1, ig2, ig3, ig4, f_fst)
        lats, lons = geo.getll()
        
        print("GRID généré")

        # Calcul des données
        data = angular_dist_0(np.radians(lons), np.radians(lats)).astype(np.float64)

        geo.write_fst(f_fst, ig1, ig2, ig3, ig4, "my_grid")

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
        rec.nj = actual_nj
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

    except Exception as e:
        print(f"Erreur durant la génération du fichier : {e}")
    
    finally:
        f_fst.close()
        print(f"Fichier {filename} généré")



#if __name__ == "__main__":
    # generate_grid(lons, lats, filename, ni, nj=None, num_elem=0, num_solpts=0, grtyp="Q", ig1=0, ig2=0, nomvar="DIST", etiket="TEST")
 #   generate_grid(np.array([0, 45]), np.array([0, 10]),"Grid_Q.fst", 180, num_elem=36, num_solpts=5, grtyp="Q", ig1=0x420000, ig2=0xa4fa00, ig3=0x660000)
  #  generate_grid(np.array([0, 45]), np.array([0, 10]),"Grid_A.fst", 180, grtyp="A", ig1=1)
   # generate_grid(np.array([0, 45]), np.array([0, 10]),"Grid_B.fst", 180, grtyp="B", ig2=1)
    #generate_grid(np.array([0, 45]), np.array([0, 10]),"Grid_G.fst", 180, grtyp="G", ig1=1, ig2=1)

