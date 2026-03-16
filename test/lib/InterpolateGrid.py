import georef
from georef.cubed_sphere import encodeig4
import numpy as np
import os
from rmn import fst24_file, fst_record, FstDataType
from GenerateGrid import generate_grid


def interpolate_grid(src_file, dest_file, dest_params):
    """
    Interpole les données d'un fichier FST vers une nouvelle grille.
    """

    if os.path.exists(dest_file):
        os.remove(dest_file)

    # Ouverture du fichier source
    try:
        f_src = fst24_file(src_file, "R")
        # 
        query = f_src.new_query() 
        for rec in query:
            rec_src = rec
            d = rec_src.data

        if rec_src is None:
            raise ValueError(f"Champ DIST introuvable dans {src_file}")
        
        src_geo = georef.GeoRef(rec_src.ni, rec_src.nj, rec_src.grtyp, rec_src.ip1, rec_src.ip2, rec_src.ip3, rec_src.ig4, f_src)
        data_src = rec_src.data.astype(np.float32)

    finally:
        f_src.close()

    # Création du fichier destination
    # On reprend la même logique que generate_grid en mettant les records de l'interpretation dans le fichier
    try:
        f_dest = fst24_file(dest_file, "RSF+R/W")
        
        # Calcul des paramètres destination 
        d = dest_params
        if d['grtyp'] == "Q":
            d_nj = d['ni'] * 6
            d_ig4 = encodeig4(d.get('num_elem', 0), d.get('num_solpts', 0))
        else:
            d_nj = d.get('nj', d['ni'] // 2)
            d_ig4 = d.get('ig4', 0)

        dest_geo = georef.GeoRef(d['ni'], d_nj, d['grtyp'], d.get('ig1', 0), d.get('ig2', 0), d.get('ig3', 0), d_ig4, f_dest)

        dest_geo.write_fst(f_dest, d.get('ig1', 0), d.get('ig2', 0), d.get('ig3', 0), d_ig4, "dest_grid")

        # Interpolation
        data_interp = dest_geo.interp(src_geo, data_src)

        # Champ interpolé
        rec_dest = fst_record()
        rec_dest.data_type = FstDataType.FST_TYPE_REAL
        rec_dest.data_bits = 32
        rec_dest.pack_bits = 32
        rec_dest.dateo = 0      
        rec_dest.deet  = 0      
        rec_dest.npas  = 0      
        rec_dest.datev = 0
        rec_dest.data = data_interp
        rec_dest.ni = d['ni']
        rec_dest.nj = d_nj
        rec_dest.nk = 1
        rec_dest.nomvar = "DIST"
        rec_dest.etiket = "INTERP"
        rec_dest.grtyp = d['grtyp']
        rec_dest.typvar = "X"
        rec_dest.ip1 = d.get('ig1', 0)
        rec_dest.ip2 = d.get('ig2', 0)
        rec_dest.ip3 = d.get('ig3', 0)
        rec_dest.ig1 = 0
        rec_dest.ig2 = 0 
        rec_dest.ig3 = 0
        rec_dest.ig4 = 0

        FST_REWRITE = 2
        f_dest.write(rec_dest, FST_REWRITE)
        print(f"Succès : Fichier interpolé créé -> {dest_file}")

    finally:
        f_dest.close()


if __name__ == "__main__":
    # Dictionnaire des fichiers sources à créer
    grids_config = [
        {
            "lons": np.array([0, 45]), "lats": np.array([0, 10]), "grtyp": "A", "ni": 180, "ig1": 1, "filename": "Grid_A.fst", "label": "Lat-Lon Equidistante"
        },
        {
            "lons": np.array([0, 45]), "lats": np.array([0, 10]), "grtyp": "B", "ni": 180, "ig2": 1, "filename": "Grid_B.fst", "label": "Lat-Lon avec Pôles"
        },
        {
            "lons": np.array([0, 45]), "lats": np.array([0, 10]), "grtyp": "G", "ni": 180, "ig1": 1, "ig2": 1, "filename": "Grid_G.fst", "label": "Gaussienne"
        },
        # quarantaine pour grille Q car pb interpolation
        #{
         #   "lons": np.array([0, 45]), "lats": np.array([0, 10]), "grtyp": "Q", "ni": 180, "num_elem": 36, "num_solpts": 5, "ig1": 0x420000, "ig2": 0xa4fa00, "ig3": 0x660000, "filename": "Grid_Q.fst", "label": "Cubed Sphere"
        #}
    ]

    # Création des fichiers sources
    for config in grids_config:
        args = {k: v for k, v in config.items() if k not in ['label', 'filename']}
        generate_grid(filename=config['filename'], **args)

    # Création des fichiers interpolés
    for src in grids_config:
        for dest in grids_config:
            output_file = f"interp_{src['grtyp']}_to_{dest['grtyp']}.fst"
            
            print(f"Calcul : {src['grtyp']} ({src['ni']} pts) -> {dest['grtyp']} ({dest['ni']} pts)")
            
            try:
                interpolate_grid(src['filename'], output_file, dest)
                
            except Exception as e:
                print(f"  [ERREUR] Échec de {src['grtyp']} vers {dest['grtyp']} : {e}")
