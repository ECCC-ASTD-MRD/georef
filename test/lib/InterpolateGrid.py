import georef
import numpy as np
from rmn import fst24_file
from GenerateGrid import generate_grid

def validate_interpolation(src_file, dest_file):
    """
    Interpole src_file vers la grille de dest_file 
    et calcule l'erreur sans créer de fichier intermédiaire.
    """
    # Lecture du fichier source
    with fst24_file(src_file, "R") as f_src:
        q_src = f_src.new_query(nomvar="DIST")
        rec_src = next(iter(q_src), None)
        if rec_src is None: return None
        src_geo = georef.GeoRef(rec_src.ni, rec_src.nj, rec_src.grtyp, rec_src.ip1, rec_src.ip2, rec_src.ip3, rec_src.ig4, f_src)
        data_src = rec_src.data.astype(np.float32)

    # Lecture du fichier destination
    with fst24_file(dest_file, "R") as f_dest:
        # On cherche la grille destination
        q = f_dest.new_query(nomvar="GRID")
        grid_record = next(iter(q), None)

        # On cherche les données pour comparer
        q_dest_data = f_dest.new_query(nomvar="DIST")
        rec_dest_data = next(iter(q_dest_data), None)

        dest_grid = georef.GeoRef.fromrecord(grid_record)
        dest_grid.shape = (grid_record.ni, grid_record.nj, grid_record.nk)
        data_dest = rec_dest_data.data.astype(np.float32)

        data_interp = dest_grid.interp(src_geo, data_src)

    # Calcul de la différence
    diff = (data_interp - data_dest).ravel()
    
    # Calcul des normes
    norm = np.linalg.norm(diff)
    max_err = np.linalg.norm(diff, ord=np.inf)

    print(f"Interpolation {src['grtyp']} vers {dest['grtyp']}")
    print(f"Norme par defaut : {norm:.6e}")
    print(f"Erreur Max       : {max_err:.6e}")
    print("-" * 30)


if __name__ == "__main__":
    # Dictionnaire des fichiers sources à créer
    grids_config = [
        {
            "lons": np.array([0, 45]), "lats": np.array([0, 10]), "grtyp": "A", "ni": 180, "nj": 90, "filename": "Grid_A.fst", "label": "Lat-Lon Equidistante"
        },
        {
            "lons": np.array([0, 45]), "lats": np.array([0, 10]), "grtyp": "B", "ni": 180, "nj": 90, "filename": "Grid_B.fst", "label": "Lat-Lon avec Pôles"
        },
        # TODO: Grille G a des données NaaN dans les fichiers .fst
        # TODO: Grille Q, la fonction cubic interpolation n’est pas implémentée pour la grille Q
    ]

    # Création des fichiers sources
    for config in grids_config:
        args = {k: v for k, v in config.items() if k not in ['label', 'filename']}
        generate_grid(filename=config['filename'], **args)

    # Création des fichiers interpolés + Comparaison fichier source et interpolés
    for src in grids_config:
        for dest in grids_config:
                        
            try:
                validate_interpolation(src['filename'], dest['filename'])
                
            except Exception as e:
                print(f"  [ERREUR] Échec de {src['grtyp']} vers {dest['grtyp']} : {e}")
