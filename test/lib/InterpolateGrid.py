import georef
import numpy as np
from rmn import fst24_file
from GenerateGrid import generate_grid
import sys

def validate_interpolation(src_file, target_file):
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

    # Lecture du fichier cible
    with fst24_file(target_file, "R") as f_target:
        # On cherche la grille cible
        q = f_target.new_query(nomvar="GRID")
        grid_record = next(iter(q), None)

        # On cherche les données pour comparer
        q_ref_data = f_target.new_query(nomvar="DIST")
        rec_ref_data = next(iter(q_ref_data), None)

        target_grid = georef.GeoRef.fromrecord(grid_record)
        target_grid.shape = (grid_record.ni, grid_record.nj, grid_record.nk)
        data_ref = rec_ref_data.data.astype(np.float32)

        options = None
        if grid_record.grtyp == "Q":
            options = georef.GeoOptions(Interp=3)

        data_interp = target_grid.interp(src_geo, data_src, options=options)

    # Calcul de la différence
    diff = (data_interp - data_ref).ravel() / np.max(np.abs(data_ref))
    
    # Calcul des normes
    norm = np.linalg.norm(diff)/diff.size 
    max_err = np.linalg.norm(diff, ord=np.inf)

    # Définition des seuils
    THRESHOLD_NORM = 6.0e-4
    THRESHOLD_MAX_ERR = 4e-2

    if norm > THRESHOLD_NORM or max_err > THRESHOLD_MAX_ERR:
        # Message d'erreur détaillé avant de stopper
        error_msg = (
            f"ERREUR : Seuils de tolérance dépassés pour {src['grtyp']} -> {dest['grtyp']} !\n"
            f"Seuils : Norme < {THRESHOLD_NORM}, MaxErr < {THRESHOLD_MAX_ERR}\n"
            f"Valeurs actuelles : Norme = {norm:.6e}, MaxErr = {max_err:.6e}"
        )
        raise ValueError(error_msg)

    print(f"Interpolation {src['grtyp']} vers {dest['grtyp']}")
    print(f"Norme par defaut : {norm:.6e}")
    print(f"Erreur Max       : {max_err:.6e}")
    print("-" * 30)

    
base_config = {
    "lons": np.array([0, 45]),
    "lats": np.array([0, 10]),
    "ni": 180,
    "nj": 90
}

if __name__ == "__main__":
    # Dictionnaire des fichiers sources à créer
    grids_config = [
        {
            **base_config, "grtyp": "A", "filename": "Grid_A.fst", "label": "Lat-Lon Equidistante"
        },
        {
            **base_config, "grtyp": "B", "filename": "Grid_B.fst", "label": "Lat-Lon avec Pôles"
        },
        {
            **base_config, "grtyp": "G", "filename": "Grid_G.fst", "label": "Gaussien"
        },

        # Problème avec l'interpolation N->N et S->S
        {
            **base_config, "grtyp": "N", "filename": "Grid_N.fst", "ig1":1, "ig2":1, "ig3":1, "ig4":1, "label": "Hemisphere Nord"
        },
        {
            **base_config, "grtyp": "S", "filename": "Grid_S.fst", "ig1":1, "ig2":1, "ig3":1, "ig4":1, "label": "Hemisphere Sud"
        },
        #{
         #   **base_config, "grtyp": "Q", "ig4": 1801, "filename": "Grid_Q.fst", "label": "Cubed Sphere"
        #},
        # TODO: Pour faire la grille U, il faut que la grille Z fonctionne (car U = 2 grilles Z concaténées)
    ]

    # Création des fichiers sources
    for config in grids_config:
        args = {k: v for k, v in config.items() if k not in ['label', 'filename']}
        generate_grid(filename=config['filename'], **args)

    # Création des fichiers interpolés + Comparaison fichier source et interpolés
    for src in grids_config:
        for dest in grids_config:
            validate_interpolation(src['filename'], dest['filename'])
                
