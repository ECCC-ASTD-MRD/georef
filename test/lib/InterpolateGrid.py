from georef import cubed_sphere
import georef
import numpy as np
from rmn import fst24_file
from GenerateGrid import generate_grid

def validate_interpolation(src_file, target_file):
    """
    Interpole src_file vers la grille de dest_file 
    et calcule l'erreur sans créer de fichier intermédiaire.
    """
    
    with fst24_file(src_file, "R") as f_src, fst24_file(target_file, "R") as f_target:
        # Lecture du fichier source 
        src_record = next(iter(f_src.new_query(nomvar="GRID")), None)
        q_src = next(iter(f_src.new_query(nomvar="DIST")), None)
        if not src_record or not q_src: return None
        src_grid = georef.GeoRef.fromrecord(src_record)
        data_src = q_src.data.astype(np.float32)
        src_type = src_record.grtyp.strip()

        # Lecture du fichier cible
        grid_record = next(iter(f_target.new_query(nomvar="GRID")), None)
        rec_ref_data = next(iter(f_target.new_query(nomvar="DIST")), None)
        if not grid_record or not rec_ref_data: return None
        target_grid = georef.GeoRef.fromrecord(grid_record)
        data_ref = rec_ref_data.data.astype(np.float32)
        tgt_type = grid_record.grtyp.strip()

        # Ajout des options
        options = georef.GeoOptions(PolarCorrect=False)
        if src_type == "Q" or tgt_type == "Q":
            options.Interp = 2

        # Interpolation 
        data_interp = target_grid.interp(src_grid, data_src, options=options)

        # Calcul de la différence
        diff = (data_interp.ravel() - data_ref.ravel()) / np.max(np.abs(data_ref))

        # Calcul des normes
        norm = np.linalg.norm(diff)/diff.size 
        max_err = np.linalg.norm(diff, ord=np.inf)

        # Définition des seuils
        THRESHOLD_NORM = 2.8e-4
        THRESHOLD_MAX_ERR = 4.2e-1

        if norm > THRESHOLD_NORM or max_err > THRESHOLD_MAX_ERR or not np.all(np.isfinite(data_interp)):
            error_msg = (
                f"\n! ERREUR : Seuils dépassés pour {src_file} -> {target_file}\n"
                f"Valeurs : Norme={norm:.6e} (max {THRESHOLD_NORM}), "
                f"MaxErr={max_err:.6e} (max {THRESHOLD_MAX_ERR})"
            )
            raise ValueError(error_msg)

        # Affichage des résultats
        print(f"[*] Succès {src_file} -> {target_file}")
        print(f"    Norme relative : {norm:.6e}")
        print(f"    Erreur Max     : {max_err:.6e}")
        print("-" * 40)

        return data_interp

    
base_config = {
    "lons": np.array([0, 45]),
    "lats": np.array([0, 10]),
    "ni": 180,
    "nj": 90
}

if __name__ == "__main__":
    #np.set_printoptions(precision=6,edgeitems=7,linewidth=256)

    # Dictionnaire des fichiers sources à créer
    grids_config = {
        # Grilles globales 
        # TODO: Pour faire la grille U, il faut que la grille Z fonctionne (car U = 2 grilles Z concaténées)
        # TODO: Grille V: erreur de segmentation fault, non défini dans DefRPNXG
        "A": {**base_config, "grtyp": "A", "filename": "Grid_A.fst", "label": "Lat-Lon Equidistante"},
        "A2": {**base_config, "ni": 170, "nj": 80, "grtyp": "A", "filename": "Grid_A2.fst", "label": "Lat-Lon Equidistante"},
        "B": {**base_config, "grtyp": "B", "filename": "Grid_B.fst", "label": "Lat-Lon avec Pôles"},
        "G": {**base_config, "grtyp": "G", "filename": "Grid_G.fst", "label": "Gaussien"},
        "Q": {**base_config, "grtyp": "Q", "ig1": 0x800000, "ig2": 0x800000, "ig3": 0x800000, "ig4": cubed_sphere.encodeig4(18,3), "filename": "Grid_Q.fst", "label": "Cubed Sphere"},
        "Q2": {**base_config, "grtyp": "Q", "ig1": 0x800000, "ig2": 0x800000, "ig3": 0x800000, "ig4": cubed_sphere.encodeig4(16,3), "filename": "Grid_Q2.fst", "label": "Cubed Sphere"},
 
        # # Grilles Hemisphères (divisé par 2)
        # # TODO: interpolation N->N et S->S (présence de NaNs)
        "N": {**base_config, "grtyp": "N", "filename": "Grid_N.fst", "ig1":1, "ig2":1, "ig3":1, "ig4":1, "label": "Hemisphere Nord"},
        "S": {**base_config, "grtyp": "S", "filename": "Grid_S.fst", "ig1":1, "ig2":1, "ig3":1, "ig4":1, "label": "Hemisphere Sud"}, 
    }

    interpolation_tasks = [
        # Grille A
        ("A", "A"),
        ("A", "A2"),
        ("A", "B"),
        ("A", "G"),
        ("A", "Q"),
        ("A", "N"),
        ("A", "S"),
        ("A", "Q2"),
        # Grille A2
        ("A2", "A"),
        ("A2", "A2"),
        ("A2", "B"),
        ("A2", "G"),
        ("A2", "Q"),
        ("A2", "N"),
        ("A2", "S"),
        ("A2", "Q2"),
        # Grille B
        ("B", "A"),
        ("B", "A2"),
        ("B", "B"),
        ("B", "G"),
        ("B", "Q"),
        ("B", "N"),
        ("B", "S"),
        ("B", "Q2"),
        # Grille G
        ("G", "A"),
        ("G", "A2"),
        ("G", "B"),
        ("G", "G"),
        ("G", "Q"),
        ("G", "N"),
        ("G", "S"),
        ("G", "Q2"),
        # Grille Q
        ("Q", "A"),
        ("Q", "A2"),
        ("Q", "B"),
        ("Q", "G"),
        ("Q", "Q"),
        ("Q", "N"),
        ("Q", "S"),
        ("Q", "Q2"),
        # Grille Q2
        ("Q2", "A"),
        ("Q2", "A2"),
        ("Q2", "B"),
        ("Q2", "G"),
        ("Q2", "Q"),
        ("Q2", "N"),
        ("Q2", "S"),
        ("Q2", "Q2"),
    ]



    # Création des fichiers sources
    for config in grids_config.values():
        args = {k: v for k, v in config.items() if k not in ['label', 'filename']}
        generate_grid(filename=config['filename'], **args)
    
    for src_key, tgt_key in interpolation_tasks:
        src_file = grids_config[src_key]['filename']
        tgt_file = grids_config[tgt_key]['filename']
        validate_interpolation(src_file,tgt_file)