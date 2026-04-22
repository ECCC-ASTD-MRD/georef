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
    # Lecture du fichier source
    with fst24_file(src_file, "R") as f_src:
        grid_src = next(iter(f_src.new_query(nomvar="GRID")), None)
        q_src = next(iter(f_src.new_query(nomvar="DIST")), None)
        if not grid_src or not q_src: return None

        src_grid = georef.GeoRef.fromrecord(grid_src)
        src_grid.shape = (grid_src.ni, grid_src.nj, grid_src.nk)
        data_src = q_src.data.astype(np.float32)
        src_type = grid_src.grtyp.strip()

        # lats_src, lons_src = src_grid.getll()
        
        # print(f"\n--- Coordonnées SOURCE ({src_type}) ---")
        # print(f"Shape des lats : {lats_src.shape}")
        # print(f"Latitude Min   : {np.min(lats_src):.4f}")
        # print(f"Latitude Max   : {np.max(lats_src):.4f}")
        # print(lats_src)
        # print(f"Shape des lons : {lons_src.shape}")
        # print(f"Longitude Min   : {np.min(lons_src):.4f}")
        # print(f"Longitude Max   : {np.max(lons_src):.4f}")
        # print(lons_src)

    # Lecture du fichier cible
    with fst24_file(target_file, "R") as f_target:
        grid_record = next(iter(f_target.new_query(nomvar="GRID")), None)
        rec_ref_data = next(iter(f_target.new_query(nomvar="DIST")), None)
        if not grid_record or not rec_ref_data: return None

        target_grid = georef.GeoRef.fromrecord(grid_record)
        target_grid.shape = (grid_record.ni, grid_record.nj, grid_record.nk)
        data_ref = rec_ref_data.data.astype(np.float32)
        tgt_type = grid_record.grtyp.strip()

        # lats_tgt, lons_tgt = target_grid.getll()
        
        # print(f"\n--- Coordonnées TARGET ({tgt_type}) ---")
        # print(f"Shape des lats : {lats_tgt.shape}")
        # print(f"Latitude Min   : {np.min(lats_tgt):.4f}")
        # print(f"Latitude Max   : {np.max(lats_tgt):.4f}")
        # print(lats_tgt)
        # print(f"Shape des lons : {lons_tgt.shape}")
        # print(f"Longitude Min   : {np.min(lons_tgt):.4f}")
        # print(f"Longitude Max   : {np.max(lons_tgt):.4f}")
        # print(lons_tgt)

    # Grilles exclues à l'interpolation
    excluded_sources = {'N', 'S'}
    is_forbidden = (src_type in excluded_sources) 
    should_interpolate = not is_forbidden

    if not should_interpolate:
        print(f"[-] Skip : Interpolation {src_type} -> {tgt_type} non supportée/voulue.")
        return np.full(target_grid.shape, np.nan, dtype=np.float32)

    # Exécution de l'interpolation
    print(f"[+] Interpolation en cours : {src_type} -> {tgt_type}")
    options = georef.GeoOptions(Interp=2) if src_type or tgt_type == "Q" else None
    data_interp = target_grid.interp(src_grid, data_src, options=options)

    # Calcul de la différence
    diff = (data_interp - data_ref).ravel() / np.max(np.abs(data_ref))
    # print(f"data_interp:{data_interp}")
    # print(f"data_ref:{data_ref}")
    # print(f"diff:{diff}")
    
    # Calcul des normes
    norm = np.linalg.norm(diff)/diff.size 
    max_err = np.linalg.norm(diff, ord=np.inf)

    # Définition des seuils
    THRESHOLD_NORM = 1.0e-3
    THRESHOLD_MAX_ERR = 8e-1

    if norm > THRESHOLD_NORM or max_err > THRESHOLD_MAX_ERR:
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
    # Dictionnaire des fichiers sources à créer
    grids_config = [
        # Grilles globales 
        {
            **base_config, "grtyp": "A", "filename": "Grid_A.fst", "label": "Lat-Lon Equidistante"
        },
        # {
        #     **base_config, "ni": 170, "nj": 80, "grtyp": "A", "filename": "Grid_A2.fst", "label": "Lat-Lon Equidistante"
        # },
        # {
        #     **base_config, "grtyp": "B", "filename": "Grid_B.fst", "label": "Lat-Lon avec Pôles"
        # },
        # {
        #     **base_config, "grtyp": "G", "filename": "Grid_G.fst", "label": "Gaussien"
        # },

        # # Grilles Hemisphères (divisé par 2)
        # # TODO: interpolation N->N et S->S (présence de NaNs)
        # {
        #     **base_config, "grtyp": "N", "filename": "Grid_N.fst", "ig1":1, "ig2":1, "ig3":1, "ig4":1, "label": "Hemisphere Nord"
        # },
        # {
        #     **base_config, "grtyp": "S", "filename": "Grid_S.fst", "ig1":1, "ig2":1, "ig3":1, "ig4":1, "label": "Hemisphere Sud"
        # },
        {
           **base_config, "grtyp": "Q", "ig1": 0x800000, "ig2": 0x800000, "ig3": 0x800000, "ig4": cubed_sphere.encodeig4(18,3), "filename": "Grid_Q.fst", "label": "Cubed Sphere"
        },
        # {
        #    **base_config, "grtyp": "Q", "ig1": 0x700000, "ig2": 0x800000, "ig3": 0x800000, "ig4": cubed_sphere.encodeig4(16,3), "filename": "Grid_Q2.fst", "label": "Cubed Sphere"
        # },
        # TODO: Pour faire la grille U, il faut que la grille Z fonctionne (car U = 2 grilles Z concaténées)
        # TODO: Grille V: erreur de segmentation fault, non défini dans DefRPNXG 
    ]

    # Création des fichiers sources
    for config in grids_config:
        args = {k: v for k, v in config.items() if k not in ['label', 'filename']}
        generate_grid(filename=config['filename'], **args)

    # Création des fichiers interpolés + Comparaison fichier source et interpolés
    for src in grids_config:
        for dest in grids_config:
            validate_interpolation(src['filename'], dest['filename'])
                
