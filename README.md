# Pour démarrer

Conda requis car Cartopy pas installé sur gpscc2.collab.science.gc.ca,
référence : https://gitlab.science.gc.ca/hpc/support/issues/5
```shell
export PATH=~kro001/miniconda3/bin:"$PATH"
. ssmuse-sh -x comm/eccc/all/opt/intelcomp/intelpsxe-cluster-19.0.3.199 \
-x rpn/libs/19.5 \
-x rpn/MIG/ENV/rpnpy/2.1.0
./demo_defGrid_L.py
```
Cette séquence devrait produire les fichiers `l_grid.fst` et `l_grid.png`.

# Références

- [Grid Types Supported by RPN Standard Files](https://collaboration.cmc.ec.gc.ca/science/si/eng/si/misc/grilles.html)
- [Python-RPN/grids](https://wiki.cmc.ec.gc.ca/wiki/Python-RPN/grids)
- [fst_missing](http://armnlib.uqam.ca/armnlib/Docs/fst_missing.html)
