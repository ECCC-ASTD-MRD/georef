L'objectif de ce stage est d'améliorer les programmathèques de
prévision numérique du temps du Centre de prévision météorologique et
environnementale du Canada (CPMEC). Plus spécifiquement, le stage
concerne les grilles utilisées pour encoder des données
météorologiques comme la température de l'air, la pression
atmosphérique, ou la salinité de l'eau.

Les grilles de données utilisées dans les modèles originaux de
prévision numérique du temps (PNT), depuis les années 1950, étaient
cartésiennes, voir la Figure 1 pour une illustration. Les grilles
cartésiennes ont une structure naturellement adaptée à la structure
matérielle des ordinateurs utilisés en PNT. En d'autre mots,
l'enregistrement de données sur des grilles cartésiennes est peu
coûteuses en ressources informatiques comme la mémoire. Cela est
important en PNT opérationnelle car ça aide à diffuser rapidement les
résultats aux usagers des prévisions, par exemple les météorologistes
aux bureaux régionaux d'Environnement et Changement climatique Canada.
Les prévisions perdent en valeur plus le temps passe.

Les avancées en PNT depuis les années 1950 ont mené à des modèles plus
complexes, pour lesquels les grilles cartésiennes ne sont pas
optimales. Les figures 2 à 4 donnent des exemples. Les
programmathèques de PNT du CPMEC ont besoin de mises-à-jour à leur
fonctions afin de supporter les nouvelles grilles. Ces grilles sont
déjà en usage dans des modèles opérationnels, mais sont traitées avec
des utilitaires spéciaux, par exemple cstintrp pour les grilles ORCA.
Ce stage vise à rapatrier les fonctionnalités éparpillées dans des
utilitaires spéciaux aux principales programmathèques du CPMEC,
notamment [rmnlib](https://github.com/armnlib/librmn) et libeerUtils.
Il est prévu que cette standardisation réduira l'effort pour
l'ensemble du CPMEC dans le traitement de ses modèles de PNT.

Figure 1 : Grille globale cartésienne. La résolution de la grille est de , une résolution grossière pour des fins d'illustration.

![Grille océanique ORCA](https://hpfx.collab.science.gc.ca/~map007/o_grid.png)
Figure 2 : Grille océanique ORCA avec résolution moyenne de 10°

Figure 3 : Grille icosaédrique

Figure 4 : Grille MESH d'éléments finis

# Exemples

## libgeoref

```shell
. r.load.dot main/opt/cmake/cmake-3.16.4
export SSM_DEV=$(mktemp --directory)
./makeit.sh
```

## Python-RPN

```shell
ssh eccc-ppp4.science.gc.ca
. r.load.dot comm/eccc/all/opt/intelcomp/intelpsxe-cluster-19.0.3.199 \
  rpn/libs/19.5 \
  rpn/MIG/ENV/rpnpy/2.1.0
cd test
./demo_a_grid.py
# Test libgeoref interpolation
./demo_a_grid.py --georef
./demo_defGrid_L.py
./demo_o_grid.py
. r.load.dot rpn/OCEAN/cstint-3.2.8 \
  cmd/cmds/apps/SPI/beta
./NEMOInterp_sinus.tcl out/O_sinus.fst out/l_grid.fst
./comparaison_interp.py
```
Cette séquence devrait produire les fichiers `a_data.png`,
`a_grid.fst`, `a_grid.png`, `l_grid.fst`, `l_grid.png`, `o_data.png`, `o_grid.png`, `O_sinus.fst`, `out.csintrp`, `out.csintrp.avg` et `out.spi` dans le dossier `out`.

La dernière commande utilise les fichiers `out.*` produits par le script `NEMOInterp_sinus.tcl`.

Le dossier `GRIDS` contient le fichier `O.fstd` qui donne les positions des points de grille de la grille ORCA.

# Références

- [Grid Types Supported by RPN Standard Files](https://science:science@collaboration.cmc.ec.gc.ca/science/si/eng/si/misc/grilles.html)
- [Python-RPN/grids](https://wiki.cmc.ec.gc.ca/wiki/Python-RPN/grids)
- [fst_missing](http://armnlib.uqam.ca/armnlib/Docs/fst_missing.html)
