# Semaine 1

- installation de MobaXterm, Visual Studio Code & Anaconda pour travail à distance
- accès par tunnel SSH à https://wiki.cmc.ec.gc.ca/ & https://gitlab.science.gc.ca/
- exécution d'exemples [`demo_a_grid.py`](https://github.com/jeixav/stage_2020/blob/5c2c86459d920a2866b46d8af58fd886be200ac3/test/demo_a_grid.py) et [Plot GIOPS Forecast Data with Basemap](https://wiki.cmc.ec.gc.ca/wiki/Talk:Python-RPN/2.1/examples#Plot_GIOPS_Forecast_Data_with_Basemap)
- recherche de documentation Numpy et Scipy

# Semaine 2

- production de graphiques Python sur grille ORCA : [graphique de grille](https://hpfx.collab.science.gc.ca/~map007/o_grid.png) et [graphique de données analytiques sinus](https://hpfx.collab.science.gc.ca/~map007/o_data.png)
- interpoler de grille ORCA à L avec [cstintrp](https://wiki.cmc.ec.gc.ca/wiki/Cstintrp_V3)
- interpoler de grille ORCA à [RPN L](https://science:science@collaboration.cmc.ec.gc.ca/science/si/eng/si/misc/grilles.html#LatLon) avec [API de SPI](https://wiki.cmc.ec.gc.ca/wiki/SPI/Documentation#Developer_documentation)

# Semaine 3

- interpoler de grille ORCA à [RPN L](https://science:science@collaboration.cmc.ec.gc.ca/science/si/eng/si/misc/grilles.html#LatLon) avec [API de SPI](https://wiki.cmc.ec.gc.ca/wiki/SPI/Documentation#Developer_documentation) (suite)
- [configurer ThinLinc pour travail à distance](https://1drv.ms/w/s!AmH_Shsw9Hrnvyo9b08sRvWJyE7v)
- remplacer [appel à `rmn.ezsint`](https://github.com/jeixav/stage_2020/blob/5c2c86459d920a2866b46d8af58fd886be200ac3/test/demo_a_grid.py#L71) par [`libgeoref.c_ezsint`](https://github.com/jeixav/stage_2020/blob/5c2c86459d920a2866b46d8af58fd886be200ac3/src/ezsint.c#L33-L145)
  - utilise [Python ctypes](https://docs.python.org/3/library/ctypes.html)
  - référence: [librmn avec python-RPN](https://github.com/meteokid/python-rpn/tree/master/lib/rpnpy/librmn)

# Semaine 4

- [correction d'une erreur dans `NEMOInterp_sinus.tcl`](https://github.com/jeixav/stage_2020/pull/8)
- vérifier égalité des résultats de cstintrp et API SPI
  - https://github.com/jeixav/stage_2020/issues/5
  - [refactoring du code de vérification des résultats](https://github.com/jeixav/stage_2020/pull/9)
- tester programme C `interpolate.c`
- [relecture de la demande de fusion de la branche feat/move_to_ppp3_4](https://github.com/jeixav/stage_2020/pull/10)
- remplacer appel à `rmn.ezsint` par `libgeoref.c_ezsint` (suite)
  - [diagnostic d'erreur de segmentation](https://github.com/jeixav/stage_2020/issues/6)

# Semaine 5

- remplacer appel à `rmn.ezsint` par `libgeoref.c_ezsint` (suite)
- générer automatiquement de la documentation avec Doxygen
- adapter `c_ezsint` pour grille ORCA
  - identifier où `c_ezsint` échoue avec grille ORCA
  - utilisation du débogueur [TotalView](https://portal.science.gc.ca/confluence/x/14Lr)

# Semaines 6 & 7

- [fusionner structure EZ _Grille à SPI TGeoRef](https://github.com/jeixav/stage_2020/commit/afa40547c2e983ee8e4dfd40ab2f52a6ee60d6ad)

# Semaine 8

- changer le paramètre `gridid` à pointeur à `TGeoRef` qui inclus les paramètres de grilles
- retire code de hachage et tableaux `Grilles[][]` maintenant inutiles

# Semaines 9 à 14

- retirer fonctions `c_addgrid`, `c_ezqkdef` et `c_ezgdef` puisque devenues inutiles
- utiliser `PTR_AS_INT` avec fonctions Fortran; adresses des structures `TGeoRef` dans des `int` que l'on re-type (recast) en `TGeoRef` pour les appels vers le C
- optimiser et nettoyer les fonctions

# Semaine 15 & 16

- présentation des résultats
- rédaction de rapport de stage
