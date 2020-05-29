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
- remplacer appel à `rmn.ezsint` par `libgeoref.c_ezsint`
  - [diagnostic d'erreur de segmentation](https://github.com/jeixav/stage_2020/issues/6)

# Semaine 5

- remplacer appel à `rmn.ezsint` par `libgeoref.c_ezsint`
- adapter `c_ezsint` pour grille ORCA
  - identifier où `c_ezsint` échoue avec grille ORCA

# Semaines 6 & 7

- adapter `c_ezsint` pour grille ORCA (suite)
  - tenir compte [des fonctions libeerutils GeoRef](https://gitlab.science.gc.ca/ECCC_CMOE_MODELS/libeerutils/blob/e279df3a8176a9c63e9f7efdae79935288715bd3/src/GeoRef.c) et [EZGrid](https://gitlab.science.gc.ca/ECCC_CMOE_MODELS/libeerutils/blob/e279df3a8176a9c63e9f7efdae79935288715bd3/src/EZGrid.c)
  - ajout du [code quadtree de libeerUtils](https://gitlab.science.gc.ca/ECCC_CMOE_MODELS/libeerutils/blob/e279df3a8176a9c63e9f7efdae79935288715bd3/src/QTree.c) dans ezscint (pointeur dans struct grille) pour index spatial
  - ajout du code d'interpolation d'un quadrilatère arbitraire [Vertex_Map](https://gitlab.science.gc.ca/ECCC_CMOE_MODELS/libeerutils/blob/e279df3a8176a9c63e9f7efdae79935288715bd3/src/Vertex.c#L36-98)
  - se baser sur [`GeoRef_BuildIndex`](https://gitlab.science.gc.ca/ECCC_CMOE_MODELS/libeerutils/blob/e279df3a8176a9c63e9f7efdae79935288715bd3/src/GeoRef.c#L1098-1204), [`GeoRef_RPNProject`](https://gitlab.science.gc.ca/ECCC_CMOE_MODELS/libeerutils/blob/e279df3a8176a9c63e9f7efdae79935288715bd3/src/GeoRef_RPN.c#L340-413) et [`GeoRef_RPNUnProject`](https://gitlab.science.gc.ca/ECCC_CMOE_MODELS/libeerutils/blob/e279df3a8176a9c63e9f7efdae79935288715bd3/src/GeoRef_RPN.c#L415-597) pour le flow

# Semaines 8 à 10

- interpoler de grille M à L avec les API de SPI
- grilles M (triangle mesh), 
  - ajout du code de gestion des triangles et interpolation barycentrique ([libeerUtils `Triangle.c`](https://gitlab.science.gc.ca/ECCC_CMOE_MODELS/libeerutils/blob/e279df3a8176a9c63e9f7efdae79935288715bd3/src/Triangle.c))
  - se baser sur `GeoRef_BuildIndex`, `GeoRef_RPNProject` et `GeoRef_RPNUnProject` pour le flow
- grilles W (WKT) et Lambert
  - utilisera [proj4](https://proj.org/) et [gdal](https://gdal.org/), à préciser

# Semaines 10 à 14
  - reengineering plus profond d'ezscint vers un « vrai » libgeoref

# Semaine 15 & 16

- présentation des résultats
- rédaction de rapport de stage
