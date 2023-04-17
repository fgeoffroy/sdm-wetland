# sdm-wetland

This repository contains R scripts for running Species Distribution Models.

- Study area: 3 departments in France (Hérault, Gard and Bouches-du-Rhône)
- Species: *Limonium girardianum* and *Althenia filiformis*. Occurence data come from [GBIF](https://www.gbif.org/fr/)
- Environmental data
  - Bioclimatic variables from [WorldClim](https://www.worldclim.org/data/index.html)
  - Elevation from [OpenTopography](https://portal.opentopography.org/raster?opentopoID=OTSDEM.092022.3035.1)
  - Land cover from [SEC OSO](https://www.theia-land.fr/en/product/land-cover-map/)
  - Wetland probability [UMR 1069 SAS INRAE](https://agroenvgeo.data.inra.fr/geonetwork/srv/api/records/518b3e0a-ee55-40cb-a3ed-da00e60505aa)
  - Pedological variables from [SoilGrids](https://www.isric.org/explore/soilgrids)

In `data_preparation.R`, spatial data are pre-processed and standardized. In `pre-analysis.R`, correlation and PCA analysis are performed. And in `model.R`, SDM models are run using the *biomod2* package.