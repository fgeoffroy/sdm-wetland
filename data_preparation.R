library(raster)
library(rgbif)
library(biomod2)
library(ggplot2)
library(pdfcrop)
# library(RColorBrewer)
library(gridExtra)
# library(geodata)
library(ade4)
library(ncdf4)
library(terra)
library(corrplot)
library(factoextra)

library(XML)
library(rgdal)
library(gdalUtils)
gdalUtils::gdal_setInstallation("C:/Program Files/QGIS 3.30.0/bin/")


setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Define study area (étang de Berre ajouté sur QGIS)
departments <- vect("../data/france_vectors/france-departments_without_holes.gpkg")
plot(departments)
departments_names = c("HÃ©rault", "Gard", "Bouches-du-RhÃ´ne")
study_area <- subset(departments, departments$name %in% departments_names)
plot(study_area)

# pdf("../results/study_area.pdf")
plot(departments, axes=F)
plot(study_area, col='#ECB176', add=T)
text(3.4, 43.6, '34', add=T, cex=0.8)
text(4.2, 44.05, '30', add=T, cex=0.8)
text(5, 43.6, '13', add=T, cex=0.8)
# dev.off()

prefecture_vect <- vect(data.frame(name=c("Montpellier", "Nimes", "Marseille"),
                                   y=c(43.62505, 43.833328, 43.296482),
                                   x=c(3.862038, 4.35, 5.36978)), geom=c("x", "y"))
prefecture_text_vect <- vect(data.frame(name=c("Montpellier", "Nimes", "Marseille"),
                                   y=c(43.68, 43.88, 43.35),
                                   x=c(3.862038, 4.35, 5.36978)), geom=c("x", "y"))



# get species taxon keys
# "Limonium girardianum"    8076932
# "Pancratium maritimum"    2853283
# "Althenia filiformis"     6386027
# "Ixobrychus minutus"
# "Netta rufina"

# spp_Girard <- name_suggest(q = "Althenia filiformis",
#                             rank = "species",
#                             limit = 10000)
# spp_Girard$data$canonicalName
# spp_Girard$data[grepl("Althenia filiformis", spp_Girard$data$canonicalName), ]

# get species occurrence data
# gbif <- occ_search(taxonKey = 6386027, country='FR',
#                    fields = c('name', 'key', 'decimalLatitude','decimalLongitude', 'coordinateUncertaintyInMeters'),
#                    # fields='all',
#                    hasCoordinate=T,
#                    limit=8000
#                    # limit=70
#                    )
# gbif$meta$count
# data <- gbif$data
# nrow(data)
# names(data)[names(data) == 'decimalLongitude'] <- 'x'
# names(data)[names(data) == 'decimalLatitude'] <- 'y'
# write.csv(data, "../data/species_occurence/GBIF/Limonium_girardianum_occurences.csv")
# write.csv(data, "../data/species_occurence/GBIF/Pancratium_maritimum_occurences.csv")
# write.csv(data, "../data/species_occurence/GBIF/Althenia_filiformis_occurences.csv")




data <- read.csv("../data/species_occurence/GBIF/Limonium_girardianum_occurences.csv")
data <- read.csv("../data/species_occurence/GBIF/Althenia_filiformis_occurences.csv")


nrow(data)
data <- data[data$coordinateUncertaintyInMeters <= 5000, ]
nrow(data)
# data <- data[data$coordinateUncertaintyInMeters <= 1000, ]
# data <- data[data$coordinateUncertaintyInMeters >= 4000 &
#                data$coordinateUncertaintyInMeters <= 5000, ]
hist(data$coordinateUncertaintyInMeters, 50)

data_vect <- vect(data, geom=c("x", "y"))
data_vect <- crop(data_vect, study_area)
data_vect <- mask(data_vect, study_area)
occurences <- as.data.frame(data_vect, geom='XY')

length(data_vect)
plot(data_vect)
plot(study_area, add=T)





# import WorldClim data
# current_clim_data <- getData('worldclim', var='bio', res=10, path='data')
# current_clim_data <- worldclim_country('FRA', var='bio', path='data')




humid_raster <- rast('../data/humidity/MPH_FR_WGS84.tif')
# humid_raster <- project(x= humid_raster, y =  "EPSG:4326", method = "bilinear")
humid_raster <- crop(humid_raster, study_area)
humid_raster <- mask(humid_raster, study_area)
humid_raster[humid_raster$MPH_FR_WGS84 == 51] <- 4
names(humid_raster) <- "wetness_index"
print(humid_raster)
cellSize(humid_raster)

# colors = c("white", rev(colorspace::sequential_hcl(4)))
colors = rev(hcl.colors(5, "Blues"))
colors[1] = "white"

# pdf("../results/zones_humides.pdf")
pdf("../results/saladelle_occurences.pdf")
# pdf("../results/althenie_occurences.pdf")
plot(humid_raster, legend=F, col=colors, axes=F)
sbar(50, xy='topleft', type='bar', below='km')
legend("topright", c("non-humide", 'proba assez forte', 'proba forte',
                     'proba très forte', "plan d'eau"), fill=colors, inset=c(0.01, 0.05))

points(data_vect, pch=20, col="red", cex=1, alpha=0.5)
# points(prefecture_vect, pch=20, col="red", cex=2)
# text(prefecture_text_vect, pch=20, col="red", cex=1, labels=prefecture_vect$name)
plot(study_area, add=T)
dev.off()

# df_humid <- na.omit(as.data.frame(humid_raster))

data_vect$humidity <- extract(humid_raster, data_vect)$wetness_index
barplot(prop.table(table(as.factor(data_vect$humidity))))
length(data_vect[data_vect$humidity == 51])
length(data_vect)
length(data_vect[is.na(data_vect$humidity)])
# 


# # Create a humidity metrics between 0 and 4 and averaged on a 1-km2 grid
# humid_raster_metrics <- humid_raster
# humid_raster_metrics[humid_raster_metrics$MPH_FR_WGS84 == 522] <- NA
# humid_raster_metrics[humid_raster_metrics$MPH_FR_WGS84 == 51] <- 4
# test_aggregate <- aggregate(humid_raster_metrics, fact=20, fun=mean)
# test_aggregate
# cellSize(test_aggregate)
# 
# data_vect$humidity_metrics <- extract(test_aggregate, data_vect)$MPH_FR_WGS84
# hist(data_vect$humidity_metrics)
# length(data_vect)
# length(data_vect[is.na(data_vect$humidity_metrics)])
# 
# # test_aggregate <- crop(test_aggregate, study_area)
# # test_aggregate <- mask(test_aggregate, study_area)
# 
# plot(test_aggregate, type='continuous',
#      col=rev(hcl.colors(50, "Blues")))
# points(data_vect, pch=20, col="red", cex=2, alpha=0.5)
# plot(study_area, add=T)
# 
# test_aggregate_resample <- resample(test_aggregate, topo_raster, method="bilinear")
# plot(test_aggregate_resample, type='continuous',
#      col=rev(hcl.colors(50, "Blues")))
# points(data_vect, pch=20, col="red", cex=2, alpha=0.5)
# plot(study_area, add=T)







# Altitude
topo_raster <- rast('../data/topology/OpenTopography/output_BE_WGS84.tif')
# imported_raster <- project(x= imported_raster, y =  "EPSG:4326", method = "bilinear")
topo_raster <- crop(topo_raster, study_area)
topo_raster <- mask(topo_raster, study_area)
names(topo_raster) <- "altitude"

print(topo_raster)
cellSize(topo_raster)
plot(topo_raster)
points(data_vect, pch=20, col="red", cex=2, alpha=0.5)
plot(study_area, add=T)
# df_humid <- na.omit(as.data.frame(imported_raster))

data_vect$altitude <- extract(topo_raster, data_vect)$altitude


  
  


# Land use
land_raster <- rast('../data/land_cover/OSO_20210101_RASTER/DATA/OCS_2021_WGS84.tif')
# imported_raster <- project(x= imported_raster, y =  "EPSG:4326", method = "bilinear")
land_raster <- crop(land_raster, study_area)
land_raster <- mask(land_raster, study_area)
land_raster[land_raster$OCS_2021_WGS84 == 0] <- NA
land_raster
cellSize(land_raster)
plot(land_raster, type="classes")
points(data_vect, pch=20, col="red", cex=1, alpha=0.5)
plot(study_area, add=T)
# # names(land_raster) = "land_cover_digit"
# # land_raster_factor <- as.factor(land_raster)
# 
# # cover_legend_all <- c('NA', 'batis denses', 'batis diffus', 'zones ind et com', 'surfaces routes',
# #                   'colza', 'cereales a pailles', 'proteagineux', 'soja', 'tournesol', 'mais', 'riz',
# #                   'tubercules/racines', 'prairies', 'vergers', 'vignes', 'forets de feuillus',
# #                   'forets de coniferes', 'pelouses', 'landes ligneuses', 'surfaces minérales',
# #                   'plages et dunes', 'glaciers ou neiges', 'eau', 'autres')
# # cover_types <- unique(land_raster)
# # cover_legend <- cover_legend_all[1 + cover_types$OCS_2021_WGS84]
# 
# # tar <- levels(land_raster_factor)[[1]]
# # tar$land_cover_name <- cover_legend
# # levels(land_raster_factor) <- data.frame(tar$ID, tar$land_cover_name)
# # names(land_raster_factor) <- "land_cover_name"
# 
# colors = rev(terrain.colors(length(cover_legend)))
# colors[match('NA', cover_legend)] = "white"
# colors[match('eau', cover_legend)] = "blue"
# 
# plot(land_raster_factor, plg=list(cex = 0.7), col=colors)
# points(data_vect, pch=20, col="red", cex=2, alpha=0.5)
# plot(study_area, add=T)
# # df_humid <- na.omit(as.data.frame(imported_raster))
# 
# 
# data_vect$land_use_fine <- extract(land_raster_factor, data_vect)$land_cover_name
# length(data_vect[data_vect$land_use_fine == "eau"])
# 
# 
# # regroup by higher catogories
# 
# # cover_legend <- c('NA', 'batis denses', 'batis diffus', 'zones ind et com', 'surfaces routes',
# #                   'colza', 'cereales a pailles', 'proteagineux', 'soja', 'tournesol', 'mais', 'riz',
# #                   'tubercules/racines', 'prairies', 'vergers', 'vignes', 'forets de feuillus',
# #                   'forets de coniferes', 'pelouses', 'landes ligneuses', 'surfaces minérales',
# #                   'plages et dunes', 'glaciers ou neiges', 'eau', 'autres')
# 
# 
# # match('surfaces routes', cover_legend_all)
# 
# # mat = rbind(c(2, 1), c(3, 1), c(4, 1),
# #             c(6, 5), c(7, 5), c(8, 5), c(9, 5), c(10, 5), c(11, 5), c(12, 5), c(13, 5), c(14, 5), c(15, 5),
# #             c(17, 16), c(19, 16)
# #             )
# mat = rbind(c(2, 1), c(3, 1), c(4, 1),
#             c(6, 5), c(7, 5), c(8, 5), c(9, 5), c(10, 5), c(12, 5), c(13, 5), c(14, 5),
#             c(17, 16), c(19, 16),
#             c(20, 18), c(21, 18)
# )
# land_raster_classified <- classify(land_raster_factor, mat)
# cover_types <- unique(land_raster_classified)
# land_raster_classified_factor <- as.factor(land_raster_classified)
# cover_legend <- cover_legend_all[1 + cover_types$land_cover_name]
# cover_legend[match('batis denses', cover_legend)] <- "impermeable"
# cover_legend[match('colza', cover_legend)] <- "agriculture"
# cover_legend[match('forets de feuillus', cover_legend)] <- "ligneux"
# 
# tar <- levels(land_raster_classified_factor)[[1]]
# tar$land_cover_name <- cover_legend
# levels(land_raster_classified_factor) <- data.frame(tar$ID, tar$land_cover_name)
# names(land_raster_classified_factor) <- "land_cover_name"
# 
# colors = c('white', 'black', 'orange', "#3EBB00", "#8BD000", 'grey',"#E6E600", 'blue')
# plot(land_raster_classified_factor, plg=list(cex = 0.7), col=colors)
# points(data_vect, pch=20, col="red", cex=1, alpha=1)
# plot(study_area, add=T)
# 
# extr <- extract(land_raster_classified_factor, data_vect)
# data_vect$land_use_coarse <- extract(land_raster_classified_factor, data_vect)$land_cover_name
# length(data_vect[data_vect$land_use_coarse == "eau"])


# freq(humid_raster)
hist(data_vect$coordinateUncertaintyInMeters, 50)
barplot(prop.table(table(as.factor(data_vect$humidity))))
hist(data_vect$altitude)
barplot(prop.table(table(as.factor(data_vect$land_use_fine))))
barplot(prop.table(table(as.factor(data_vect$land_use_coarse))))



# land_raster_classified_factor
# cellSize(land_raster_classified_factor)
# land_raster
# cellSize(land_raster)
# humid_raster
# cellSize(humid_raster)
# 
# 
# 
# aggregate_percentage <- function(x) {
#   freq <- sum(x)
#   return(freq)
# }
# 
# test_segregate <- segregate(land_raster_classified_factor)
# test_segregate
# test_aggregate <- aggregate(test_segregate, fact=10, fun=sum)
# test_aggregate
# test_aggregate$`1`
# cellSize(test_aggregate)
# 
# # land_percentage <- zonal(land_raster_classified_factor, humid_raster)




# Aggregate by custom higher categories

mat = rbind(c(2, 1), c(3, 1), c(4, 1),
            c(6, 5), c(7, 5), c(8, 5), c(9, 5), c(10, 5), c(12, 5), c(13, 5), c(14, 5),
            c(19, 17),    # c(17, 16), c(19, 16),
            c(20, 18), c(21, 18)
)
land_raster_classified <- classify(land_raster, mat)
land_raster_classified
unique(land_raster_classified)
colors <- c('black', 'orange', 'grey', 'darkviolet', "#5D8700", "#A5D721", '#EEFFBA', 'blue')
labels <- c('impermeable', 'agriculture_autre', 'riz', 'vigne', 'feuillus', "autres_ligneux", "non_ligneux", 'eau')

pdf("../results/land_use_classified.pdf")
plot(land_raster_classified, type="classes", col=colors, legend=F, axes=F)
sbar(50, xy='topleft', type='bar', below='km')
legend("topright", labels, fill=colors, inset=c(0.01, 0.05))
# points(data_vect, pch=20, col="red", cex=1, alpha=0.5)
plot(study_area, add=T)
dev.off()



# Create percentage in 1km² grid

aggregate_percentage <- function(x) {
  freq <- sum(x) / 100
  # Careful here the denominator depends on the initial resolution and aggregation factor
  return(freq)
}

land_raster_segregate <- segregate(land_raster_classified)
land_raster_segregate
cellSize(land_raster_segregate)
land_raster_aggregate <- aggregate(land_raster_segregate, fact=100, fun=aggregate_percentage)
land_raster_aggregate
cellSize(land_raster_aggregate)
names(land_raster_aggregate) <- labels

pdf("../results/land_use_percentage_water.pdf")
plot(land_raster_aggregate$eau, axes=F, main="Pourcentage de sols occupés en eau")
sbar(50, xy='topleft', type='bar', below='km')
# points(data_vect, pch=20, col="red", cex=1, alpha=0.5)
plot(study_area, add=T)
dev.off()

land_raster_aggregate$riz








# labels <- c('batis denses', 'batis diffus', 'zones ind et com', 'surfaces routes',
#                   'colza', 'cereales a pailles', 'proteagineux', 'soja', 'tournesol', 'mais', 'riz',
#                   'tubercules/racines', 'prairies', 'vergers', 'vignes', 'forets de feuillus',
#                   'forets de coniferes', 'pelouses', 'landes ligneuses', 'surfaces minérales',
#                   'plages et dunes', 'eau')
# 
# land_raster_segregate_fine <- segregate(land_raster)
# land_raster_segregate_fine
# cellSize(land_raster_segregate_fine)
# land_raster_aggregate_fine <- aggregate(land_raster_segregate_fine, fact=100, fun=aggregate_percentage)
# land_raster_aggregate_fine
# cellSize(land_raster_aggregate_fine)
# names(land_raster_aggregate_fine)
# names(land_raster_aggregate_fine) <- labels
# plot(land_raster_aggregate_fine$vignes)
# points(data_vect, pch=20, col="red", cex=1, alpha=0.5)
# plot(study_area, add=T)
# land_raster_aggregate_fine$riz






# Worldclim data
list_files <- list.files("../data/climatic/worldclim/historical/", full.names = TRUE)
bioclim_world <- rast(list_files)
bioclim_world <- crop(bioclim_world, study_area)
bioclim_world <- mask(bioclim_world, study_area)
cellSize(bioclim_world)
names(bioclim_world) <- lapply(names(bioclim_world), function(x) sub(".*m_", "", x))


print(bioclim_world)
plot(bioclim_world)
plot(bioclim_[[1]])
points(data_vect, pch=20, col="red", cex=2, alpha=0.5)
plot(study_area, add=T)





# Bioclim forecast

list_files <- list.files("../data/climatic/worldclim/forecast/", full.names = TRUE)
bioclim_forecast <- c()
for (file in list_files) {
  params <- strsplit(file, "_")[[1]]
  model <- params[4]
  scenario <- params[5]
  horizon <- sub(".tif.*", "", params[6])
  r <- rast(file)
  bio_factors <- lapply(names(r), function(x) sub(".*wc2_", "", x))
  names(r) <- lapply(bio_factors, function(x) paste('bio-', x, '_', model, '_', scenario, '_', horizon, sep=''))
  bioclim_forecast <- append(bioclim_forecast, r)
}
bioclim_forecast <- crop(bioclim_forecast, study_area)
bioclim_forecast <- mask(bioclim_forecast, study_area)
bioclim_forecast
names(bioclim_forecast)






# # DRIAS daily data .nc
# daily_clim_rast = rast('data/DRIAS/prtotAdjust_France_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CNRM-ALADIN63_v2_MF-ADAMONT-SAFRAN-1980-2011_day_19510101-20051231.nc')
# daily_clim_rast_proj = project(daily_clim_rast, study_area)
# # crs(daily_clim_rast) <- crs(humid_raster)
# cellSize(daily_clim_rast)
# daily_clim_rast <- crop(daily_clim_rast, study_area)
# daily_clim_rast <- mask(daily_clim_rast, study_area)
# daily_clim_rast
# plot(daily_clim_rast[[1]])





# # gribdata = nc_open('data/data.nc')
# # gribdata = brick('data/data.nc')
# gribdata = rast('data/data.nc')
# gribdata = rast('../data/adaptor.mars.internal-1680883832.5103214-3700-14-5665f8d5-5fab-4cfb-aacc-1556037bbdee.grib')
# # gribdata = readGDAL('data/adaptor.mars.internal-1680883832.5103214-3700-14-5665f8d5-5fab-4cfb-aacc-1556037bbdee.grib')
# # imported_raster <- raster('data/')
# 
# gribdata <- crop(gribdata, study_area)
# gribdata <- mask(gribdata, study_area)
# 
# print(gribdata)
# # plot(gribdata)
# plot(gribdata$`lai_hv_expver=1_1`)
# points(data_vect, pch=20, col="red", cex=2, alpha=0.5)
# plot(study_area, add=T)
# 
# df_clim <- na.omit(as.data.frame(gribdata))
# colnames(df_clim)[1]
# 
# info <- gdalinfo('data/adaptor.mars.internal-1680883832.5103214-3700-14-5665f8d5-5fab-4cfb-aacc-1556037bbdee.grib')






# rastlist <- list.files(path = "data/2_5m_mean_00s/", pattern='.tif$', all.files= T, full.names= T)
# # rastlist <- list.files(path = "data/2_5m_min_00s/", pattern='.tif$', all.files= T, full.names= T)
# 
# MERRAclim = rast(rastlist)
# MERRAclim <- crop(MERRAclim, study_area)
# MERRAclim <- mask(MERRAclim, study_area)
# # plot(MERRAclim)
# plot(MERRAclim[[1]])
# points(data_vect, pch=20, col="red", cex=2, alpha=0.5)
# plot(study_area, add=T)
# 
# 
# MERRAclim <- resample(MERRAclim, humid_raster, method="bilinear")
# topo_raster <- resample(topo_raster, humid_raster, method="bilinear")
# humid_raster <- resample(humid_raster, humid_raster, method="bilinear")
# MERRAclim <- extend(MERRAclim, study_area)
# topo_raster <- extend(topo_raster, study_area)
# humid_raster <- extend(humid_raster, study_area)
# 
# plot(MERRAclim[[1]])
# points(data_vect, pch=20, col="red", cex=2, alpha=0.5)
# plot(study_area, add=T)
# 
# plot(topo_raster)
# points(data_vect, pch=20, col="red", cex=2, alpha=0.5)
# plot(study_area, add=T)
# 
# plot(humid_raster)
# points(data_vect, pch=20, col="red", cex=2, alpha=0.5)
# plot(study_area, add=T)
# 
# expl_var = rast(list(humid_raster, MERRAclim))
# # expl_var <- MERRAclim
# # plot(expl_var)
# plot(expl_var[[1]])
# points(data_vect, pch=20, col="red", cex=2, alpha=0.5)
# plot(study_area, add=T)












# Pedological data
list_files <- list.files("../data/pedological/soilgrids/", pattern = "\\.tif$",
                         full.names = TRUE)
soil_raster <- rast(list_files)
soil_raster <- crop(soil_raster, study_area)
soil_raster <- mask(soil_raster, study_area)
names(soil_raster) <- lapply(names(soil_raster), function(x) sub(".*grid_", "", x))
names(soil_raster) <- lapply(names(soil_raster), function(x) sub("_5.*", "", x))
# names(soil_raster) <- lapply(names(soil_raster), function(x) gsub("-", "_", x))
names(soil_raster)
soil_raster

cellSize(soil_raster)
plot(soil_raster)
plot(soil_raster$bdod)
points(data_vect, pch=20, col="red", cex=2, alpha=0.5)
plot(study_area, add=T)

length(is.na(soil_raster$`soil_grid_bdod_5-15cm_mean`))
hist(soil_raster$`soil_grid_ocd_5-15cm_mean`)

soil_raster$`soil_grid_bdod_5-15cm_mean`
data_vect$`soil_grid_bdod_5-15cm_mean` <- extract(soil_raster, data_vect)$`soil_grid_bdod_5-15cm_mean`
hist(data_vect$`soil_grid_bdod_5-15cm_mean`)















# Harmonize rasters and save

soil_raster
humid_raster
land_raster_aggregate
topo_raster
bioclim_world
bioclim_forecast

cellSize(soil_raster)
cellSize(humid_raster)
cellSize(land_raster_aggregate)
cellSize(topo_raster)
cellSize(bioclim_world)
cellSize(bioclim_forecast)

ref_sample <- humid_raster



soil_resampled <- resample(soil_raster, ref_sample, method="bilinear")
humid_resampled <- resample(humid_raster, ref_sample, method="bilinear")
land_resampled <- resample(land_raster_aggregate, ref_sample, method="bilinear")
# land_fine_resampled <- resample(land_raster_aggregate_fine, ref_sample, method="bilinear")
topo_resampled <- resample(topo_raster, ref_sample, method="bilinear")
bioclim_resampled <- resample(bioclim_world, ref_sample, method="bilinear")
# bioclim_forecast_resampled <- resample(bioclim_forecast, ref_sample, method="bilinear")

# soil_extended <- extend(soil_resampled, study_area)
# humid_extended <- extend(humid_resampled, study_area)
# land_extended <- extend(land_resampled, study_area)
# topo_extended <- extend(topo_resampled, study_area)
# bioclim_extended <- extend(bioclim_resampled, study_area)

# expl_raster <- c(soil_resampled, humid_resampled, land_resampled, topo_resampled, bioclim_resampled)
# expl_raster
# names(expl_raster)

# writeVector(study_area, '../data/clean/study_area.gpkg', overwrite=T)

# write.csv(occurences, '../data/clean/Limonium_girardianum_occurences.csv')
# write.csv(occurences, '../data/clean/Pancratium_maritimum_occurences.csv')
# write.csv(occurences, '../data/clean/Althenia_filiformis_occurences.csv')

writeRaster(soil_resampled, '../data/clean/soil.tif', overwrite=T)
writeRaster(humid_resampled, '../data/clean/humidity.tif', overwrite=T)
writeRaster(land_resampled, '../data/clean/land_use_percentage.tif', overwrite=T)
# writeRaster(land_fine_resampled, '../data/clean/land_use_fine_percentage.tif', overwrite=T)
writeRaster(topo_resampled, '../data/clean/altitude.tif', overwrite=T)
writeRaster(bioclim_resampled, '../data/clean/bioclimatic.tif', overwrite=T)
writeRaster(bioclim_forecast, '../data/clean/forecast/bioclimatic_forecast.tif', overwrite=T)







