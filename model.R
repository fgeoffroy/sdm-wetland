# library(raster)
library(biomod2)
library(ggplot2)
library(viridisLite)
# library(RColorBrewer)
library(gridExtra)
# library(geodata)
library(ade4)
library(ncdf4)
library(terra)
library(corrplot)
library(factoextra)
# library(parallel)
library(doParallel)
# detectCores(logical = TRUE)
# registerDoParallel(makeCluster(6))
# registerDoParallel(cores=6)



setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Import clean data
study_area <- vect('../data/clean/study_area.gpkg')
occurences <- read.csv('../data/clean/Limonium_girardianum_occurences.csv')
occurences <- read.csv('../data/clean/Althenia_filiformis_occurences.csv')
occurences$Limonium.girardianum = 1
occurences$Althenia.filiformis = 1
list_files <- list.files("../data/clean/", pattern = "\\.tif$",
                         full.names = TRUE)
list_files <- list_files[! grepl("fine", list_files)]
all_expl_var <- rast(list_files)
names(all_expl_var)

select_expl <- c('bio_1', 'bio_7', 'bio_18', 'bio_13', 'bio_15',
                'phh2o', 'ocd',
                # 'wetness_index',
                'eau',
                'feuillus')
expl_var <- all_expl_var[[ select_expl ]]
# expl_var$wetness_index = as.factor(expl_var$wetness_index)
plot(expl_var$eau, colNA='black')

occurences_vect <- vect(occurences, geom=c("x", "y"))
# plot(expl_var$eau[expl_var$eau < 25], colNA='black')
# points(occurences_vect, add=T)

hist(expl_var$eau)

# Nb of occurence 
length(occurences_vect)
# Nb of occurences on different cells
length(unique(terra::extract(expl_var, occurences_vect, ID=F, cells=T)$cell))
# Nb of occurences on different cells with no NA
length(unique(na.omit(terra::extract(expl_var, occurences_vect, ID=F, cells=T))$cell))



species_name <- "Limonium.girardianum"
species_name <- "Althenia.filiformis"

# format data
bm_data <- 
  BIOMOD_FormatingData(
    resp.name = species_name,
    resp.var = occurences[species_name],
    resp.xy = occurences[, c('x', 'y')],
    expl.var = expl_var,
    PA.nb.rep = 2,
    PA.nb.absences = 500,
    # PA.nb.absences = nrow(occurences),
    filter.raster = T,
    PA.strategy = 'random'
  )
# pdf("../results/bm_data_saladelle.pdf")
plot(bm_data)
# dev.off()
bm_data

# models options
model_opt <- 
  BIOMOD_ModelingOptions(
    GLM = list(type = 'quadratic', interaction.level = 1),
    GBM = list(n.trees = 1000),
    GAM = list(algo = 'GAM_mgcv', k = 3, interaction.level = 1),
    # GAM = list(algo = 'BAM_mgcv'),
    MARS = list(type = 'quadratic', interaction.level = 1),
    ANN = list(NbCV = 10, maxit = 1000),
  )


# run models

algo_names <- c("GLM", "GBM", "GAM", "CTA", "ANN", "SRE", "MARS", "RF")  # , "FDA"
# model_names <- c("GAM")
registerDoSEQ()
registerDoParallel(cores=min(4, length(algo_names)))
# registerDoParallel(2)

models <- 
  BIOMOD_Modeling(
    bm.format = bm_data,
    # models = c("GLM", "GBM", "RF", "GAM"),
    # models = c("GLM", "GBM", "RF"),
    models = algo_names,
    bm.options = model_opt,
    metric.eval = c("TSS", "ROC"),
    nb.rep = 2,
    data.split.perc = 80,
    var.import = 3,
    modeling.id = "test_1",
    # nb.cpu = 4,
    do.full.models = F,
  )
registerDoSEQ()

save.image("../results/.RData")

load("../results/.RData")



# evaluation metrics
(models_scores <- get_evaluations(models))

dim(models_scores)
dimnames(models_scores)

bm_PlotEvalMean(
  models, 
  group.by = "algo",
  metric.eval = c("ROC", "TSS"),
  dataset = "validation",
  xlim = c(0.5,1),
  ylim = c(0.5,1),
  do.plot = F
)
# ggsave("../results/evaluation_all_models_calibration_saladelle.pdf")
# ggsave("../results/evaluation_all_models_validation_saladelle.pdf")
# ggsave("../results/evaluation_all_models_calibration_althenie.pdf")
# ggsave("../results/evaluation_all_models_validation_althenie.pdf")

# models_scores_wo_all_run <- models_scores[models_scores$run != "allRun", ]
mean_metrics_eval <- aggregate(models_scores$validation,
                               list(models_scores$metric.eval, models_scores$algo),
                               FUN=function(x) mean(x, na.rm=TRUE))
xtabs(x~., mean_metrics_eval)

bm_PlotEvalBoxplot(
  models,
  group.by = c("algo", "PA"),
  dataset = "validation",
)

bm_PlotEvalBoxplot(
  models, 
  group.by = c("algo", "run"),
  dataset = "validation",
)


bm_PlotEvalMean(
  models, 
  group.by = "run",
  metric.eval = c("ROC","TSS"),
  # dataset = "validation",
  xlim = c(0.5,1),
  ylim = c(0.5,1)
)

bm_PlotEvalMean(
  models, 
  group.by = "PA",
  metric.eval = c("ROC","TSS"),
  # dataset = "validation",
  xlim = c(0.5,1),
  ylim = c(0.5,1)
)

bm_PlotEvalMean(
  models, 
  group.by = "full.name",
  metric.eval = c("ROC","TSS"),
  dataset = "validation",
  xlim = c(0.5,1),
  ylim = c(0.5,1)
)



## variable importance
(var_import <- get_variables_importance(models))

## make the mean of variable importance by algorithm
mean_var_import <- aggregate(var_import$var.imp, list(var_import$expl.var, var_import$algo),
                             FUN=mean)
xtabs(x~., mean_var_import)
reordered_var <- reorder(var_import$expl.var, var_import$var.imp, function(x) -mean(x))
ggplot(var_import, aes(x=reordered_var, y=var.imp, fill=reordered_var)) +
  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=20, size=10, color="red", fill="red") +
  theme(legend.position="none", axis.title.x = element_blank()) +
  ylab('Importance')
# ggsave("../results/var_importance_all_models_saladelle.pdf")
# ggsave("../results/var_importance_all_models_althenie.pdf")


bm_PlotVarImpBoxplot (
  models, 
  # group.by = c("expl.var", "run", "algo"),
  do.plot = F
)


# Response curves

algo_names <- c("GLM", "GBM")
algo_names <- c("GLM", "GBM", "GAM", "CTA", "ANN", "SRE", "FDA", "MARS", "RF")
for (algo_name in algo_names) {
  models_names <- BIOMOD_LoadModels(models, algo=algo_name)
  bm_PlotResponseCurves(
    bm.out = models,
    models.chosen = models_names,
    new.env = get_formal_data(models,'expl.var'), 
    show.variables= get_formal_data(models,'expl.var.names'),
    do.bivariate = FALSE,
    fixed.var = 'mean',
    data_species = get_formal_data(models,'resp.var'),
    main = algo_name
  )
}



# ensemble models
# models_names <- get_built_models(models)
# models_names_wo_all <- models_names[ !grepl("all", models_names) ]

ensemble_models <- 
  BIOMOD_EnsembleModeling(
    bm.mod = models,
    # models.chosen = models_names_wo_all,
    # em.by = 'PA+run',
    em.by = 'all',
    metric.eval = c('ROC'),
    em.algo = c('EMwmean', 'EMca', 'EMcv'),
    metric.select =  c('ROC'),
    metric.select.thresh = c(0.95),
    # models.eval.meth = c('TSS','ROC'),
    var.import = 3,
    # nb.cpu = 4,
  )
ensemble_models@em.models_kept


# ensemble variables importance
(ensemble_var_import <- get_variables_importance(ensemble_models))
ensemble_mean_var_import <- aggregate(ensemble_var_import$var.imp, list(ensemble_var_import$expl.var, ensemble_var_import$algo),
                             FUN=mean)
xtabs(x~., ensemble_mean_var_import)
# corrplot(xtabs(x~., ensemble_mean_var_import), is.corr=F)

ensemble_mean_var_import$x <- round(ensemble_mean_var_import$x, 3)
# pdf("../results/var_importance_ensemble_saladelle.pdf")
# pdf("../results/var_importance_ensemble_althenie.pdf")
grid.table(xtabs(x~., ensemble_mean_var_import))
# dev.off()

# ensemble evaluation metrics
(ensemble_models_scores <- get_evaluations(ensemble_models))

# pdf("../results/evaluation_ensemble_saladelle.pdf")
# pdf("../results/evaluation_ensemble_althenie.pdf")
grid.table(ensemble_models_scores[c(6, 7, 8, 9, 10, 11)])
# dev.off()

# only if 2 or more metrics have been chosen to filter models when creating the ensemble:
bm_PlotEvalMean(
  ensemble_models, 
  group.by = "full.name",
  metric.eval = c("ROC","TSS"),
  # dataset = "validation",
  xlim = c(0.5,1),
  ylim = c(0.5,1)
)


# ensemble response curves
gam_eval_strip <- 
  bm_PlotResponseCurves(
    bm.out = ensemble_models,
    models = BIOMOD_LoadModels(ensemble_models, algo="EMwmean"),    # EMca   EMcv
    new.env = get_formal_data(models,'expl.var'), 
    show.variables= get_formal_data(models,'expl.var.names'),
    do.bivariate = F,
    fixed.var = 'mean',
    data_species = get_formal_data(models,'resp.var'),
    main = ''
  )
# ggsave("../results/response_ensemble_saladelle.pdf")
# ggsave("../results/response_ensemble_althenie.pdf")

gam_eval_strip <- 
  bm_PlotResponseCurves(
    bm.out = ensemble_models,
    models = BIOMOD_LoadModels(ensemble_models, algo="EMwmean"),    # EMca   EMcv
    new.env = get_formal_data(models,'expl.var'), 
    show.variables= c("bio_15", "eau"),        # 1  7  13  15  18
    # show.variables= c("bio_18", "feuillus"),        # 1  7  13  15  18
    do.bivariate = T,
    fixed.var = 'mean',
    data_species = get_formal_data(models,'resp.var')
  )


df_inter_raw <- as.data.frame(gam_eval_strip$tab)
c_eau = c()
c_bio = c()
c_EMwmean = c()
for (i in 1:max(df_inter_raw$id)) {
  c_eau <- append(c_eau, df_inter_raw[df_inter_raw$id == i &
                                        df_inter_raw$expl.name == "eau", "expl.val"])
  c_bio <- append(c_bio, df_inter_raw[df_inter_raw$id == i &
                                        df_inter_raw$expl.name == "bio_15", "expl.val"])
  c_EMwmean <- append(c_EMwmean, df_inter_raw[df_inter_raw$id == i &
                                        df_inter_raw$expl.name == "eau", "pred.val"])
}
df_inter <- data.frame(eau=c_eau, bio=c_bio, EMwmean=c_EMwmean)
df_inter_plot <- df_inter[df_inter$eau <= 50 &
                            df_inter$bio >= 30, ]
ggplot(df_inter_plot, aes(x=bio, y=EMwmean, color=eau, group=eau)) +
  geom_line() +
  scale_colour_viridis_c() +
  xlab("BIO15 Coefficient de variation de la précipitation") +
  ylab("Probabilité de présence") +
  labs(colour = "% eau")
# ggsave("../results/response_ensemble_bio15_eau_saladelle.pdf")
# ggsave("../results/response_ensemble_bio15_eau_althenie.pdf")

# ggplot(df_inter_plot, aes(x=eau, y=EMwmean, color=bio, group=bio)) +
#   geom_line() +
#   scale_colour_viridis_c()








# current projections
cellSize(expl_var)
plot(expl_var$eau)
expl_var_proj <- aggregate(expl_var, fact=10, fun=mean)
cellSize(expl_var_proj)
# plot(expl_var_proj$wetness_index)
plot(expl_var_proj$eau, range=c(0, 100))

get_built_models(models)

# !!!!!!!! If we assume strong interventions on water
# expl_var_proj$eau <- expl_var_proj$eau + 10
# expl_var_proj[expl_var_proj$eau > 100] <- 100

registerDoSEQ()

models_proj_current <- 
  BIOMOD_Projection(
    bm.mod = models,
    proj.name = 'Current',
    new.env = expl_var_proj,
    # models.chosen = "all",
    # models.chosen = c("Limonium.girardianum_PA1_RUN1_GLM", "Limonium.girardianum_PA1_RUN1_GBM", "Limonium.girardianum_PA1_RUN1_RF"),
    # metric.binary = "ROC",
    # output.format = ".img",
    # do.stack = FALSE
    on_0_1000 = F,
    build.clamping.mask = T,
  )
# plot(models_proj_current)
plot(unwrap(models_proj_current@proj.out@val)[[1:9]])
plot(rast('./Limonium.girardianum/proj_Current/proj_Current_ClampingMask.tif'))

ensemble_models_proj_current <- 
  BIOMOD_EnsembleForecasting(
    bm.em = ensemble_models,
    bm.proj = models_proj_current,
    metric.binary = "ROC",
    # output.format = ".img",
    # do.stack = FALSE
    on_0_1000 = F,
    build.clamping.mask = T,
  )
plot(ensemble_models_proj_current)
plot(unwrap(ensemble_models_proj_current@proj.out@val)[[1]], range=c(0, 1))
plot(unwrap(ensemble_models_proj_current@proj.out@val)[[2]], range=c(0, 1))
plot(unwrap(ensemble_models_proj_current@proj.out@val)[[3]])
# plot(unwrap(ensemble_models_proj_current@proj.out@val)[[c(1, 3)]])
# plot(unwrap(ensemble_models_proj_current@proj.out@val)[[c(2, 4)]])
# plot(rast('./Limonium.girardianum/proj_Current/proj_Current_ClampingMask.tif'))


# pdf("../results/projection_current_ensemble_EMwmean_saladelle.pdf")
# pdf("../results/projection_current_ensemble_EMwmean_althenie.pdf")
plot(unwrap(ensemble_models_proj_current@proj.out@val)[[1]], range=c(0, 1), axes=F,
     col=viridis(100))
sbar(50, xy='topleft', type='bar', below='km')
plot(study_area, add=T)
# dev.off()


# pdf("../results/projection_current_ensemble_EMca_saladelle.pdf")
# pdf("../results/projection_current_ensemble_EMca_althenie.pdf")
plot(unwrap(ensemble_models_proj_current@proj.out@val)[[2]], range=c(0, 1), axes=F,
     col=viridis(100))
sbar(50, xy='topleft', type='bar', below='km')
plot(study_area, add=T)
# dev.off()


# pdf("../results/projection_current_ensemble_EMcv_saladelle.pdf")
# pdf("../results/projection_current_ensemble_EMcv_althenie.pdf")
plot(unwrap(ensemble_models_proj_current@proj.out@val)[[3]], axes=F,
     col=cividis(100))
sbar(50, xy='topleft', type='bar', below='km')
plot(study_area, add=T)
# dev.off()








link <- ensemble_models_proj_current@proj.out@link
bin_current <- rast(Filter(function(x) grepl('bin', x), link))
# writeRaster(bin_current, '../results/projection_bin_current_ensemble_saladelle.tif', overwrite=T)
# writeRaster(bin_current, '../results/projection_bin_current_ensemble_intervention_saladelle.tif', overwrite=T)
# writeRaster(bin_current, '../results/projection_bin_current_ensemble_althenie.tif', overwrite=T)
# writeRaster(bin_current, '../results/projection_bin_current_ensemble_intervention_althenie.tif', overwrite=T)
names(bin_current)
plot(bin_current)

# pdf("../results/projection_bin_current_ensemble_EMwmean_saladelle.pdf")
# pdf("../results/projection_bin_current_ensemble_EMwmean_althenie.pdf")
plot(bin_current[[1]], axes=F)
sbar(50, xy='topleft', type='bar', below='km')
plot(study_area, add=T)
# dev.off()


# pdf("../results/projection_bin_current_ensemble_EMca_saladelle.pdf")
# pdf("../results/projection_bin_current_ensemble_EMca_althenie.pdf")
plot(bin_current[[2]], axes=F)
sbar(50, xy='topleft', type='bar', below='km')
plot(study_area, add=T)
# dev.off()




# link_intervention <- ensemble_models_proj_current@proj.out@link
# bin_current_intervention <- rast(Filter(function(x) grepl('bin', x), link_intervention))
# writeRaster(bin_current_intervention, '../results/projection_bin_current_ensemble_intervention_althenie.tif', overwrite=T)










# forecast projections
expl_var_proj_forecast_wo_clim <- expl_var_proj[[ which(! grepl('bio', names(expl_var_proj))) ]]


clim_forecast <- rast("../data/clean/forecast/bioclimatic_forecast.tif")

GCMs <- c('CNRM-CM6-1')
horizons <- c('2021-2040', '2041-2060', '2061-2080', '2081-2100')
scenarios <- c('ssp585', 'ssp370', 'ssp245', 'ssp126')

# scenarios <- c('ssp585')
# horizons <- c('2081-2100')
forecasts <- expand.grid(GCMs, horizons, scenarios)

registerDoSEQ()
rast_proj_forecast <- c()
rast_bin_proj_forecast <- c()
titles <- c()
for(i in 1:nrow(forecasts)) {
  GCM <- forecasts[i, "Var1"]
  horizon <- forecasts[i, "Var2"]
  scenario <- forecasts[i, "Var3"]
  name_forecast <- paste(GCM, scenario, horizon, sep="_")
  titles <- append(titles, paste(scenario, horizon, sep=' '))

  clim_var <- clim_forecast[[which(grepl(GCM, names(clim_forecast)) &
                                grepl(scenario, names(clim_forecast)) &
                                grepl(horizon, names(clim_forecast)))]]
  names(clim_var) <- lapply(names(clim_var), function(x) gsub("_.*", "", x))
  names(clim_var) <- lapply(names(clim_var), function(x) gsub("-", "_", x))
  clim_var <- clim_var[[ which(names(clim_var) %in% names(expl_var_proj)) ]]
  clim_var <- resample(clim_var, expl_var_proj_forecast_wo_clim, method="bilinear")
  expl_var_proj_forecast <- c(clim_var, expl_var_proj_forecast_wo_clim)
  
  models_proj_forecast <- 
    BIOMOD_Projection(
      bm.mod = models,
      new.env = expl_var_proj_forecast,
      proj.name = name_forecast,
      # models.chosen = "all",
      # models.chosen = c("Limonium.girardianum_PA1_RUN1_GLM", "Limonium.girardianum_PA1_RUN1_GBM", "Limonium.girardianum_PA1_RUN1_RF"),
      # metric.binary = "ROC",
      # output.format = ".img",
      # do.stack = FALSE
      on_0_1000 = F,
      build.clamping.mask = T,
    )
  # plot(unwrap(models_proj_forecast@proj.out@val)[[1:9]], colNA='black')
  # 
  ensemble_models_proj_forecast <- 
    BIOMOD_EnsembleForecasting(
      bm.em = ensemble_models,
      bm.proj = models_proj_forecast,
      metric.binary = "ROC",
      # output.format = ".img",
      # do.stack = FALSE
      on_0_1000 = F,
      build.clamping.mask = T,
    )
  
  # path <- paste('./Limonium.girardianum/proj_CNRM-CM6-1_ssp126_2021-2040/', name_forecast,
  #               '/proj_', name_forecast, '_Limonium.girardianum.tif', sep='')
  rast <- unwrap(ensemble_models_proj_forecast@proj.out@val)
  names(rast) <- lapply(names(rast), function(x) paste(x, name_forecast, sep="_"))
  rast_proj_forecast <- append(rast_proj_forecast, rast)
  
  link <- ensemble_models_proj_forecast@proj.out@link
  bin_rast <- rast(Filter(function(x) grepl('bin', x), link))
  names(bin_rast) <- lapply(names(bin_rast), function(x) paste(x, name_forecast, sep="_"))
  rast_bin_proj_forecast <- append(rast_bin_proj_forecast, bin_rast)
}

# writeRaster(rast_proj_forecast, '../results/forecast.tif', overwrite=T)
# writeRaster(rast_proj_forecast, '../results/forecast_wet_intervention.tif', overwrite=T)

# writeRaster(rast_bin_proj_forecast, '../results/projection_bin_forecast_ensemble_saladelle.tif', overwrite=T)
# writeRaster(rast_bin_proj_forecast, '../results/projection_bin_forecast_ensemble_intervention_saladelle.tif', overwrite=T)
# writeRaster(rast_bin_proj_forecast, '../results/projection_bin_forecast_ensemble_althenie.tif', overwrite=T)
# writeRaster(rast_bin_proj_forecast, '../results/projection_bin_forecast_ensemble_intervention_althenie.tif', overwrite=T)

# probabilistic projections
rast_proj_forecast
names(rast_proj_forecast)
# plot(rast_proj_forecast)

rast_proj_forecast_EMwmean <- rast_proj_forecast[[ which(grepl('EMwmean', names(rast_proj_forecast))) ]]
names(rast_proj_forecast_EMwmean)
plot(rast_proj_forecast_EMwmean, range=c(0, 1))

# pdf('../results/projection_forecast_ensemble_EMwmean_saladelle.pdf')
# pdf('../results/projection_forecast_ensemble_EMwmean_althenie.pdf')
plot(rast_proj_forecast_EMwmean, range=c(0, 1), col=viridis(100), main=titles, axes=F, 
     fun=function() plot(study_area, add=T), legend=F, mar=c(0, 0, 0, 0))
# dev.off()
# pdf('../results/projection_forecast_legend.pdf')
image(1, seq(from=0, to=1, by=0.01), t(seq_along(seq(from=0, to=1, by=0.01))), col=viridis(101), axes=FALSE)
axis(4, c(0, 1), cex=15)
# dev.off()


rast_proj_forecast_EMca <- rast_proj_forecast[[ which(grepl('EMca', names(rast_proj_forecast))) ]]
names(rast_proj_forecast_EMca)
plot(rast_proj_forecast_EMca, range=c(0, 1))

rast_proj_forecast_EMwcv <- rast_proj_forecast[[ which(grepl('EMcv', names(rast_proj_forecast))) ]]
names(rast_proj_forecast_EMwcv)
v_min = min(rast_proj_forecast_EMwcv)@ptr$range_min
v_max = max(rast_proj_forecast_EMwcv)@ptr$range_max
# pdf('../results/projection_forecast_ensemble_EMcv_saladelle.pdf')
# pdf('../results/projection_forecast_ensemble_EMcv_althenie.pdf')
plot(rast_proj_forecast_EMwcv, range=c(v_min, v_max), col=cividis(100), main=titles, axes=F)
# dev.off()



# binary projections
rast_bin_proj_forecast
names(rast_bin_proj_forecast)

diff_pixel_rasters = c()
titles_EMwmean <- c()
titles_EMca <- c()
for (i in 1:nrow(forecasts)) {
  GCM <- forecasts[i, "Var1"]
  horizon <- forecasts[i, "Var2"]
  scenario <- forecasts[i, "Var3"]
  
  range_change <- 
    BIOMOD_RangeSize(
      bin_current,
      rast_bin_proj_forecast[[ which(grepl(GCM, names(rast_bin_proj_forecast)) &
                                       grepl(scenario, names(rast_bin_proj_forecast)) &
                                       grepl(horizon, names(rast_bin_proj_forecast))) ]]
    )
  r <- range_change$Diff.By.Pixel
  # trick to have all levels on plots
  r[1] <- -2
  r[2] <- -1
  r[3] <- 0
  r[4] <- 1
  r[5] <- -1
  r[6] <- -1
  r[7] <- -1
  r <- as.factor(r)
  for (j in 1:nlyr(r)) {
    cats <- data.frame(ID=-2:1, label=c("lost", "pres", "abs", "gain"))
    levels(r[[j]]) <- cats
  }
  names(r) <- c(paste(GCM, horizon, scenario, "binary-EMwmean", sep='_'),
                paste(GCM, horizon, scenario, "binary-EMca", sep='_'))
  
  diff_pixel_rasters <- append(diff_pixel_rasters, r)
  titles_EMwmean <- append(titles_EMwmean, paste(scenario, ' ', horizon, '\n', round(range_change$Compt.By.Models[1, 7], 0), '%', sep=''))
  titles_EMca <- append(titles_EMca, paste(scenario, ' ', horizon, '\n', round(range_change$Compt.By.Models[2, 7], 0), '%', sep=''))
  
  # print(range_change$Compt.By.Models)
}


# diff_pixel_rasters <- as.factor(diff_pixel_rasters)
# tar <- levels(diff_pixel_rasters)[[1]]
# tar$ <- cover_legend
# levels(land_raster_factor) <- data.frame(tar$ID, tar$land_cover_name)
# names(land_raster_factor) <- "land_cover_name"
diff_pixel_rasters
names(diff_pixel_rasters)


diff_pixel_rasters_EMwmean <- diff_pixel_rasters[[ which(grepl('EMwmean', names(diff_pixel_rasters))) ]]
names(diff_pixel_rasters_EMwmean)
# plot(diff_pixel_rasters_EMwmean, col=rainbow(4))
# colors <- data.frame(id=-2:1,
#                      col=c("red", "orange", "yellow", "green"))
colors <- c('red', '#A5D721', '#F2F2F2', '#5D8700')
# labels <- c("lost", "pres", "abs", "gain")
plot(diff_pixel_rasters_EMwmean)
# pdf('../results/projection_bin_forecast_ensemble_EMwmean_saladelle.pdf')
# pdf('../results/projection_bin_forecast_ensemble_EMwmean_althenie.pdf')
plot(diff_pixel_rasters_EMwmean, type="classes", axes=F, legend=F, main=titles_EMwmean,
     col=colors, fun=function() plot(study_area, add=T), mar=c(0, 0, 2.6, 0))
# dev.off()

pdf('../results/projection_bin_forecast_legend.pdf')
colors <- c('#5D8700', '#A5D721', '#F2F2F2', 'red')
labels <- c("gain", "prés", "abs", "perdu")
# plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", labels, fill=colors)
# dev.off()


diff_pixel_rasters_EMca <- diff_pixel_rasters[[ which(grepl('EMca', names(diff_pixel_rasters))) ]]
names(diff_pixel_rasters_EMca)
# plot(diff_pixel_rasters_EMca, col=rainbow(4))
colors <- data.frame(id=-2:1,
                     col=c("red", "orange", "yellow", "green"))
colors <- c('red', 'orange', 'yellow', 'green')
labels <- c("lost", "pres", "abs", "gain")
plot(diff_pixel_rasters_EMca)
plot(diff_pixel_rasters_EMca, type="classes", col=colors, legend=F, main=titles_EMca)
legend("topright", labels, fill=colors, inset=c(0.01, 0.05))





# Intervention on water



# current binary projections

range_change <- 
  BIOMOD_RangeSize(
    rast('../results/projection_bin_current_ensemble_saladelle.tif'),
    rast('../results/projection_bin_current_ensemble_intervention_saladelle.tif')
    # rast('../results/projection_bin_current_ensemble_althenie.tif'),
    # rast('../results/projection_bin_current_ensemble_intervention_althenie.tif')
  )
range_change$Compt.By.Models

# pdf('../results/projection_bin_current_intervention_ensemble_EMwmean_saladelle.pdf')
# pdf('../results/projection_bin_current_intervention_ensemble_EMwmean_althenie.pdf')
colors <- c('red', '#A5D721', '#F2F2F2', '#5D8700')
plot(range_change$Diff.By.Pixel[[1]], col=colors, legend=F, axes=F,
     main=paste(round(range_change$Compt.By.Models[1, 7], 0), "%"))
sbar(50, xy='topleft', type='bar', below='km')
colors <- c('#5D8700', '#A5D721', '#F2F2F2', 'red')
labels <- c("gain", "prés", "abs", "perdu")
legend("topright", labels, fill=colors)
plot(study_area, add=T)
# dev.off()




# forecast binary projections

bin_forecast <- rast('../results/projection_bin_forecast_ensemble_saladelle.tif')
bin_forecast_intervention <- rast('../results/projection_bin_forecast_ensemble_intervention_saladelle.tif')
# bin_forecast <- rast('../results/projection_bin_forecast_ensemble_althenie.tif')
# bin_forecast_intervention <- rast('../results/projection_bin_forecast_ensemble_intervention_althenie.tif')

diff_pixel_rasters <- c()
titles_EMwmean <- c()
titles_EMca <- c()
for (i in 1:nrow(forecasts)) {
  GCM <- forecasts[i, "Var1"]
  horizon <- forecasts[i, "Var2"]
  scenario <- forecasts[i, "Var3"]
  
  
  
  r1 <-
    bin_forecast[[ which(grepl(GCM, names(bin_forecast)) &
                                       grepl(scenario, names(bin_forecast)) &
                                       grepl(horizon, names(bin_forecast))) ]]
  
  r2 <-
    bin_forecast_intervention[[ which(grepl(GCM, names(bin_forecast_intervention)) &
                           grepl(scenario, names(bin_forecast_intervention)) &
                           grepl(horizon, names(bin_forecast_intervention))) ]]
  # trick to have all levels on plots
  r1[1] <- 0
  r1[2] <- 0
  r1[3] <- 1
  r1[4] <- 1
  r2[1] <- 0
  r2[2] <- 0
  r2[3] <- 1
  r2[4] <- 1
  
  range_change <- 
    BIOMOD_RangeSize(
      r1,
      r2
    )
  r <- range_change$Diff.By.Pixel
  # trick to have all levels on plots
  r[1] <- -2
  r[2] <- -1
  r[3] <- 0
  r[4] <- 1
  r[5] <- -1
  r[6] <- -1
  r[7] <- -1
  
  r <- as.factor(r)
  for (j in 1:nlyr(r)) {
    cats <- data.frame(ID=-2:1, label=c("lost", "pres", "abs", "gain"))
    levels(r[[j]]) <- cats
  }
  names(r) <- c(paste(GCM, horizon, scenario, "binary-EMwmean", sep='_'),
                paste(GCM, horizon, scenario, "binary-EMca", sep='_'))
  
  diff_pixel_rasters <- append(diff_pixel_rasters, r)
  titles_EMwmean <- append(titles_EMwmean, paste(scenario, ' ', horizon, '\n', round(range_change$Compt.By.Models[1, 7], 0), '%', sep=''))
  titles_EMca <- append(titles_EMca, paste(scenario, ' ', horizon, '\n', round(range_change$Compt.By.Models[2, 7], 0), '%', sep=''))
  
  # print(range_change$Compt.By.Models)
}

diff_pixel_rasters
names(diff_pixel_rasters)


diff_pixel_rasters_EMwmean <- diff_pixel_rasters[[ which(grepl('EMwmean', names(diff_pixel_rasters))) ]]
names(diff_pixel_rasters_EMwmean)
# plot(diff_pixel_rasters_EMwmean, col=rainbow(4))
# colors <- data.frame(id=-2:1,
#                      col=c("red", "orange", "yellow", "green"))
colors <- c('red', '#A5D721', '#F2F2F2', '#5D8700')
# labels <- c("lost", "pres", "abs", "gain")
plot(diff_pixel_rasters_EMwmean)
# pdf('../results/projection_bin_forecast_intervention_ensemble_EMwmean_saladelle.pdf')
# pdf('../results/projection_bin_forecast_intervention_ensemble_EMwmean_althenie.pdf')
plot(diff_pixel_rasters_EMwmean, type="classes", axes=F, legend=F, main=titles_EMwmean,
     col=colors, fun=function() plot(study_area, add=T), mar=c(0, 0, 2.6, 0))
# dev.off()
