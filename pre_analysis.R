# library(raster)
library(tidyr)
library(biomod2)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(gridExtra)
library(reshape2)
# library(geodata)
library(ade4)
library(ncdf4)
library(terra)
library(corrplot)
library(factoextra)
library(readxl)
library(energy)
library(furrr)
# library(purrr)
library(Matrix)


setwd(dirname(rstudioapi::getSourceEditorContext()$path))


compute_distance_corr_matrix <- function(df){
  df_cross <- expand_grid(x=names(df), y=names(df))
  plan(multisession, workers = 6)
  start_time <- Sys.time()
  df_map <- future_map2_dfr(df_cross$x, df_cross$y,
                            .f=function(x,y) {
                              cor <- dcor(df[, x], df[, y])
                              data.frame(x=x, y=y, cor=cor)
                            }, .options=furrr_options(seed = 1))
  end_time <- Sys.time()
  print(end_time - start_time)
  cor_mat <- acast(df_map, y ~ x, value.var="cor")
  return(cor_mat)
}




# Import clean data
study_area <- vect('../data/clean/study_area.gpkg')
occurences <- read.csv('../data/clean/Limonium_girardianum_occurences.csv')
# occurences <- read.csv('../data/clean/Pancratium_maritimum_occurences.csv')
# occurences <- read.csv('../data/clean/Althenia_filiformis_occurences.csv')
list_files <- c(list.files("../data/clean/", pattern = "\\.tif$",
                         full.names = TRUE))
list_files <- list_files[! grepl("fine", list_files)]
# list_files <- list_files[! grepl("land_use_percentage", list_files)]
# list_files <- list_files[! grepl("altitude", list_files)]
expl_var <- rast(list_files)
expl_var
names(expl_var)



occurences_vect <- vect(occurences, geom=c("x", "y"))
# length(unique((extract(expl_var, occurences_vect, ID=F, cells=T))$cell))
# length(unique((na.omit(extract(expl_var, occurences_vect, ID=F, cells=T)))$cell))

df_presence <- terra::extract(expl_var, occurences_vect, ID=F, cells=T)
nrow(df_presence)
df_presence <- df_presence[!duplicated(df_presence$cell), ]
# df_presence <- subset(df_presence, select = -c(cell) )
nrow(df_presence)
colSums(is.na(df_presence))
df_presence <- na.omit(df_presence)
nrow(df_presence)
# df_presence <- df_presence[sample(nrow(df_presence), 100), ]
df_presence$presence <- 1

df_background <- na.omit(spatSample(expl_var, 800, method='regular', cells=T))   # 800 1600
# df_background <- spatSample(expl_var, 340, method='random', na.rm=T, cells=T)
df_background <- df_background[!df_background$cell %in% df_presence$cell, ]
df_background$presence <- 0
nrow(df_background)

df_sample <- rbind(df_presence, df_background)
nrow(df_sample)





# Correlation matrix
corr_mat <- cor(df_sample[, -which(names(df_sample) == 'cell')])
# pdf("../results/correlation.pdf")
corrplot(corr_mat, order = 'hclust')  # FPC
# dev.off()

dist_corr_mat <- compute_distance_corr_matrix(df_sample[, -which(names(df_sample) == 'cell')])
# pdf("../results/correlation_dist.pdf")
corrplot(dist_corr_mat, is.corr=F, col.lim=c(0,1), order = 'hclust')
# dev.off()

corr_df <- as.data.frame(as.table(corr_mat))
corr_df$Var1 <- as.character(corr_df$Var1)
corr_df$Var2 <- as.character(corr_df$Var2)
corr_df <- corr_df[with(corr_df, order(Var1, Var2)), ]
dist_corr_df <- as.data.frame(as.table(dist_corr_mat))
dist_corr_df$Var1 <- as.character(dist_corr_df$Var1)
dist_corr_df$Var2 <- as.character(dist_corr_df$Var2)
dist_corr_df <- dist_corr_df[with(dist_corr_df, order(Var1, Var2)), ]
plot(abs(corr_df$Freq), dist_corr_df$Freq)
abline(lm(dist_corr_df$Freq ~ abs(corr_df$Freq)), col="red")

all_corr_df <- corr_df
names(all_corr_df)[names(all_corr_df) == 'Freq'] <- 'pearson'
all_corr_df$distance <- dist_corr_df$Freq
ggplot(all_corr_df, aes(pearson, distance, color=Var1)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(intercept = 0, slope = -1) +
  xlim(-1, 1) + ylim(0, 1) +
  coord_fixed()



# PCA
# df_sample[df_sample$presence == 0, 'presence'] <- "background"
# df_sample[df_sample$presence == 1, 'presence'] <- "presence"
pca <- prcomp(df_sample[, -which(names(df_sample) %in% c('presence', 'cell'))], scale.=T)

# outliers?
fviz_pca_ind(pca, label='ind')
fviz_pca_ind(pca, col.ind = "contrib",
             gradient.cols = c("white", "#2E9FDF", "#FC4E07" ))
# outliers <- c('545536', '545537')
outliers <- c('562', '563')
outliers_cells <- df_sample[outliers, ]$cell
(outliers_df <- extract(expl_var, outliers_cells, xy=T)[, c('x', 'y')])
outliers_vect <- vect(outliers_df, geom=c("x", "y"))
plot(expl_var$bio_1)
points(outliers_vect, pch=20, col="red", cex=2, alpha=0.5)

# df_sample <- df_sample[!rownames(df_sample) %in% outliers, ]

fviz_pca_ind(pca, col.ind = "contrib",
             gradient.cols = c("white", "#2E9FDF", "#FC4E07" ))

fviz_pca_ind(pca, col.ind="cos2", geom = "point",
             gradient.cols = c("white", "#2E9FDF", "#FC4E07" ))
fviz_pca_ind(pca, habillage=df_sample$presence, label="none",
             addEllipses=TRUE, ellipse.level=0.99, palette = "Dark2")
# ggsave("../results/PCA_indiv.pdf")
fviz_pca_var(pca, col.var = "steelblue", repel = T)
fviz_pca_var(pca, col.var = "contrib",
             gradient.cols = c("white", "blue", "red"),
             ggtheme = theme_minimal(), repel=T)
# ggsave("../results/PCA_var.pdf")
fviz_pca_biplot(pca, label = "var", habillage=df_sample$presence,
                addEllipses=TRUE, ellipse.level=0.95,
                ggtheme = theme_minimal(), repel=T,
                axes = c(1, 2)
                )
# ggsave("../results/PCA_biplot.pdf")


fviz_eig(pca, addlabels = TRUE)
# ggsave("../results/PCA_eig.pdf")
fviz_contrib(pca, choice = "var", axes = 1:2)
# ggsave("../results/PCA_contrib_1-2.pdf")
fviz_contrib(pca, choice = "var", axes = 1:3)
# ggsave("../results/PCA_contrib_1-3.pdf")
fviz_contrib(pca, choice = "var", axes = 1)
fviz_contrib(pca, choice = "var", axes = 2)
fviz_contrib(pca, choice = "var", axes = 3)
fviz_pca_var(pca, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)




plot(df_sample$bdod, df_sample$presence)




# Compare historical and forecast climatic data

# Plot CO2 emissions for SSP scenarios

CO2_forecast <- read_xlsx('../data/climatic/SSP Public Database/iamc_db.xlsx', n_max=4)
CO2_forecast <- as.data.frame(CO2_forecast)
CO2_forecast <- subset(CO2_forecast, select = -c(Notes))
CO2_forecast <- melt(CO2_forecast, id.vars=c("Model", "Scenario", "Region", "Variable", "Unit"))

CO2_historical <- read_xlsx('../data/climatic/SSP Public Database/iamc_db.xlsx', skip=5, n_max=1)
CO2_historical <- as.data.frame(CO2_historical)
names(CO2_historical)[names(CO2_historical) == 'Scenario (History)'] <- 'Scenario'
CO2_historical <- melt(CO2_historical, id.vars=c("Model", "Scenario", "Region", "Variable", "Unit"))

CO2 <- rbind(CO2_historical, CO2_forecast)
CO2$value <- CO2$value / 1000
CO2$variable <- as.numeric(as.character(CO2$variable))
ggplot(CO2, aes(x=variable, y=value, color=Scenario, group=Scenario)) +
  geom_line(size=1) +
  scale_colour_manual(values=c('#33A43A', '#DB9011', '#D4101A', '#742370', 'grey'),
                      labels=c('SSP1-26', 'SSP2-45', 'SSP3-70', 'SSP5-85', 'Historique'),
                      guide = guide_legend(reverse = TRUE),
                      name = "Scénario") +
  theme(axis.title.x=element_blank(),
        text = element_text(size = 20)) +
  ylab("Gt CO2/année") +
  ggtitle("Emissions totales de CO2")
ggsave("../results/SSP_CO2.pdf")






# historical
list_files <- list.files("../data/climatic/worldclim/historical/", full.names = TRUE)
clim_historical <- rast(list_files)
clim_historical <- crop(clim_historical, study_area)
clim_historical <- mask(clim_historical, study_area)
names(clim_historical) <- lapply(names(clim_historical), function(x) sub(".*m_", "", x))
names(clim_historical) <- lapply(names(clim_historical), function(x) gsub("_", "-", x))

# forecast
clim_forecast <- rast("../data/clean/forecast/bioclimatic_forecast.tif")
clim_forecast


mean_historical <- global(clim_historical, fun=function(x) mean(x, na.rm=T))
colnames(mean_historical) <- c('mean')
mean_historical$factor <- rownames(mean_historical)
rownames(mean_historical) <- NULL
mean_historical[c('model', 'scenario')] = NA
mean_historical$horizon <- 'historical'


mean_forecast <- global(clim_forecast, fun=function(x) mean(x, na.rm=T))
mean_forecast <- cbind.data.frame(mean_forecast, do.call(rbind, strsplit(rownames(mean_forecast), '_')))
rownames(mean_forecast) <- NULL
colnames(mean_forecast) <- c('mean', 'factor', 'model', 'scenario', 'horizon')
mean_forecast

models_scenarios <- expand.grid(unique(mean_forecast$model), unique(mean_forecast$scenario))
mean_clim <- mean_forecast
for (row in 1:nrow(models_scenarios)) {
  mean_historical$model <- models_scenarios[row, "Var1"]
  mean_historical$scenario <- models_scenarios[row, "Var2"]
  mean_clim <- rbind(mean_clim, mean_historical)
}


relative_variations = c()
for (row in 1:nrow(mean_clim)) {
  initial_value <- mean_clim[mean_clim$factor == mean_clim[row, 'factor'] &
                               mean_clim$model == mean_clim[row, 'model'] &
                               mean_clim$scenario == mean_clim[row, 'scenario'] &
                               mean_clim$horizon == 'historical', 'mean']
  relative_variation <- (mean_clim[row, 'mean'] - initial_value) / initial_value
  relative_variations <- append(relative_variations, relative_variation)
}
mean_clim$relative_variation <- relative_variations


mean_clim$horizon <- factor(mean_clim$horizon, levels=c('historical', '2021-2040', '2041-2060', '2061-2080', '2081-2100'))
mean_clim$factor <- factor(mean_clim$factor, levels=c('bio-1', 'bio-2', 'bio-3', 'bio-4', 'bio-5', 'bio-6', 'bio-7', 'bio-8', 'bio-9', 'bio-10', 'bio-11', 'bio-12', 'bio-13', 'bio-14', 'bio-15', 'bio-16', 'bio-17', 'bio-18', 'bio-19'))

mean_clim$label <- ifelse(mean_clim$horizon == '2081-2100' & mean_clim$scenario == 'ssp126',
                          as.character(mean_clim$factor), NA_character_)

ggplot(mean_clim,
       aes(x=horizon, y=mean, color=factor, linetype=scenario, group=interaction(factor, scenario))) +
  geom_line() +
  # geom_point() +
  scale_linetype_manual(values=c("solid", "longdash", "dashed", "dotted"), guide = guide_legend(reverse = TRUE)) +
  # guides(color=guide_legend(ncol=2)) +
  scale_colour_discrete(guide = "none") +
  geom_label_repel(aes(label = label),
                   nudge_x = 1,
                   na.rm = TRUE) +
  scale_y_log10()
ggsave("../results/climatic_evolution.pdf")

ggplot(mean_clim,
       aes(x=horizon, y=relative_variation, color=factor, linetype=scenario, group=interaction(factor, scenario))) +
  geom_line() +
  # geom_point() +
  scale_linetype_manual(values=c("solid", "longdash", "dashed", "dotted"), guide = guide_legend(reverse = TRUE)) +
  # guides(color=guide_legend(ncol=2)) +
  scale_colour_discrete(guide = "none") +
  geom_label_repel(aes(label = label),
                   nudge_x = 1,
                   na.rm = TRUE)
ggsave("../results/climatic_relative_variation.pdf")



get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


factor_names <- c('BIO1 Annual Mean Temperature',
                  'BIO2 Mean Diurnal Range',
                  'BIO3 Isothermality',
                  'BIO4 Temperature Seasonality',
                  'BIO5 Max Temperature of Warmest Month',
                  'BIO6 Min Temperature of Coldest Month',
                  'BIO7 Temperature Annual Range',
                  'BIO8 Mean Temperature of Wettest Quarter',
                  'BIO9 Mean Temperature of Driest Quarter',
                  'BIO10 Mean Temperature of Warmest Quarter',
                  'BIO11 Mean Temperature of Coldest Quarter',
                  'BIO12 Annual Precipitation',
                  'BIO13 Precipitation of Wettest Month',
                  'BIO14 Precipitation of Driest Month',
                  'BIO15 Precipitation Seasonality',
                  'BIO16 Precipitation of Wettest Quarter',
                  'BIO17 Precipitation of Driest Quarter',
                  'BIO18 Precipitation of Warmest Quarter',
                  'BIO19 Precipitation of Coldest Quarter'
)

for (i in 1:19) {
  factor <- paste("bio-", i, sep="")
  p <- ggplot(mean_clim[mean_clim$factor == factor,],
              aes(x=horizon, y=mean, color=scenario, group=scenario)) +
    geom_line() +
    geom_point() +
    scale_colour_manual(values=c('#33A43A', '#DB9011', '#D4101A', '#742370'), guide = guide_legend(reverse = TRUE)) +
    ggtitle(factor_names[i])
  legend <- get_legend(p)
  p <- p + theme(legend.position="none")
  assign(paste("p", i, sep=""), p)
  
}
grid <- grid.arrange(p1, p2, p3, legend, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19,
             nrow=5, ncol=4)
plot(grid)
ggsave("../results/climatic_evolution_by_factor.pdf", grid, scale=2)


for (i in 1:19) {
  factor <- paste("bio-", i, sep="")
  p <- ggplot(mean_clim[mean_clim$factor == factor,],
              aes(x=horizon, y=relative_variation, color=scenario, group=scenario)) +
    geom_line() +
    geom_point() +
    scale_colour_manual(values=c('#33A43A', '#DB9011', '#D4101A', '#742370'), guide = guide_legend(reverse = TRUE)) +
    ggtitle(factor_names[i]) +
    ylim(min(mean_clim$relative_variation), max(mean_clim$relative_variation))
  legend <- get_legend(p)
  p <- p + theme(legend.position="none")
  assign(paste("p", i, sep=""), p)
  
}
grid <- grid.arrange(p1, p2, p3, legend, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19,
                     nrow=5, ncol=4)
plot(grid)
ggsave("../results/climatic_relative_evolution_by_factor.pdf", grid, scale=2)





# Plot forecast
plot(clim_historical$`bio-13`)
clim_forecast_bio1 <- clim_forecast[[ which(grepl('bio-13_', names(clim_forecast))) ]]
v_min = min(clim_forecast_bio1)@ptr$range_min
v_max = max(clim_forecast_bio1)@ptr$range_max
plot(clim_forecast_bio1, range=c(v_min, v_max))
# 






# Correlations in climatic data

# historical
corr_mat_historical <- cor(as.data.frame(clim_historical))
corrplot(corr_mat_historical, order = 'hclust')

# Correlation clusters :
# 1 11 9 10 5 6 8?        <- "temperature absolute"
# 2 7 3 4                 <- "temperature variation"
# 13? 16 14 17 12 18 19   <- "precipitation absolute"
# 15                      <- "precipitation variation"




df_historical <- as.data.frame(clim_historical)
dist_cor_mat_historical <- compute_distance_corr_matrix(df_historical)
corrplot(dist_cor_mat_historical, is.corr=F, col.lim=c(0,1), order='hclust')



# forecast
df_forecast <- as.data.frame(clim_forecast)
names(df_forecast)

df_CNRM.CM6.1_ssp585_2081.2100 <- df_forecast[, grepl("CNRM-CM6-1", names(df_forecast)) &
                                                grepl("ssp585", names(df_forecast)) &
                                                grepl("2081-2100", names(df_forecast))]
names(df_CNRM.CM6.1_ssp585_2081.2100) <- lapply(names(df_CNRM.CM6.1_ssp585_2081.2100),
                                                function(x) gsub("_.*", "", x))

corr_mat_CNRM.CM6.1_ssp585_2081.2100 <- cor(df_CNRM.CM6.1_ssp585_2081.2100)
corrplot(corr_mat_CNRM.CM6.1_ssp585_2081.2100, order = 'hclust')

dist_corr_mat_CNRM.CM6.1_ssp585_2081.2100 <- compute_distance_corr_matrix(df_CNRM.CM6.1_ssp585_2081.2100)
corrplot(dist_corr_mat_CNRM.CM6.1_ssp585_2081.2100, is.corr=F, col.lim=c(0,1), order='hclust')



df_CNRM.CM6.1_ssp245_2041.2060 <- df_forecast[, grepl("CNRM-CM6-1", names(df_forecast)) &
                                                grepl("ssp245", names(df_forecast)) &
                                                grepl("2041-2060", names(df_forecast))]
names(df_CNRM.CM6.1_ssp245_2041.2060) <- lapply(names(df_CNRM.CM6.1_ssp245_2041.2060),
                                                function(x) gsub("_.*", "", x))

corr_mat_CNRM.CM6.1_ssp245_2041.2060 <- cor(df_CNRM.CM6.1_ssp245_2041.2060)
corrplot(corr_mat_CNRM.CM6.1_ssp245_2041.2060, order = 'hclust')

dist_corr_mat_df_CNRM.CM6.1_ssp245_2041.2060 <- compute_distance_corr_matrix(df_CNRM.CM6.1_ssp245_2041.2060)
corrplot(dist_corr_mat_df_CNRM.CM6.1_ssp245_2041.2060, is.corr=F, col.lim=c(0,1), order='hclust')






# Selected variables
selet_expl <- c('bio_1', 'bio_7', 'bio_18', 'bio_13', 'bio_15',
                'phh2o',
                'wetness_index',
                'feuillus')

corr_mat <- cor(df_sample[, which(names(df_sample) %in% selet_expl)])
corr_mat
# pdf("../results/correlation.pdf")
corrplot(corr_mat, order = 'hclust')  # FPC
# dev.off()

dist_corr_mat <- compute_distance_corr_matrix(df_sample[, which(names(df_sample) %in% selet_expl)])
dist_corr_mat
# pdf("../results/correlation_dist.pdf")
corrplot(dist_corr_mat, is.corr=F, col.lim=c(0,1), order = 'hclust')
