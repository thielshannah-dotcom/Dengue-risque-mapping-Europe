# BING-F432 - Spatial and molecular epidemiology (ULB, 2025-2026)
# Project 1: Risk mapping of an arbovirus (here, Dengue virus) in Europe
# Pair 10 (Alexis Biefnot and Hannah Thiels)

# This script performs ecological niche modeling (ENM) using the
# Boosted Regression Trees (BRT) approach to map the risk of local
# Dengue virus circulation mostly transmitted by Aedes albopictus in Europe.
#
# Method inspired by:
#   - Elith et al. (2008) "A working guide to boosted regression trees"
#   - Valavi et al. (2019) "blockCV: spatial and environmental blocking for k-fold CV"
#   - Klitting et al. (2022) - LASV analyses (TP4-5 from the BING-F432 course)
#   - Serres et al. (2025) "Escalating human exposure to tropical mosquito-borne viruses in Europe"
#
# Code adapted from the TP4-5 course tutorial (S. Dellicour, 2025-2026)
# Professor's GitHub: https://github.com/sdellicour


# Create the outputs/ folder if needed after setting the project working directory

# 0. Package installation and loading 

# (uncomment the install.packages lines if needed)

if (!require(raster))       install.packages("raster")
if (!require(sf))           install.packages("sf")
if (!require(sp))           install.packages("sp")
if (!require(RColorBrewer)) install.packages("RColorBrewer")
if (!require(blockCV))      install.packages("blockCV")
if (!require(dismo))        install.packages("dismo")
if (!require(gbm))          install.packages("gbm")
if (!require(ncf))          install.packages("ncf")

library(raster)
library(sf)
library(sp)
library(RColorBrewer)
library(blockCV)
library(dismo)
library(gbm)
library(ncf)


# 1. Data loading 

setwd("~/Desktop/BING-F432_Travail_Groupe/Projet1_Risk_Mapping") # adapt to the user's directory
if (!dir.exists("outputs")) dir.create("outputs")

# 1.1 Load presence points (arbovirus occurrences)
occurrences <- read.csv("data/occurrences/Arbovirus_risk_mapping_simulated_dataset_2-3.csv", header = TRUE)
head(occurrences)
str(occurrences)
cat("Nombre de points de présence :", nrow(occurrences), "\n")

# 1.2 Load the coastline shapefile for the study area
#     This file will be used as a basemap for visualizations
coasts <- shapefile("data/shapefiles/Coasts_study_area_shapefile.shp")

# 1.3 Load the 19 environmental variable rasters
#     Climatic variables (4 seasons x 3 variables = 12 rasters):
#       - Temperature, precipitation, relative humidity
#     Land-cover variables (6 rasters):
#       - Primary forests, secondary forests, primary non-forest areas,
#         secondary non-forest areas, croplands, pastures/rangeland
#     Demographic variable (1 raster):
#       - Human population density (log10)

envVariableNames <- c(
  "temperature_spring", "temperature_summer", "temperature_inFall", "temperature_winter",
  "precipitation_spring", "precipitation_summer", "precipitation_inFall", "precipitation_winter",
  "relative_humidity_spring", "relative_humidity_summer", "relative_humidity_inFall", "relative_humidity_winter",
  "primary_forest_areas", "primary_non-forest_areas",
  "secondary_forest_areas", "secondary_non-forest_areas",
  "croplands_all_categories", "managed_pasture_and_rangeland",
  "human_pop_density_log10"
)

envVariables <- list()
for (i in 1:length(envVariableNames)) {
  fileName <- paste0("data/rasters/Raster_", envVariableNames[i], ".asc")
  envVariables[[i]] <- raster(fileName)
  names(envVariables[[i]]) <- envVariableNames[i]
} # environmental variable rasters are loaded one by one

cat("Nombre de variables environnementales chargées :", length(envVariables), "\n") # Check: 19 variables?
cat("Résolution :", res(envVariables[[1]]), "degrés\n")
cat("Étendue :", as.vector(extent(envVariables[[1]])), "\n")

# Check whether, and how many, points fall outside the raster
vals_check <- raster::extract(envVariables[[1]], 
                              occurrences[, c("longitude", "latitude")])

n_na <- sum(is.na(vals_check))

cat("Points hors raster (NA) :", n_na, "/", nrow(occurrences), "\n")

if (n_na > 0) {
  warning("Certains points de présence tombent hors du raster — vérification nécessaire.")
  
  # Identify the indices of problematic points
  idx_na <- which(is.na(vals_check))
  
  cat("Indices des points concernés :\n")
  print(idx_na)
  
  cat("Coordonnées des points hors raster :\n")
  print(occurrences[idx_na, ])
  
  # Cleaning: remove points outside the raster from the original dataframe
  occurrences <- occurrences[!is.na(vals_check), ]
  
  cat("Nombre de points conservés après nettoyage :", 
      nrow(occurrences), "\n")
}

# 2. Exploratory data visualization 
# 2.1 Spatial exploration of occurrences: map and coordinate distribution 

pdf("outputs/Fig1_exploration_occurrences.pdf", width = 14, height = 5)
par(mfrow = c(1, 3), oma = c(0, 0, 2, 0), mar = c(4, 4, 2, 1),
    lwd = 0.3, col = "gray30", col.axis = "gray30", fg = "gray30")
plot(coasts, col = "gray95", border = "gray50", lwd = 0.5,
     main = "Localisation des présences")
points(occurrences$longitude, occurrences$latitude,
       pch = 16, cex = 1.0, col = rgb(1, 0, 0, 0.6))
points(occurrences$longitude, occurrences$latitude,
       pch = 1, cex = 1.0, col = "black", lwd = 0.3)
hist(occurrences$longitude, breaks = 25, col = "steelblue",
     border = "white", main = "Longitudes",
     xlab = "Longitude (°)", ylab = "Fréquence")
abline(v = mean(occurrences$longitude), col = "red", lty = 2, lwd = 1.5)
hist(occurrences$latitude, breaks = 25, col = "darkorange",
     border = "white", main = "Latitudes",
     xlab = "Latitude (°)", ylab = "Fréquence")
abline(v = mean(occurrences$latitude), col = "red", lty = 2, lwd = 1.5)
mtext(paste0("Exploration des données (n = ", nrow(occurrences), ")"),
      side = 3, outer = TRUE, cex = 1.1, font = 2)
dev.off()

# 2.2 Map of presence points over the 19 environmental variables 

# Function to choose the color palette according to variable type
get_cols <- function(name) {
  if (grepl("temperature", name))      return(colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100))  # blue = cold, red = warm
  if (grepl("precipitation", name))    return(colorRampPalette(brewer.pal(9, "YlGnBu"))(100))       # hydrological convention
  if (grepl("humidity", name))         return(colorRampPalette(brewer.pal(9, "PuBu"))(100))         # distinct from precipitation
  if (grepl("primary_forest", name))   return(colorRampPalette(brewer.pal(9, "Greens"))(100))       # dense canopy green
  if (grepl("primary_non", name))      return(colorRampPalette(brewer.pal(9, "YlGn"))(100))         # open vegetation
  if (grepl("secondary_forest", name)) return(colorRampPalette(brewer.pal(9, "BuGn"))(100))         # secondary forest (blue-green hue)
  if (grepl("secondary_non", name))    return(colorRampPalette(brewer.pal(9, "YlGn"))(100))         # open vegetation (variant)
  if (grepl("cropland", name))         return(colorRampPalette(brewer.pal(9, "YlOrBr"))(100))       # agricultural hues (cultivated land)
  if (grepl("pasture", name))          return(colorRampPalette(brewer.pal(9, "Oranges"))(100))      # pastures (warm hue distinct from croplands)
  if (grepl("human", name))            return(colorRampPalette(brewer.pal(9, "Purples"))(100))      # urban demographic convention, contrasting with red points
  return(colorRampPalette(brewer.pal(9, "Greys"))(100))                                             # neutral fallback
}

# Readable titles for each variable
envVariableTitles <- c(
  "Temp. printemps", "Temp. été",
  "Temp. automne", "Temp. hiver",
  "Précip. printemps", "Précip. été",
  "Précip. automne", "Précip. hiver",
  "Humidité printemps", "Humidité été",
  "Humidité automne", "Humidité hiver",
  "Forêts primaires", "Non-forêt primaire",
  "Forêts secondaires", "Non-forêt secondaire",
  "Croplands", "Pâturages & rangeland",
  "Pop. humaine (log10)"
)

pdf("outputs/Fig2_19variables.pdf", width = 22, height = 14)
par(mfrow = c(4, 5), oma = c(0, 0, 3, 0), mar = c(0.3, 0.3, 1.5, 2.5),
    lwd = 0.2, col = "gray30", col.axis = "gray30", fg = "gray30")

for (i in seq_along(envVariableNames)) {
  plot(envVariables[[i]], col = get_cols(envVariableNames[i]),
       ann = FALSE, legend = FALSE, axes = FALSE, box = FALSE)
  plot(coasts, add = TRUE, border = "gray50", lwd = 0.3)
  # Occurrence points in red with a black outline to stand out on any palette
  points(occurrences$longitude, occurrences$latitude,
         pch = 21, cex = 0.5, col = "black", bg = "red", lwd = 0.3)
  mtext(envVariableTitles[i], side = 3, line = 0, cex = 0.7, col = "gray30", font = 2)
  plot(envVariables[[i]], col = get_cols(envVariableNames[i]),
       legend.only = TRUE, add = TRUE, legend.width = 0.5, legend.shrink = 0.3,
       smallplot = c(0.83, 0.85, 0.1, 0.9),
       axis.args = list(cex.axis = 0.55, lwd = 0, lwd.tick = 0.2,
                        col.tick = "gray30", tck = -0.8, col.axis = "gray30",
                        line = 0, mgp = c(0, 0.2, 0)))
}
# The 20th panel remains empty (4x5 grid = 20 panels, but 19 variables)
plot.new()

mtext("19 variables environnementales - Points rouges = occurrences",
      side = 3, outer = TRUE, cex = 1.3, font = 2, line = 0.5)
dev.off()
cat("Figure 2 sauvegardée : outputs/Fig2_19variables.pdf\n")

# 2.3 Exploration of environmental conditions at presence sites

env_at_presences <- raster::extract(stack(envVariables),
                                    occurrences[, c("longitude", "latitude")])
env_at_presences <- as.data.frame(env_at_presences)

pdf("outputs/Fig3_valeurs_env_aux_presences.pdf", width = 22, height = 10)
par(mfrow = c(4, 5), oma = c(0, 0, 3, 0), mar = c(3, 3, 2, 0.5),
    lwd = 0.3, col = "gray30", col.axis = "gray30", fg = "gray30")
for (i in seq_along(envVariableNames)) {
  vals <- env_at_presences[, i]
  vals <- vals[!is.na(vals)]
  if (length(vals) > 0) {
    hist(vals, breaks = 20, col = "steelblue", border = "white",
         main = envVariableNames[i], xlab = "", ylab = "", cex.main = 0.75)
    abline(v = median(vals), col = "red", lty = 2, lwd = 1.5)
  } else { plot.new() }
}
plot.new()
mtext("Valeurs aux sites de présence (rouge = médiane)",
      side = 3, outer = TRUE, cex = 1.2, font = 2, line = 0.5)
dev.off()

# visualization of environmental variables only at presence points
# Used to check the consistency and quality of extracted data (plausible values, no anomalies)
# For example, this shows that the primary forest environmental variable will probably not be explanatory. 

# 3. Background and pseudo-absence preparation 

# 3.1 Create a template raster for the study area
#     Each terrestrial cell (non-NA) receives the value 1
template <- envVariables[[1]]
template[!is.na(template[])] <- 1

# 3.2 Identify presence cells (cells containing at least 1 occurrence)
cells_presence <- unique(
  raster::extract(template, occurrences[, c("longitude", "latitude")], cellnumbers = TRUE)
)[, "cells"]
cat("Nombre de cellules de présence uniques :", length(cells_presence), "\n")

# 3.3 Define the background: all terrestrial cells EXCEPT presence cells
#     The background represents all cells from which pseudo-absences (PA)
#     will be sampled randomly.
#     Note: in an ideal approach, the background could be restricted to cells
#     where the vector (mosquito) is present, but this information is not available here.
background <- template
background[(1:length(background[])) %in% cells_presence] <- NA
cells_background <- which(!is.na(background[]))
cat("Nombre de cellules dans le background :", length(cells_background), "\n")

# 4. PA generation

# Justification for the number of replicates:
# Each replicate involves a distinct random PA sampling.
# 30 replicates make it possible to assess result robustness to the stochasticity
# of the BRT algorithm and pseudo-absence sampling.

number_of_replicates <- 30
dataframes <- list() # will be used to generate 1 dataframe per replicate

set.seed(42)  # for reproducibility regardless of the user

for (i in 1:number_of_replicates) {
  
  # Coordinates of presence cells (centroids)
  presences <- xyFromCell(template, cells_presence)
  
  # Random pseudo-absence sampling (1:1 ratio with presences)
  cells_pseudo_absence <- sample(cells_background, length(cells_presence), replace = FALSE) # random sampling without replacement
  pseudo_absences <- xyFromCell(template, cells_pseudo_absence)
  
  # Data frame construction: response (1/0), coordinates, environmental values
  data <- rbind(
    cbind(rep(1, nrow(presences)), presences),
    cbind(rep(0, nrow(pseudo_absences)), pseudo_absences)
  )
  colnames(data) <- c("response", "longitude", "latitude")
  
  # Extraction of environmental values at each location
  data <- cbind(data, raster::extract(stack(envVariables), data[, c(2, 3)]))
  dataframes[[i]] <- as.data.frame(data)
  
  if (i == 1) {
    cat("\nStructure du premier data frame :\n")
    str(dataframes[[1]])
    head(dataframes[[1]])
  }
}

cat("\n", number_of_replicates, "data frames générés avec succès.\n")

# Fix: R replaces "-" with "." in column names
# Update envVariableNames to match the actual columns
envVariableNames <- make.names(envVariableNames)
for (i in seq_along(envVariables)) {
  names(envVariables[[i]]) <- envVariableNames[i]
}

# 5. Spatial autocorrelation analysis (for cross-validation) 

# Spatial autocorrelation can bias model performance assessment
# when using standard cross-validation. Therefore, spatial cross-validation
# (Valavi et al. 2019) is used with hexagonal blocks whose size is
# determined by the spatial correlogram analysis.

# Reference: Valavi, R. et al. (2019). "blockCV: An R package for generating
#   spatially or environmentally separated folds for k-fold cross-validation..."

pdf("outputs/Fig4_correlogramme_spatial.pdf", width = 7, height = 5)
par(oma = c(0, 0, 0, 0), mar = c(3, 3, 1, 1))

for (i in 1:number_of_replicates) {
  data <- dataframes[[i]]
  correlogram <- ncf::correlog(
    data[, "longitude"], data[, "latitude"], data[, "response"],
    na.rm = TRUE, increment = 50, resamp = 0, latlon = TRUE
  )
  if (i == 1) {
    plot(correlogram$mean.of.class, correlogram$correlation,
         ann = FALSE, axes = FALSE, lwd = 0.2, cex = 0.5, col = NA,
         ylim = c(-1, 1.03), xlim = c(0, 3000))
    abline(h = 0, lwd = 0.5, col = "red", lty = 2)
    # x-axis
    axis(side = 1, lwd.tick = 0.2, cex.axis = 0.8, lwd = 0.2,  
         at = seq(0,3000, by = 250))
    # y-axis
    axis(side = 2, lwd.tick = 0.2, cex.axis = 0.8, lwd = 0.2)
    title(xlab = "Distance (km)", cex.lab = 1, mgp = c(1.8, 0, 0))
    title(ylab = "Corrélation", cex.lab = 1, mgp = c(1.8, 0, 0))
  }
  lines(correlogram$mean.of.class[-1], correlogram$correlation[-1],
        lwd = 0.2, col = "gray60")
}
dev.off()
cat("Corrélogramme spatial sauvegardé : Figure_spatial_correlogram.pdf\n")

# IMPORTANT: Examine the correlogram to determine the distance from
# which the correlation is no longer positive (crossing zero means that there is no longer spatial correlation). This distance will be used to
# define the size of the spatial cross-validation blocks.
# After inspection, adjust the value below accordingly:

theRanges <- c(1400, 1400) * 1000  # distance in meters (to be adjusted after inspection)
# For example: if correlation is ~0 from 500 km onward, theRanges = c(500,500)*1000
# ---> Reading the spatial correlogram (Fig4):
# The 30 curves show that spatial correlation between observations
# (presence = 1, pseudo-absence = 0) is positive at short distances (~0.6),
# gradually decreases, and drops below zero (red dashed line)
# around 1250-1500 km. Beyond this distance, observations
# are no longer positively autocorrelated.
# A block size of ~1400 km is therefore retained for spatial cross-validation,
# so that the training and test data are
# spatially independent (Valavi et al. 2019).


# 6. BRT model training (Boosted Regression Trees)

# BRT parameters (see Elith et al. 2008 for recommendations)

gbm.x <- envVariableNames          # predictor variables
gbm.y <- "response"                # response variable (presence=1 / pseudo-absence=0)
offset <- NULL                     # no known sampling bias 
tree.complexity <- 5               # number of nodes per tree (interactions between variables)
learning.rate <- 0.001             # contribution of each tree to the model 
bag.fraction <- 0.80               # proportion of observations used at each iteration
n.folds <- 5                       # number of folds for spatial cross-validation
family <- "bernoulli"              # distribution for binary data (presence/absence)
n.trees <- 10                      # initial number of trees (will be incremented by gbm.step)
step.size <- 5                     # tree-number increment step
max.trees <- 10000                 # maximum allowed number of trees
tolerance.method <- "auto"         # method used to decide when to stop the algorithm 
tolerance <- 0.001                 # minimum deviance decrease threshold below which the algorithm stops
prev.stratify <- TRUE              # Could have been omitted because folds are a priori defined by blockCV
plot.main <- TRUE                  # Displays the deviance plot as a function of the number of trees during training (to check whether the model converges)
plot.folds <- FALSE                # No need to generate the folds plot because it was already done with blockCV
verbose <- TRUE                    # gives information about what the model is doing
silent <- FALSE                    # messages are displayed (goes together with verbose)
keep.fold.models <- FALSE          # not necessary here 
keep.fold.vector <- FALSE          # not necessary here 
keep.fold.fit <- FALSE             # not necessary here 

# Lists to store results
brt_models <- list()              
predictions <- list()    
AUCs <- matrix(nrow = number_of_replicates, ncol = 1)
colnames(AUCs) <- "AUC"

cat("\n--- Début de l'entraînement des", number_of_replicates, "modèles BRT ---\n\n")

for (i in 1:number_of_replicates) {
  cat("=== Réplicat", i, "/", number_of_replicates, "===\n")
  
  data <- dataframes[[i]]
  
  # Generate a SpatialPointsDataFrame for spatial cross-validation
  spdf <- SpatialPointsDataFrame(
    data[c("longitude", "latitude")],
    data[, c("response", envVariableNames)],
    proj4string = crs(template)          # use the same projection as the initial template (therefore the same resolution as the environmental variable rasters)
  )
  
  # Generation of hexagonal spatial blocks/folds for cross-validation
  # (use of blockCV to avoid bias due to spatial autocorrelation)
  myblocks <- cv_spatial(
    spdf, column = "response", k = n.folds,
    size = theRanges[1], selection = "random"
  )
  fold.vector <- myblocks$folds_ids
  
  # Creation of weight and constraint vectors
  site.weights <- rep(1, nrow(data))  # here all presences and PAs have the same weight 
  var.monotone <- rep(0, length(gbm.x)) # no constraint imposed on variables
  
  # Reset the number of trees before each analysis
  n.trees <- 10
  
  # BRT model training via gbm.step()
  brt_models[[i]] <- gbm.step(
    data, gbm.x, gbm.y, offset, fold.vector, tree.complexity,
    learning.rate, bag.fraction, site.weights, var.monotone,
    n.folds, prev.stratify, family, n.trees, step.size,
    max.trees, tolerance.method, tolerance, plot.main,
    plot.folds, verbose, silent, keep.fold.models,
    keep.fold.vector, keep.fold.fit
  )
  
  # Retrieve the AUC metric (mean over cross-validation folds)
  AUCs[i, "AUC"] <- brt_models[[i]]$cv.statistics$discrimination.mean
  
  # Prediction of ecological suitability across the whole study area
  dataframe_pred <- as.data.frame(stack(envVariables)) # transforms the 19 rasters into a large table (1 row = 1 raster cell and 1 column = 1 environmental variable)
  dataframe_pred <- dataframe_pred[which(!is.na(rowMeans(dataframe_pred))), ] # removes all cells that are not fully populated
  n.trees_final <- brt_models[[i]]$gbm.call$best.trees # retrieves the optimal number of trees (value at which spatial CV deviance is minimal)
  prediction <- predict.gbm(brt_models[[i]], dataframe_pred, n.trees_final,
                            type = "response", single.tree = FALSE) # runs the BRT model on all cells to obtain their suitability probability
  buffer <- envVariables[[1]] # creates an "empty" raster used as the geographic container for predictions
  buffer[!is.na(buffer[])] <- prediction
  predictions[[i]] <- buffer
  
  cat("  → AUC =", round(AUCs[i, "AUC"], 3), "| Arbres =", n.trees_final, "\n\n")
}

# Save results
write.csv(AUCs, "outputs/AUC_value_replicates.csv", row.names = FALSE, quote = FALSE)
saveRDS(brt_models, "outputs/BRT_model_replicates.rds")
saveRDS(predictions, "outputs/BRT_model_predictions.rds")

cat("\n--- Entraînement terminé ---\n")
cat("AUC moyen :", round(mean(AUCs[, 1]), 2),
    ", AUC range : [", round(min(AUCs), 2), ",", round(max(AUCs), 2), "]\n")

# 7. Risk map generation 

# 7.1 Calculate the mean and standard deviation of predictions across replicates
# Objective: aggregate the 30 prediction rasters into two final maps -
# consensus suitability (mean) and its spatial uncertainty (standard deviation)

# Reloads the list of 30 prediction rasters saved at the end of step 6.
# Allows this section to be resumed without rerunning BRT training.
predictions <- readRDS("outputs/BRT_model_predictions.rds")

# Stacks the 30 rasters into one aligned multi-layer object (RasterStack)
# spatially. Required to apply calc() cell by cell across
# the 30 replicates.
prediction_stack <- stack(predictions)

# Calculates the cell-by-cell mean of the 30 predictions -> risk map
# consensus. na.rm = TRUE ignores NAs (marine cells) to avoid
# propagating missing values over valid terrestrial cells.
mean_prediction <- calc(prediction_stack, fun = mean, na.rm = TRUE)

# Calculates the cell-by-cell standard deviation -> spatial uncertainty map.
# A high standard deviation indicates disagreement between replicates, therefore a prediction
# that is not very robust to stochasticity (BRT + pseudo-absence sampling).
sd_prediction <- calc(prediction_stack, fun = sd, na.rm = TRUE)

# Exports the two rasters as GeoTIFF for external use (QGIS, ArcGIS, etc.)
writeRaster(mean_prediction, "outputs/Risk_map_mean.tif", overwrite = TRUE)
writeRaster(sd_prediction, "outputs/Risk_map_sd.tif", overwrite = TRUE)

# 7.2 Visualize the risk map: mean ecological suitability + coefficient of variation

# Load prediction rasters from the 30 BRT replicates
predictions <- readRDS("outputs/BRT_model_predictions.rds")

# Convert the list of prediction rasters into a RasterStack
prediction_stack <- stack(predictions)

# Compute the mean ecological suitability across BRT replicates
mean_risk_map <- calc(prediction_stack, fun = mean, na.rm = TRUE)

# Compute the standard deviation across BRT replicates
sd_risk_map <- calc(prediction_stack, fun = sd, na.rm = TRUE)

# Compute the coefficient of variation (CV)
# CV = standard deviation / mean.
# It measures relative uncertainty, i.e. uncertainty scaled by the average
# predicted suitability. This is more informative than raw SD when SD values
# are low and poorly contrasted visually.
# The small constant avoids division by zero in cells with near-zero mean suitability.
cv_risk_map <- sd_risk_map / (mean_risk_map + 0.0001)

# Cap the CV scale at the 95th percentile for visual readability.
# This prevents a small number of extreme cells from dominating the colour scale.
cv_values <- getValues(cv_risk_map)
cv_max_coherent <- round(quantile(cv_values, 0.95, na.rm = TRUE), 2)
cv_map_final <- clamp(cv_risk_map, lower = 0, upper = cv_max_coherent)

# Colour palettes
# For mean suitability: blue = low suitability, red = high suitability.
# For CV: same colour convention, where blue = low relative uncertainty
# and red = high relative uncertainty.
cols_risk <- rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100))
cols_cv   <- rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100))

# Export the figure
pdf("outputs/Fig5_cartographie_risque.pdf", width = 10, height = 5)

par(mfrow = c(1, 2),
    oma = c(1, 0.5, 1.5, 0.5),
    mar = c(1, 0.5, 1.5, 3.5),
    lwd = 0.3,
    col = "gray20",
    col.axis = "gray20",
    fg = "gray20",
    family = "serif")

# Panel A: mean ecological suitability

plot(mean_risk_map,
     col = cols_risk,
     zlim = c(0, 1),
     ann = FALSE,
     legend = FALSE,
     axes = FALSE,
     box = FALSE)

plot(coasts, add = TRUE, border = "gray30", lwd = 0.6)

points(occurrences$longitude,
       occurrences$latitude,
       pch = 21,
       cex = 0.45,
       col = "black",
       bg = "red",
       lwd = 0.25)

mtext("A", side = 3, line = 0.2, adj = 0,
      cex = 0.9, font = 2, col = "gray10")

mtext("Suitabilité écologique moyenne",
      side = 1, line = 0.2, cex = 0.78,
      col = "gray30")

plot(mean_risk_map,
     col = cols_risk,
     zlim = c(0, 1),
     legend.only = TRUE,
     add = TRUE,
     legend.width = 1.2,
     legend.shrink = 0.6,
     smallplot = c(0.88, 0.92, 0.15, 0.85),
     axis.args = list(
       cex.axis = 0.65,
       lwd = 0,
       lwd.tick = 0.3,
       col.tick = "gray30",
       tck = -0.5,
       col = "gray30",
       col.axis = "gray30",
       mgp = c(0, 0.4, 0),
       at = seq(0, 1, 0.2)
     ),
     alpha = 1,
     side = 4)

# Panel B: coefficient of variation

plot(cv_map_final,
     col = cols_cv,
     zlim = c(0, cv_max_coherent),
     ann = FALSE,
     legend = FALSE,
     axes = FALSE,
     box = FALSE)

plot(coasts, add = TRUE, border = "gray30", lwd = 0.6)

mtext("B", side = 3, line = 0.2, adj = 0,
      cex = 0.9, font = 2, col = "gray10")

mtext("Coefficient de variation",
      side = 1, line = 0.2, cex = 0.78,
      col = "gray30")

plot(cv_map_final,
     col = cols_cv,
     zlim = c(0, cv_max_coherent),
     legend.only = TRUE,
     add = TRUE,
     legend.width = 1.2,
     legend.shrink = 0.6,
     smallplot = c(0.88, 0.92, 0.15, 0.85),
     axis.args = list(
       cex.axis = 0.65,
       lwd = 0,
       lwd.tick = 0.3,
       col.tick = "gray30",
       tck = -0.5,
       col = "gray30",
       col.axis = "gray30",
       mgp = c(0, 0.4, 0)
     ),
     alpha = 1,
     side = 4)

dev.off()

cat("Cartographie de risque sauvegardée : outputs/Fig5_cartographie_risque.pdf\n")
cat("Échelle CV bornée au 95e percentile :", cv_max_coherent, "\n")


# 8. Performance assessment (AUC + SIppc) 
# Objective: quantify the predictive quality of the 30 BRT models using two
# complementary metrics - AUC (standard ENM) and SIppc (more robust
# to prevalence imbalance and the use of pseudo-absences).

# Reloads the 30 BRT models and the AUC table saved in step 6.
# Allows this section to be resumed without rerunning training.
brt_models <- readRDS("outputs/BRT_model_replicates.rds")
AUCs <- read.csv("outputs/AUC_value_replicates.csv", header = TRUE)

# Summary display of the mean AUC and its [min, max] interval across
# the 30 replicates. The interval reflects variability due to stochasticity
# of the BRT algorithm and pseudo-absence sampling.
cat("\n--- Métriques de performance ---\n")
cat("AUC moyen :", round(mean(AUCs[, 1]), 2),
    ", range : [", round(min(AUCs[, 1]), 2), ",", round(max(AUCs[, 1]), 2), "]\n")

# SIppc calculation
# SIppc = Sorensen's Index Prevalence-Pseudoabsence-Calibrated (Leroy et al. 2018).
# This is an adjusted F1-score measure that accounts for two biases specific to ENM:
#   (1) prevalence is not natural (it is set by the 1:1 ratio that we
#       impose between presences and pseudo-absences)
#   (2) the "absences" are actually pseudo-absences, potentially
#       false negatives (the virus could circulate there without being detected)
# SIppc corrects these biases through an x factor (see below) that weights
# false positives from pseudo-absences. It is also independent of the
# chosen classification threshold: the threshold that maximizes it is found for each model.

# Pre-allocation of vectors that will store, for each replicate, the maximum SIppc
# and the classification threshold (0 to 1) that produces it.
SI_ppcs <- rep(NA, length(brt_models))
thresholds <- rep(NA, length(brt_models))

# Loop over the 30 BRT models: SIppc calculation for each one.
for (i in 1:length(brt_models)) {
  tmp <- matrix(nrow = 101, ncol = 2)   # Temporary matrix that stores, for each tested threshold (0, 0.01, 0.02, ..., 1),
  # the corresponding SIppc value. Used to identify the optimal threshold.
  colnames(tmp) <- c("threshold", "SIppc")
  tmp[, "threshold"] <- seq(0, 1, 0.01)
  
  # Retrieves the model training dataframe (300 rows: 150 presences
  # + 150 pseudo-absences) and separates the response column from the 19
  # environmental variables. Note: evaluation is done here on training data
  dataframe <- brt_models[[i]]$gbm.call$dataframe
  responses <- dataframe$response
  data_pred <- dataframe[, 4:ncol(dataframe)]
  
  # Uses the optimal number of trees (best iteration) selected by gbm.step.
  n.trees_final <- brt_models[[i]]$gbm.call$best.trees
  
  # Predicts the suitability probability for each observation (between 0 and 1).
  prediction <- predict.gbm(brt_models[[i]], data_pred, n.trees_final,
                            type = "response", single.tree = FALSE)
  
  # Counts used by the SIppc formula:
  P <- sum(responses == 1)  # number of presences
  A <- sum(responses == 0)  # number of pseudo-absences
  prev <- P / (P + A)       # prev = observed prevalence (proportion of presences in the dataset)
  x <- (P / A) * ((1 - prev) / prev)
  SI_ppc <- 0 # Initialize the maximum SIppc at 0; it will be updated in the loop.
  
  # Loop over 101 possible classification thresholds (0, 0.01, ..., 1).
  # For each threshold, calculate TP/FN/FP and evaluate SIppc.
  for (threshold in seq(0, 1, 0.01)) {
    TP <- length(which((responses == 1) & (prediction >= threshold))) # TP: true presences correctly predicted (response = 1 and pred >= threshold)
    FN <- length(which((responses == 1) & (prediction < threshold)))  # FN: missed presences (response = 1 but pred < threshold)
    FP_pa <- length(which((responses == 0) & (prediction >= threshold)))    # FP_pa: pseudo-absences incorrectly predicted as presences
    # (response = 0 but pred >= threshold). Note: "FP_pa" emphasizes that these FP
    
    # SIppc formula (Leroy et al. 2018):
    # SIppc = 2×TP / (2×TP + x×FP_pa + FN)
    # This is a modified F1-score where x weights FP from pseudo-absences.
    # The closer SIppc is to 1, the better the sensitivity/precision trade-off.                                                                        # come from pseudo-absences, not confirmed absences.
    SI_ppc_tmp <- (2 * TP) / ((2 * TP) + (x * FP_pa) + FN)
    tmp[which(tmp[, "threshold"] == threshold), "SIppc"] <- SI_ppc_tmp
    
    # Keeps track of the threshold that maximizes SIppc for this replicate.
    # This is the threshold that would be used to binarize the risk map
    # into risk areas / no-risk areas.
    if (SI_ppc < SI_ppc_tmp) {
      SI_ppc <- SI_ppc_tmp
      optimised_threshold <- threshold
    }
  }
  
  # Stores the maximum SIppc and its optimal threshold for this replicate
  SI_ppcs[i] <- SI_ppc
  thresholds[i] <- optimised_threshold
}

# Summary display: mean SIppc and range across the 30 replicates, plus the
# mean threshold that maximizes SIppc. A threshold around 0.5 suggests a well-calibrated model;
# a very low or very high threshold would indicate an imbalance between sensitivity
# (ability to detect presences) and specificity (ability to exclude
# pseudo-absences).
cat("SIppc moyen :", round(mean(SI_ppcs), 2),
    ", range : [", round(min(SI_ppcs), 2), ",", round(max(SI_ppcs), 2), "]\n")
cat("Seuil optimal moyen :", round(mean(thresholds), 2), "\n")

# 9. Relative importance of environmental variables 

brt_models <- readRDS("outputs/BRT_model_replicates.rds")

relative_importances <- matrix(nrow = length(brt_models), ncol = length(envVariables))
colnames(relative_importances) <- envVariableNames

for (i in 1:length(brt_models)) {
  for (j in 1:length(envVariables)) {
    relative_importances[i, j] <- summary(brt_models[[i]])[envVariableNames[j], "rel.inf"]
  }
}

write.csv(relative_importances, "outputs/Relative_importances.csv", row.names = FALSE, quote = FALSE)

# Display a summary of relative importances (mean + 95% CI)
cat("\n--- Importance relative des variables ---\n")
RI_summary <- matrix(nrow = length(envVariables), ncol = 1)
rownames(RI_summary) <- envVariableNames
colnames(RI_summary) <- "mean RI [95% CI]"

for (i in 1:ncol(relative_importances)) {
  mean_RI <- round(mean(relative_importances[, i]), 1)
  ci95_RI <- round(t.test(relative_importances[, i])$conf.int[1:2], 1)
  RI_summary[i, 1] <- paste0(mean_RI, " [", ci95_RI[1], ", ", ci95_RI[2], "]")
}
print(RI_summary)

# 10. Response curves

# Response curves show how ecological suitability varies as a function
# of one environmental variable, while the others are fixed at their median value.

brt_models <- readRDS("outputs/BRT_model_replicates.rds")

# Calculate the median, min and max values of each variable at presence sites
sp_pres <- SpatialPoints(dataframes[[1]][which(dataframes[[1]][, "response"] == 1), 2:3])
envVariableValues <- matrix(nrow = 3, ncol = length(envVariables))
colnames(envVariableValues) <- envVariableNames
rownames(envVariableValues) <- c("median", "minV", "maxV")

for (i in 1:length(envVariables)) {
  pts <- rasterize(sp_pres, envVariables[[i]])
  rast <- envVariables[[i]]
  rast[is.na(pts)] <- NA
  envVariableValues[, i] <- c(
    median(rast[], na.rm = TRUE),
    min(rast[], na.rm = TRUE),
    max(rast[], na.rm = TRUE)
  )
}

# Generate response curves for a selection of variables
# (the most influential ones - to be adjusted after RI inspection)
vars_for_curves <- 1:length(envVariableNames)  # all variables

pdf("outputs/Fig6_courbes_reponse.pdf", width = 20, height = 8)
par(mfrow = c(2, ceiling(length(vars_for_curves) / 2)),
    oma = c(0, 0, 2, 0), mar = c(2.5, 2, 1, 0.5), lwd = 0.2, col = "gray30")

for (i in vars_for_curves) {
  valuesInterval <- (envVariableValues["maxV", i] - envVariableValues["minV", i]) / 100
  if (valuesInterval == 0) next  # skip if there is no variation
  
  dataframe_rc <- data.frame(matrix(
    nrow = length(seq(envVariableValues["minV", i], envVariableValues["maxV", i], valuesInterval)),
    ncol = length(envVariables)
  ))
  colnames(dataframe_rc) <- envVariableNames
  
  for (j in 1:length(envVariables)) {
    interval_j <- (envVariableValues["maxV", j] - envVariableValues["minV", j]) / 100
    if (i == j) {
      dataframe_rc[, envVariableNames[j]] <- seq(
        envVariableValues["minV", j], envVariableValues["maxV", j], interval_j
      )[1:nrow(dataframe_rc)]
    } else {
      dataframe_rc[, envVariableNames[j]] <- rep(envVariableValues["median", j], nrow(dataframe_rc))
    }
  }
  
  preds_rc <- list()
  for (j in 1:length(brt_models)) {
    n.trees_final <- brt_models[[j]]$gbm.call$best.trees
    preds_rc[[j]] <- predict.gbm(brt_models[[j]], newdata = dataframe_rc,
                                 n.trees_final, type = "response", single.tree = FALSE)
  }
  
  col_line <- rgb(222, 67, 39, 255, maxColorValue = 255)
  for (j in 1:length(brt_models)) {
    if (j == 1) {
      plot(dataframe_rc[, envVariableNames[i]], preds_rc[[j]], col = col_line,
           ann = FALSE, axes = FALSE, lwd = 0.3, type = "l", ylim = c(0, 1))
    } else {
      lines(dataframe_rc[, envVariableNames[i]], preds_rc[[j]], col = col_line, lwd = 0.3)
    }
  }
  axis(side = 1, lwd.tick = 0.2, cex.axis = 0.65, lwd = 0, tck = -0.04,
       col.axis = "gray30", mgp = c(0, 0.15, 0))
  axis(side = 2, lwd.tick = 0.2, cex.axis = 0.65, lwd = 0, tck = -0.04,
       col.axis = "gray30", mgp = c(0, 0.4, 0))
  title(xlab = gsub("_", " ", envVariableNames[i]), cex.lab = 0.7,
        mgp = c(1.2, 0, 0), col.lab = "gray30")
  box(lwd = 0.2, col = "gray30")
}

mtext("Courbes de réponse - Variables environnementales",
      side = 3, outer = TRUE, cex = 1.2, font = 2)
dev.off()
cat("Courbes de réponse sauvegardées : Figure_response_curves.pdf\n")

# 11. Final summary
cat("\n========================================\n")
cat("  RÉSUMÉ DES RÉSULTATS - PROJET 1\n")
cat("========================================\n")
cat("Nombre de réplicats BRT :", number_of_replicates, "\n")
cat("AUC moyen :", round(mean(AUCs[, 1]), 3), "\n")
cat("SIppc moyen :", round(mean(SI_ppcs), 3), "\n")
cat("Taille des blocs spatiaux :", theRanges[1] / 1000, "km\n")
cat("Fichiers de sortie :\n")
cat("  - Figure_exploration_variables.pdf\n")
cat("  - Figure_spatial_correlogram.pdf\n")
cat("  - Figure_risk_map.pdf\n")
cat("  - Figure_response_curves.pdf\n")
cat("  - Risk_map_mean.tif / Risk_map_sd.tif\n")
cat("  - AUC_value_replicates.csv\n")
cat("  - Relative_importances.csv\n")
cat("  - BRT_model_replicates.rds / BRT_model_predictions.rds\n")
cat("========================================\n")

# END OF SCRIPT

cat("END OF PRINCIPAL SCRIPT !")




#===============================================================================

 #PHASE OF DIAGNOSTIC : Why is RH winter's RI = 59,4% ? 
diagnostic_dir <- "/Users/hannahthiels/Desktop/BING-F432_Travail_Groupe/Projet1_Risk_Mapping/outputs/Diagnostic des résultats (notamment RI du HR hivernal)"

if (!dir.exists(diagnostic_dir)) {
  dir.create(diagnostic_dir, recursive = TRUE)
}

# 1) Does the winter HR system heavily discriminate against actual absences and pseudo-absences on its own?

#The univariate AUC analysis showed that relative_humidity_winter was the strongest single predictor
#for separating presences from pseudo-absences, with a mean AUC of 0.855 across the 30 replicates. 
#However, several temperature predictors also showed strong univariate discrimination, 
#notably temperature_inFall (mean AUC = 0.817), temperature_winter (mean AUC = 0.816), temperature_spring (mean AUC = 0.807), 
#and relative_humidity_inFall (mean AUC = 0.801). Therefore, the low relative importance assigned to 
#temperature variables in the final BRT models should not be interpreted as evidence that temperature carries little ecological information. 
#Rather, it suggests that multiple climatic predictors contain overlapping discriminatory information, 
#and that the BRT may preferentially allocate split-based importance to relative_humidity_winter.


#PS : For each predictor, a univariate AUC was computed separately for each of the 30 pseudo-absence replicates. 
#Because the objective was to quantify discriminatory power regardless of direction, absolute AUC was used, defined as max(AUC, 1 - AUC).
#Variables marked with an asterisk had a mean raw AUC below 0.5, meaning that they separated presences from pseudo-absences in the inverse direction. 
#Bars represent mean absolute AUC across replicates, and error bars represent ±1 standard deviation.

univ_auc_raw <- matrix(NA, nrow = length(dataframes), ncol = length(envVariableNames))
univ_auc_abs <- matrix(NA, nrow = length(dataframes), ncol = length(envVariableNames))

colnames(univ_auc_raw) <- envVariableNames
colnames(univ_auc_abs) <- envVariableNames

for (i in seq_along(dataframes)) {
  d <- dataframes[[i]]
  
  for (v in envVariableNames) {
    presence_values <- d[d$response == 1, v]
    pseudoabsence_values <- d[d$response == 0, v]
    
    auc_raw <- dismo::evaluate(
      p = presence_values,
      a = pseudoabsence_values
    )@auc
    
    univ_auc_raw[i, v] <- auc_raw
    univ_auc_abs[i, v] <- max(auc_raw, 1 - auc_raw)
  }
}

univ_auc_summary <- data.frame(
  variable = envVariableNames,
  mean_AUC_raw = apply(univ_auc_raw, 2, mean, na.rm = TRUE),
  sd_AUC_raw = apply(univ_auc_raw, 2, sd, na.rm = TRUE),
  mean_AUC_abs = apply(univ_auc_abs, 2, mean, na.rm = TRUE),
  sd_AUC_abs = apply(univ_auc_abs, 2, sd, na.rm = TRUE),
  min_AUC_abs = apply(univ_auc_abs, 2, min, na.rm = TRUE),
  max_AUC_abs = apply(univ_auc_abs, 2, max, na.rm = TRUE)
)

univ_auc_summary$direction <- ifelse(
  univ_auc_summary$mean_AUC_raw >= 0.5,
  "positive",
  "inverse"
)

univ_auc_summary <- univ_auc_summary[order(-univ_auc_summary$mean_AUC_abs), ]

print(univ_auc_summary)
write.csv(univ_auc_summary,
          file.path(diagnostic_dir, "Diagnostic_univariate_AUCs_summary.csv"),
          row.names = FALSE,
          quote = FALSE)



#PLOT STEP 1 OF DIAGNOSTIC - Mean absolute univariate AUC with SD and inverse-direction marks

auc_plot <- univ_auc_summary

auc_plot$label_abs <- sprintf("%.3f", auc_plot$mean_AUC_abs)
auc_plot$star <- ifelse(auc_plot$direction == "inverse", "*", "")

pdf(file.path(diagnostic_dir, "Diagnostic_univariate_AUCs_clean_with_SD.pdf"),
    width = 10, height = 7)

par(mar = c(5, 10, 4, 4),
    lwd = 0.3,
    col = "gray30",
    col.axis = "gray30",
    fg = "gray30")

bar_cols <- ifelse(
  auc_plot$variable == "relative_humidity_winter",
  "#D95F02",
  ifelse(auc_plot$mean_AUC_abs >= 0.80, "#7570B3", "gray75")
)

bp <- barplot(
  rev(auc_plot$mean_AUC_abs),
  names.arg = rev(auc_plot$variable),
  horiz = TRUE,
  las = 1,
  cex.names = 0.72,
  xlim = c(0.5, 0.93),
  col = rev(bar_cols),
  border = NA,
  xlab = "Mean absolute univariate AUC across 30 replicates",
  main = "Univariate discriminative power of predictors"
)

# Thin SD error bars
arrows(
  x0 = rev(auc_plot$mean_AUC_abs - auc_plot$sd_AUC_abs),
  y0 = bp,
  x1 = rev(auc_plot$mean_AUC_abs + auc_plot$sd_AUC_abs),
  y1 = bp,
  angle = 90,
  code = 3,
  length = 0.025,
  col = "gray25",
  lwd = 0.7
)

# Reference lines
abline(v = 0.5, lty = 2, col = "gray70")
abline(v = 0.8, lty = 3, col = "gray70")

# Absolute AUC labels
text(
  x = rev(auc_plot$mean_AUC_abs) + 0.012,
  y = bp,
  labels = rev(auc_plot$label_abs),
  cex = 0.78,
  font = 2,
  col = "gray15",
  xpd = TRUE,
  adj = 0
)

# Visible stars for inverse raw AUC
text(
  x = rev(auc_plot$mean_AUC_abs) + 0.047,
  y = bp,
  labels = rev(auc_plot$star),
  cex = 1.25,
  font = 2,
  col = "#B2182B",
  xpd = TRUE,
  adj = 0
)

legend("bottomright",
       legend = c("Winter relative humidity",
                  "Other predictors with absolute AUC ≥ 0.80",
                  "Other predictors",
                  "Thin bars = ±1 SD"),
       fill = c("#D95F02", "#7570B3", "gray75", NA),
       border = NA,
       bty = "n",
       cex = 0.72)

mtext("* mean raw AUC < 0.5; absolute AUC = 1 - raw AUC",
      side = 1,
      line = 3.5,
      cex = 0.75,
      col = "#B2182B")

dev.off()



# 2) Effect of removing winter relative humidity from the BRT

# Objective:
#   The univariate AUC analysis showed that relative_humidity_winter has the
#   strongest standalone discriminatory power. This step tests whether the BRT
#   model truly depends on this variable, or whether other climatic predictors
#   can compensate for its removal.
#
#   For each replicate, two BRT models are fitted using the same spatial folds:
#     (1) a full model with all predictors,
#     (2) a reduced model excluding relative_humidity_winter.
#
#   We then compare:
#     - cross-validated AUC,
#     - relative importance redistribution among predictors.

var_removed <- "relative_humidity_winter"

gbm.x_full <- envVariableNames
gbm.x_reduced <- setdiff(envVariableNames, var_removed)

diagnostic_replicates <- 1:number_of_replicates

auc_comparison <- data.frame(
  replicate = diagnostic_replicates,
  AUC_full_refit = NA_real_,
  AUC_without_RHwinter = NA_real_,
  delta_AUC = NA_real_
)

RI_full_diag <- matrix(NA,
                       nrow = length(diagnostic_replicates),
                       ncol = length(gbm.x_full))
colnames(RI_full_diag) <- gbm.x_full

RI_reduced_diag <- matrix(NA,
                          nrow = length(diagnostic_replicates),
                          ncol = length(gbm.x_reduced))
colnames(RI_reduced_diag) <- gbm.x_reduced

extract_RI <- function(model) {
  ri <- summary(model, plotit = FALSE)
  out <- setNames(ri$rel.inf, ri$var)
  return(out)
}

for (idx in seq_along(diagnostic_replicates)) {
  
  i <- diagnostic_replicates[idx]
  cat("\n=== Diagnostic replicate", i, "/", length(diagnostic_replicates), "===\n")
  
  d <- dataframes[[i]]
  
  spdf_diag <- SpatialPointsDataFrame(
    d[c("longitude", "latitude")],
    d[, c("response", gbm.x_full)],
    proj4string = crs(template)
  )
  
  myblocks_diag <- cv_spatial(
    spdf_diag,
    column = "response",
    k = n.folds,
    size = theRanges[1],
    selection = "random"
  )
  
  fold.vector_diag <- myblocks_diag$folds_ids
  
  site.weights.full <- rep(1, nrow(d))
  var.monotone.full <- rep(0, length(gbm.x_full))
  
  site.weights.reduced <- rep(1, nrow(d))
  var.monotone.reduced <- rep(0, length(gbm.x_reduced))
  
  cat("Fitting full model...\n")
  brt_full_diag <- gbm.step(
    d, gbm.x_full, gbm.y, offset, fold.vector_diag,
    tree.complexity, learning.rate, bag.fraction,
    site.weights.full, var.monotone.full,
    n.folds, prev.stratify, family,
    n.trees = 10, step.size, max.trees,
    tolerance.method, tolerance,
    plot.main = FALSE, plot.folds = FALSE,
    verbose = FALSE, silent = TRUE,
    keep.fold.models = FALSE,
    keep.fold.vector = FALSE,
    keep.fold.fit = FALSE
  )
  
  cat("Fitting reduced model without relative_humidity_winter...\n")
  brt_reduced_diag <- gbm.step(
    d, gbm.x_reduced, gbm.y, offset, fold.vector_diag,
    tree.complexity, learning.rate, bag.fraction,
    site.weights.reduced, var.monotone.reduced,
    n.folds, prev.stratify, family,
    n.trees = 10, step.size, max.trees,
    tolerance.method, tolerance,
    plot.main = FALSE, plot.folds = FALSE,
    verbose = FALSE, silent = TRUE,
    keep.fold.models = FALSE,
    keep.fold.vector = FALSE,
    keep.fold.fit = FALSE
  )
  
  auc_full <- brt_full_diag$cv.statistics$discrimination.mean
  auc_reduced <- brt_reduced_diag$cv.statistics$discrimination.mean
  
  auc_comparison$AUC_full_refit[idx] <- auc_full
  auc_comparison$AUC_without_RHwinter[idx] <- auc_reduced
  auc_comparison$delta_AUC[idx] <- auc_full - auc_reduced
  
  ri_full <- extract_RI(brt_full_diag)
  RI_full_diag[idx, names(ri_full)] <- ri_full
  
  ri_reduced <- extract_RI(brt_reduced_diag)
  RI_reduced_diag[idx, names(ri_reduced)] <- ri_reduced
  
  cat("AUC full:", round(auc_full, 3),
      "| AUC without RH_winter:", round(auc_reduced, 3),
      "| delta:", round(auc_full - auc_reduced, 3), "\n")
}

write.csv(auc_comparison,
          file.path(diagnostic_dir, "Diagnostic_AUC_with_vs_without_RHwinter.csv"),
          row.names = FALSE,
          quote = FALSE)

write.csv(RI_full_diag,
          file.path(diagnostic_dir, "Diagnostic_RI_full_refit.csv"),
          row.names = FALSE,
          quote = FALSE)

write.csv(RI_reduced_diag,
          file.path(diagnostic_dir, "Diagnostic_RI_without_RHwinter.csv"),
          row.names = FALSE,
          quote = FALSE)

cat("\n--- Diagnostic AUC comparison ---\n")
print(summary(auc_comparison))

cat("\nMean AUC full refit:",
    round(mean(auc_comparison$AUC_full_refit, na.rm = TRUE), 3), "\n")
cat("Mean AUC without RH_winter:",
    round(mean(auc_comparison$AUC_without_RHwinter, na.rm = TRUE), 3), "\n")
cat("Mean delta AUC:",
    round(mean(auc_comparison$delta_AUC, na.rm = TRUE), 3), "\n")

cat("\nDiagnostic outputs saved in:\n")
cat(diagnostic_dir, "\n")


# PLOTS STEP 2) : effect of removing relative_humidity_winter

#   Fig 1  Paired dot plot of cross-validated AUC (full vs reduced)
#   Fig 2  Sorted dot plot of delta_AUC per replicate, mean +/- SD
#   Fig 3  Dumbbell plot of mean relative-importance redistribution
#          for a curated subset of variables

auc <- read.csv(file.path(diagnostic_dir,
                          "Diagnostic_AUC_with_vs_without_RHwinter.csv"))

RI_full <- read.csv(file.path(diagnostic_dir,
                              "Diagnostic_RI_full_refit.csv"),
                    check.names = FALSE)

RI_red  <- read.csv(file.path(diagnostic_dir,
                              "Diagnostic_RI_without_RHwinter.csv"),
                    check.names = FALSE)

# FIGURE 1 - Paired dot plot of AUC (full vs reduced)


pdf(file.path(diagnostic_dir, "diagnostic_paired_AUC_full_vs_reduced.pdf"),
    width = 6.5, height = 5.2)
par(mar = c(4.5, 4.5, 3, 1.5))

x_pos <- c(1, 2)
ylim  <- range(c(auc$AUC_full_refit, auc$AUC_without_RHwinter),
               na.rm = TRUE) + c(-0.015, 0.030)

plot(NA, xlim = c(0.55, 2.45), ylim = ylim,
     xaxt = "n", xlab = "", ylab = "Cross-validated AUC",
     main = "Paired AUC: full model vs model without relative_humidity_winter",
     bty = "l", cex.main = 0.95, las = 1)
axis(1, at = x_pos, labels = c("Full model", "Without RH_winter"))

segments(x0 = rep(x_pos[1], nrow(auc)),
         y0 = auc$AUC_full_refit,
         x1 = rep(x_pos[2], nrow(auc)),
         y1 = auc$AUC_without_RHwinter,
         col = adjustcolor("grey55", alpha.f = 0.55), lwd = 0.8)

set.seed(1)
points(jitter(rep(x_pos[1], nrow(auc)), amount = 0.06),
       auc$AUC_full_refit,
       pch = 21, bg = adjustcolor("steelblue", alpha.f = 0.7),
       col = "steelblue4", cex = 1.05)
points(jitter(rep(x_pos[2], nrow(auc)), amount = 0.06),
       auc$AUC_without_RHwinter,
       pch = 21, bg = adjustcolor("firebrick", alpha.f = 0.7),
       col = "firebrick4", cex = 1.05)

mean_full <- mean(auc$AUC_full_refit,      na.rm = TRUE)
mean_red  <- mean(auc$AUC_without_RHwinter, na.rm = TRUE)
segments(x_pos[1] - 0.18, mean_full, x_pos[1] + 0.18, mean_full,
         col = "steelblue4", lwd = 3)
segments(x_pos[2] - 0.18, mean_red,  x_pos[2] + 0.18, mean_red,
         col = "firebrick4", lwd = 3)
text(x_pos[1], ylim[2] - 0.005,
     sprintf("Mean = %.3f", mean_full),
     col = "steelblue4", cex = 0.85, font = 2)
text(x_pos[2], ylim[2] - 0.005,
     sprintf("Mean = %.3f", mean_red),
     col = "firebrick4", cex = 0.85, font = 2)

dev.off()

# FIGURE 2 - Horizontal strip plot of delta_AUC

# Each blue point = one replicate (vertical jitter only to avoid
# overplotting). The firebrick bar above the dots shows the mean
# +/- 1 SD across replicates. The dotted grey vertical line at 0
# is the no-effect reference.

pdf(file.path(diagnostic_dir, "diagnostic_delta_AUC_distribution.pdf"),
    width = 6.5, height = 3.2)
par(mar = c(4.5, 1.5, 3.5, 1.5))

n_rep      <- nrow(auc)
mean_delta <- mean(auc$delta_AUC, na.rm = TRUE)
sd_delta   <- sd(auc$delta_AUC,   na.rm = TRUE)

data_rng <- range(c(0, auc$delta_AUC,
                    mean_delta - sd_delta, mean_delta + sd_delta),
                  na.rm = TRUE)
xlim     <- data_rng + c(-1, 1) * 0.08 * diff(data_rng)

plot(NA, xlim = xlim, ylim = c(0, 1.8),
     xlab = expression(Delta * "AUC = AUC"["full"] - "AUC"["without RH_winter"]),
     ylab = "", yaxt = "n", bty = "n",
     main = "Distribution of AUC loss across replicates",
     cex.main = 0.95, las = 1)

# Zero reference
abline(v = 0, col = "grey60", lty = 3)

# Individual replicate points with vertical jitter
set.seed(1)
y_pts <- 0.55 + runif(n_rep, -0.22, 0.22)
points(auc$delta_AUC, y_pts,
       pch = 21, bg = adjustcolor("steelblue", alpha.f = 0.75),
       col = "steelblue4", cex = 1.2)

# Mean +/- SD bar above the dots
y_sum <- 1.30
arrows(mean_delta - sd_delta, y_sum,
       mean_delta + sd_delta, y_sum,
       code = 3, angle = 90, length = 0.04,
       col = "firebrick", lwd = 2.2)
points(mean_delta, y_sum, pch = 18, col = "firebrick", cex = 2.4)

# Numerical annotation
text(mean_delta, 1.65,
     sprintf("Mean AUC loss = %.3f \u00b1 %.3f", mean_delta, sd_delta),
     col = "firebrick", cex = 0.95, font = 2)

# Compact legend
legend("topleft",
       legend = c(sprintf("Replicates (n = %d)", n_rep), "Mean \u00b1 SD"),
       pch = c(21, 18),
       pt.bg = c("steelblue", NA),
       col = c("steelblue4", "firebrick"),
       pt.cex = c(1.1, 1.8),
       bty = "n", cex = 0.78, inset = c(0, -0.02))

dev.off()

# FIGURE 3 Dumbbell plot of RI redistribution


make_RI_dumbbell <- function(var_list,
                             output_pdf,
                             main_title,
                             pdf_width = 9,
                             pdf_height = 5.5,
                             label_cex = 0.85,
                             value_cex = 0.72) {
  
  # Mean RI in each model, restricted to var_list
  vars_full <- var_list[var_list %in% colnames(RI_full)]
  mean_RI_full <- colMeans(RI_full[, vars_full, drop = FALSE], na.rm = TRUE)
  
  vars_red <- var_list[var_list %in% colnames(RI_red)]
  mean_RI_red <- colMeans(RI_red[, vars_red, drop = FALSE], na.rm = TRUE)
  
  comp <- data.frame(
    variable = var_list,
    RI_full  = as.numeric(mean_RI_full[var_list]),
    RI_red   = NA_real_,
    stringsAsFactors = FALSE
  )
  comp$RI_red[comp$variable %in% vars_red] <-
    as.numeric(mean_RI_red[comp$variable[comp$variable %in% vars_red]])
  comp$delta_RI <- comp$RI_red - comp$RI_full
  
  # Order by full-model RI (largest at top)
  comp <- comp[order(-comp$RI_full), ]
  
  pdf(output_pdf, width = pdf_width, height = pdf_height)
  par(mar = c(4.5, 13, 3, 2))
  
  n    <- nrow(comp)
  xmax <- max(c(comp$RI_full, comp$RI_red), na.rm = TRUE) * 1.20
  
  plot(NA, xlim = c(0, xmax), ylim = c(0.5, n + 0.5),
       yaxt = "n", xlab = "Mean relative importance (%)", ylab = "",
       main = main_title,
       bty = "l", cex.main = 0.95, las = 1)
  axis(2, at = n:1, labels = comp$variable,
       las = 2, cex.axis = label_cex)
  abline(v = pretty(c(0, xmax)), col = "grey94", lty = 1)
  
  for (i in seq_len(n)) {
    y  <- n - i + 1
    rf <- comp$RI_full[i]
    rr <- comp$RI_red[i]
    
    if (!is.na(rr)) {
      col_seg <- if (rr > rf) "forestgreen" else "firebrick"
      segments(rf, y, rr, y, col = col_seg, lwd = 2.5)
    }
    
    points(rf, y, pch = 21, bg = "steelblue",
           col = "steelblue4", cex = 1.6)
    
    if (!is.na(rr)) {
      points(rr, y, pch = 21, bg = "firebrick",
             col = "firebrick4", cex = 1.6)
      text(rf, y, sprintf("%.1f", rf),
           pos = 1, cex = value_cex, col = "steelblue4", offset = 0.5)
      text(rr, y, sprintf("%.1f", rr),
           pos = 3, cex = value_cex, col = "firebrick4", offset = 0.5)
    } else {
      text(rf, y, sprintf("%.1f  (removed)", rf),
           pos = 4, cex = value_cex + 0.05, col = "steelblue4",
           offset = 0.6, font = 3)
    }
  }
  
  legend("bottomright",
         legend = c("Full model RI",
                    "Reduced model RI",
                    "Gain after removal",
                    "Loss after removal"),
         pch   = c(21, 21, NA, NA),
         pt.bg = c("steelblue", "firebrick", NA, NA),
         col   = c("steelblue4", "firebrick4", "forestgreen", "firebrick"),
         lwd   = c(NA, NA, 2.5, 2.5),
         lty   = c(NA, NA, 1, 1),
         bty = "n", cex = 0.78)
  
  dev.off()
  
  return(comp)
}

#  plot 3: All environmental variables for RI redistribution after removing relative_humidity_winter 

ri_means_all <- colMeans(RI_full, na.rm = TRUE)
ranked_vars <- names(sort(ri_means_all, decreasing = TRUE))
all_vars <- ranked_vars  # already ordered by RI_full

comp_all <- make_RI_dumbbell(
  var_list   = all_vars,
  output_pdf = file.path(diagnostic_dir,
                         "diagnostic_RI_redistribution_all_variables.pdf"),
  main_title = "All predictors: RI redistribution after removing relative_humidity_winter",
  pdf_width  = 9,
  pdf_height = max(5.5, 0.38 * length(all_vars) + 2),
  label_cex  = 0.75,
  value_cex  = 0.65
)

write.csv(comp_all,
          file.path(diagnostic_dir,
                    "Diagnostic_RI_all_variables_summary.csv"),
          row.names = FALSE)


cat("All diagnostic figures saved in:\n", diagnostic_dir, "\n")


# 3) Correlation structure among environmental predictors

# Objective:
#   The previous diagnostic steps suggested that several climatic predictors
#   contain overlapping discriminatory information. This step quantifies the
#   correlation structure among all environmental predictors to assess whether
#   the high RI of relative_humidity_winter may be partly related to redundancy
#   among climatic variables.

#   Fig   Spearman correlation heatmap between all environmental predictors
#          to assess redundancy among climatic variables and support the
#          interpretation of RI redistribution after removing relative_humidity_winter

# PS : Spearman correlation is used because relationships between environmental
#   predictors may be monotonic but not strictly linear.

diagnostic_dir <- "/Users/hannahthiels/Desktop/BING-F432_Travail_Groupe/Projet1_Risk_Mapping/outputs/Diagnostic des résultats (notamment RI du HR hivernal)"

if (!dir.exists(diagnostic_dir)) {
  dir.create(diagnostic_dir, recursive = TRUE)
}

# Use all terrestrial raster cells with complete environmental information
env_data <- as.data.frame(stack(envVariables))
env_data <- env_data[complete.cases(env_data), ]

# Keep only environmental variables used in the BRT
env_data <- env_data[, envVariableNames]

# Spearman correlation matrix
cor_matrix <- cor(env_data,
                  method = "spearman",
                  use = "complete.obs")

write.csv(round(cor_matrix, 3),
          file.path(diagnostic_dir, "Diagnostic_spearman_correlation_matrix.csv"),
          row.names = TRUE,
          quote = FALSE)

# Extract correlations with relative_humidity_winter
cor_RHwinter <- sort(cor_matrix["relative_humidity_winter", ],
                     decreasing = TRUE)

write.csv(data.frame(
  variable = names(cor_RHwinter),
  spearman_rho = as.numeric(cor_RHwinter)
),
file.path(diagnostic_dir, "Diagnostic_correlations_with_RHwinter.csv"),
row.names = FALSE,
quote = FALSE)

cat("\n--- Spearman correlations with relative_humidity_winter ---\n")
print(round(cor_RHwinter, 3))


# Correlation heatmap

pdf(file.path(diagnostic_dir, "diagnostic_spearman_correlation_heatmap.pdf"),
    width = 9, height = 9)

par(mar = c(10, 10, 4, 2),
    lwd = 0.3,
    col = "gray30",
    col.axis = "gray30",
    fg = "gray30")

# Reorder variables by hierarchical clustering for readability
dist_cor <- as.dist(1 - abs(cor_matrix))
hc <- hclust(dist_cor, method = "average")
ord <- hc$order

cor_ord <- cor_matrix[ord, ord]
var_ord <- colnames(cor_ord)

# Color palette: blue = negative, white = no correlation, red = positive
cols <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(200)

image(
  x = 1:ncol(cor_ord),
  y = 1:nrow(cor_ord),
  z = t(cor_ord[nrow(cor_ord):1, ]),
  col = cols,
  breaks = seq(-1, 1, length.out = 201),
  axes = FALSE,
  xlab = "",
  ylab = "",
  main = "Spearman correlation between environmental predictors"
)

axis(1,
     at = 1:ncol(cor_ord),
     labels = var_ord,
     las = 2,
     cex.axis = 0.65)

axis(2,
     at = 1:nrow(cor_ord),
     labels = rev(var_ord),
     las = 2,
     cex.axis = 0.65)

box()

# Add correlation values inside cells
for (i in 1:nrow(cor_ord)) {
  for (j in 1:ncol(cor_ord)) {
    val <- cor_ord[i, j]
    text(j,
         nrow(cor_ord) - i + 1,
         labels = sprintf("%.2f", val),
         cex = 0.45,
         col = ifelse(abs(val) > 0.65, "white", "gray20"))
  }
}

# Color legend
legend_x <- ncol(cor_ord) + 1.2
legend_y <- seq(1, nrow(cor_ord), length.out = 200)

par(xpd = TRUE)

rect_x_left <- ncol(cor_ord) + 0.7
rect_x_right <- ncol(cor_ord) + 1.0

for (k in 1:199) {
  rect(rect_x_left,
       legend_y[k],
       rect_x_right,
       legend_y[k + 1],
       col = cols[k],
       border = NA)
}

text(ncol(cor_ord) + 1.35, 1, "-1", cex = 0.7)
text(ncol(cor_ord) + 1.35, nrow(cor_ord) / 2, "0", cex = 0.7)
text(ncol(cor_ord) + 1.35, nrow(cor_ord), "1", cex = 0.7)

mtext("Spearman rho",
      side = 4,
      line = 4,
      cex = 0.8)

par(xpd = FALSE)

dev.off()

cat("Correlation heatmap saved: ",
    file.path(diagnostic_dir, "diagnostic_spearman_correlation_heatmap.pdf"),
    "\n")

# END OF RI DIAGNOSTIC

cat("END OF RI DIAGNOSTIC !")


#PLOTS FOR PAPER WORK#
# FIGURE 7 : TOP 6 : environmental maps + response curves + RI + Coefficient of variation

if (!require(png)) install.packages("png")
library(png)

# 1. Select the top 6 variables by mean relative importance

ri_means  <- colMeans(relative_importances, na.rm = TRUE)
top6_vars <- names(sort(ri_means, decreasing = TRUE))[1:6]
top6_RI   <- sort(ri_means, decreasing = TRUE)[1:6]

cat("Top 6 variables by RI:\n")
print(data.frame(variable = top6_vars, mean_RI = round(top6_RI, 2)))


# 2. Variable titles

title_lookup <- list(
  temperature_spring             = c("Air temperature",   "spring (°C)"),
  temperature_summer             = c("Air temperature",   "summer (°C)"),
  temperature_inFall             = c("Air temperature",   "autumn (°C)"),
  temperature_winter             = c("Air temperature",   "winter (°C)"),
  precipitation_spring           = c("Precipitation",     "spring (kg/m²/day)"),
  precipitation_summer           = c("Precipitation",     "summer (kg/m²/day)"),
  precipitation_inFall           = c("Precipitation",     "autumn (kg/m²/day)"),
  precipitation_winter           = c("Precipitation",     "winter (kg/m²/day)"),
  relative_humidity_spring       = c("Relative humidity", "spring (%)"),
  relative_humidity_summer       = c("Relative humidity", "summer (%)"),
  relative_humidity_inFall       = c("Relative humidity", "autumn (%)"),
  relative_humidity_winter       = c("Relative humidity", "winter (%)"),
  primary_forest_areas           = c("Forested",          "primary land"),
  primary_non.forest_areas       = c("Primary",           "non-forest land"),
  secondary_forest_areas         = c("Forested",          "secondary land"),
  secondary_non.forest_areas     = c("Secondary",         "non-forest land"),
  croplands_all_categories       = c("Croplands",         "(all categories)"),
  managed_pasture_and_rangeland  = c("Pastures",          "and rangeland"),
  human_pop_density_log10        = c("Human population",  "(log10)")
)


# 3. Color palettes

get_cols_top6 <- function(name) {
  if (grepl("temperature", name)) {
    return(colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100))
  }
  if (grepl("precipitation", name)) {
    return(colorRampPalette(brewer.pal(9, "YlGnBu"))(100))
  }
  if (grepl("humidity", name)) {
    return(colorRampPalette(brewer.pal(9, "PuBu"))(100))
  }
  if (grepl("primary_forest", name)) {
    return(colorRampPalette(brewer.pal(9, "Greens"))(100))
  }
  if (grepl("primary_non", name)) {
    return(colorRampPalette(brewer.pal(9, "YlGn"))(100))
  }
  if (grepl("secondary_forest", name)) {
    return(colorRampPalette(brewer.pal(9, "BuGn"))(100))
  }
  if (grepl("secondary_non", name)) {
    return(colorRampPalette(brewer.pal(9, "YlGn"))(100))
  }
  if (grepl("cropland", name)) {
    return(colorRampPalette(brewer.pal(9, "YlOrBr"))(100))
  }
  if (grepl("pasture", name)) {
    return(colorRampPalette(brewer.pal(9, "Oranges"))(100))
  }
  if (grepl("human", name)) {
    return(colorRampPalette(brewer.pal(9, "Purples"))(100))
  }
  return(colorRampPalette(brewer.pal(9, "Greys"))(100))
}


# 4. Helper function: response curves

compute_rc <- function(brt_models, varname, n_grid = 100) {
  
  all_train <- do.call(rbind, lapply(brt_models, function(m) {
    m$gbm.call$dataframe[, m$gbm.call$predictor.names]
  }))
  
  preds <- brt_models[[1]]$gbm.call$predictor.names
  
  v_range  <- range(all_train[, varname], na.rm = TRUE)
  grid_val <- seq(v_range[1], v_range[2], length.out = n_grid)
  
  fixed <- apply(all_train[, preds], 2, median, na.rm = TRUE)
  
  pred_df <- as.data.frame(matrix(
    rep(fixed, each = n_grid),
    nrow = n_grid,
    ncol = length(fixed)
  ))
  
  colnames(pred_df) <- preds
  pred_df[, varname] <- grid_val
  
  y_list <- lapply(brt_models, function(m) {
    n_t <- m$gbm.call$best.trees
    predict.gbm(
      m,
      pred_df,
      n.trees = n_t,
      type = "response",
      single.tree = FALSE
    )
  })
  
  list(x = grid_val, y = y_list)
}


# 5. Helper function: variable-specific but readable y-axis limits

response_ylim <- function(y_list) {
  
  y <- unlist(y_list)
  y <- y[is.finite(y)]
  
  rng <- range(y, na.rm = TRUE)
  pad <- diff(rng) * 0.12
  
  if (!is.finite(pad) || pad == 0) {
    pad <- 0.05
  }
  
  ylim <- c(rng[1] - pad, rng[2] + pad)
  ylim[1] <- max(0, ylim[1])
  ylim[2] <- min(1, ylim[2])
  
  if (diff(ylim) < 0.15) {
    mid <- mean(ylim)
    ylim <- c(max(0, mid - 0.075), min(1, mid + 0.075))
  }
  
  return(ylim)
}


# 6. Helper function: side labels in boxes

draw_side_box <- function(label, cex = 0.8) {
  par(mar = c(0, 0, 0, 0))
  plot.new()
  rect(0.16, 0.05, 0.84, 0.95, border = "gray55", lwd = 0.8)
  text(0.50, 0.50, label, srt = 90, cex = cex, col = "gray20")
}


# 7. Render each environmental map as a temporary high-resolution PNG
#    This avoids raster overflow and keeps each map clipped to its own panel.

tmp_dir <- tempfile("top6_maps_")
dir.create(tmp_dir)

map_pngs <- character(length(top6_vars))

for (k in seq_along(top6_vars)) {
  
  v <- top6_vars[k]
  idx <- which(envVariableNames == v)
  cols <- get_cols_top6(v)
  r <- envVariables[[idx]]
  
  map_pngs[k] <- file.path(tmp_dir, paste0(v, ".png"))
  
  png(map_pngs[k],
      width = 1200,
      height = 720,
      res = 180,
      bg = "white")
  
  par(mar = c(2.2, 0.2, 0.2, 0.2),
      lwd = 0.3,
      col = "gray30",
      col.axis = "gray30",
      fg = "gray30")
  
  plot(r,
       col = cols,
       ann = FALSE,
       legend = FALSE,
       axes = FALSE,
       box = FALSE,
       asp = NA)
  
  plot(coasts, add = TRUE, border = "gray50", lwd = 0.30)
  
  points(occurrences$longitude,
         occurrences$latitude,
         pch = 21,
         cex = 0.42,
         col = "black",
         bg = "red",
         lwd = 0.25)
  
  plot(r,
       col = cols,
       legend.only = TRUE,
       add = TRUE,
       horizontal = TRUE,
       legend.width = 0.45,
       legend.shrink = 0.72,
       smallplot = c(0.18, 0.82, 0.055, 0.095),
       axis.args = list(
         cex.axis = 0.58,
         lwd = 0,
         lwd.tick = 0.22,
         col.tick = "gray30",
         tck = -0.45,
         col.axis = "gray30",
         line = 0,
         mgp = c(0, 0.12, 0)
       ))
  
  dev.off()
}


# 8. Build final composite figure

if (!dir.exists("outputs")) {
  dir.create("outputs", recursive = TRUE)
}

pdf("outputs/Fig7_top6_maps_curves_RI.pdf",
    width = 17,
    height = 6.8)

layout(
  matrix(1:21, nrow = 3, byrow = TRUE),
  widths  = c(0.52, rep(1, 6)),
  heights = c(1.35, 0.90, 0.13)
)

par(oma = c(0.2, 0.2, 0.2, 0.2),
    lwd = 0.3,
    col = "gray30",
    col.axis = "gray30",
    fg = "gray30")

red_curve <- rgb(222, 67, 39, maxColorValue = 255)


# Row 1: environmental maps

draw_side_box("Environmental factors used for ENM", cex = 1.60)

for (k in seq_along(top6_vars)) {
  
  ttl <- title_lookup[[top6_vars[k]]]
  img <- png::readPNG(map_pngs[k])
  
  par(mar = c(0.1, 0.1, 3.0, 0.1))
  plot.new()
  
  rasterImage(img, 0, 0, 1, 1)
  
  mtext(ttl[1],
        side = 3,
        line = 1.20,
        cex = 1.05,
        col = "gray15")
  
  mtext(ttl[2],
        side = 3,
        line = 0.05,
        cex = 1.05,
        col = "gray15")
}

# Row 2: response curves

draw_side_box("Response curves", cex = 1.60)

for (k in seq_along(top6_vars)) {
  
  v <- top6_vars[k]
  
  rc <- compute_rc(brt_models, v, n_grid = 100)
  ylim_v <- response_ylim(rc$y)
  
  par(mar = c(3.2, 3.2, 0.6, 0.8),
      pin = c(1.55, 1.55))
  
  plot(NA,
       xlim = range(rc$x, na.rm = TRUE),
       ylim = ylim_v,
       xlab = "",
       ylab = "",
       bty = "n",
       las = 1,
       cex.axis = 0.68,
       tck = -0.03,
       mgp = c(2, 0.42, 0),
       xaxs = "i",
       yaxs = "i")
  
  box(lwd = 0.4, col = "gray55")
  
  for (j in seq_along(rc$y)) {
    lines(rc$x,
          rc$y[[j]],
          col = adjustcolor(red_curve, alpha.f = 0.55),
          lwd = 0.45)
  }
  
  title(xlab = "environmental values",
        cex.lab = 1.0,
        mgp = c(1.55, 0, 0),
        col.lab = "gray20")
  
  title(ylab = "predicted values",
        cex.lab = 1.0,
        mgp = c(2.05, 0, 0),
        col.lab = "gray20")
}

# Row 3: RI values

par(mar = c(0, 0, 0, 0))
plot.new()
text(0.50, 0.72, "(RI:)", cex = 1.70, col = "black", font = 2)

for (k in seq_along(top6_vars)) {
  
  par(mar = c(0, 0, 0, 0))
  plot.new()
  
  text(0.50,
       0.72,
       sprintf("%.1f%%", top6_RI[k]),
       col = red_curve,
       font = 2,
       cex = 1.30)
}

dev.off()

cat("\nFigure saved: outputs/Fig7_top6_maps_curves_RI.pdf\n")


# FIGURE 8) Example of presence and pseudo-absence distribution


# This figure shows the spatial distribution of presence and pseudo-absence
# points for one BRT replicate. It is useful to document the training design
# used by the ENM, since the model learns ecological suitability by contrasting
# presence locations with pseudo-absence/background locations.

# Select one replicate as an example
replicate_to_plot <- 1
d_example <- dataframes[[replicate_to_plot]]

pdf("outputs/Fig8_presence_pseudoabsence_distribution.pdf",
    width = 7, height = 7)

par(mar = c(0.5, 0.5, 2.5, 0.5),
    lwd = 0.3,
    col = "gray30",
    col.axis = "gray30",
    fg = "gray30")

plot(coasts,
     col = "gray95",
     border = "gray50",
     lwd = 0.5,
     main = "")

points(d_example$longitude[d_example$response == 0],
       d_example$latitude[d_example$response == 0],
       pch = 21,
       cex = 0.55,
       col = "gray30",
       bg = adjustcolor("#2C7FB8", alpha.f = 0.55),
       lwd = 0.25)

points(d_example$longitude[d_example$response == 1],
       d_example$latitude[d_example$response == 1],
       pch = 21,
       cex = 0.65,
       col = "gray20",
       bg = adjustcolor("#D7301F", alpha.f = 0.85),
       lwd = 0.3)

legend("bottomleft",
       legend = c("Presences", "Pseudo-absences"),
       pch = 21,
       pt.bg = c(adjustcolor("#D7301F", alpha.f = 0.85),
                 adjustcolor("#2C7FB8", alpha.f = 0.55)),
       col = c("gray20", "gray30"),
       pt.cex = c(1.1, 1.0),
       bty = "n",
       cex = 0.85)

dev.off()

cat("Presence/pseudo-absence distribution figure saved: outputs/Fig8_presence_pseudoabsence_distribution.pdf\n")

#END OF PAPER PLOTS 

cat("END OF PAPER PLOTS ! THANK YOU FOR READING !")


