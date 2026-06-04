# Dengue-risque-mapping-Europe

This repository contains the R script, input files and selected outputs used for ecological suitability and risk mapping of dengue virus circulation in Europe.

## Repository structure

- `Projet1_Risk_Mapping_Code.R`: main R script used for data loading, exploratory analyses, ecological niche modelling, model evaluation and risk map production.
- `data/occurrences/`: simulated arbovirus occurrence dataset.
- `data/rasters/`: environmental raster predictors used for ecological niche modelling.
- `data/shapefiles/`: coastline shapefile used for map visualisation.
- `outputs/`: selected model outputs, figures and exported rasters.
- `outputs/diagnostic/`: diagnostic outputs used to interpret model behaviour and predictor importance.

## Main analyses

The workflow includes:

- exploration of occurrence data and environmental predictors;
- ecological niche modelling using Boosted Regression Trees (BRT);
- pseudo-absence sampling and replicated model training;
- spatial cross-validation using spatial blocks;
- model performance assessment using AUC and SIppc;
- production of a mean ecological suitability map and uncertainty map;
- diagnostic analyses of predictor effects, relative importance and environmental correlations.

## Ecological suitability map

The final risk map summarises the mean predicted ecological suitability for local dengue virus circulation across Europe, together with spatial uncertainty across BRT replicates.

![Ecological suitability map](outputs/Fig5_cartographie_risque.png) 

The PDF version is available in:

- `outputs/Fig5_cartographie_risque.pdf`

## Key outputs

Selected outputs include:

- `outputs/Fig1_exploration_occurrences.pdf`: spatial distribution of occurrence points and coordinate distributions.
- `outputs/Fig2_19variables.pdf`: maps of the 19 environmental predictors with occurrence points.
- `outputs/Fig3_valeurs_env_aux_presences.pdf`: environmental values observed at presence locations.
- `outputs/Fig4_correlogramme_spatial.pdf`: spatial correlogram used to guide spatial cross-validation.
- `outputs/Fig5_cartographie_risque.pdf`: final ecological suitability / risk map.
- `outputs/Fig5_cartographie_risque.png`: PNG version of the final risk map displayed in this README.
- `outputs/Fig6_courbes_reponse.pdf`: response curves showing how predicted suitability varies with environmental predictors.
- `outputs/Fig7_top6_maps_curves_RI.pdf`: maps, response curves and relative importance for the main predictors.
- `outputs/Fig8_presence_pseudoabsence_distribution.pdf`: distribution of presences and pseudo-absences used for model training.
- `outputs/BRT_model_predictions.rds`: saved BRT prediction rasters used to compute the final suitability map.

Diagnostic outputs include:

- `outputs/diagnostic/diagnostic_univariate_AUCs_clean_with_SD.pdf`: univariate discriminative power of predictors.
- `outputs/diagnostic/diagnostic_paired_AUC_full_vs_reduced.pdf`: comparison between the full model and the model without winter relative humidity.
- `outputs/diagnostic/diagnostic_delta_AUC_distribution.pdf`: distribution of AUC loss after removing winter relative humidity.
- `outputs/diagnostic/diagnostic_RI_redistribution_all_variables.pdf`: redistribution of relative importance after removing winter relative humidity.
- `outputs/diagnostic/diagnostic_spearman_correlation_heatmap.pdf`: Spearman correlation structure among environmental predictors.



## How to run the R script

Before running the script, set the working directory to the repository root folder, which must contain `data/` and `outputs/`.

In R:

```r
source("Projet1_Risk_Mapping_Code.R")
