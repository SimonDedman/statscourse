# TidyModels for Species Distribution Modelling: A BRT/xgBoost Workflow ####
# Simon Dedman — Florida International University — March 04, 2026
#
# This script contains all audience-facing R code examples from the
# presentation "TidyModels for Species Distribution Modelling".
# Code is grouped by slide and presented in presentation order.
#
# Sections are standalone concept demos (SLIDE 9-37) and a complete
# runnable workflow (SLIDE 40-42). The concept demos each illustrate
# one package; they are NOT meant to be run sequentially.
#
# Data: samples.rds (2,244 survey points, training data)
#       grids.rds (378,570 prediction grid cells)
#
# Packages required:
#   tidymodels, rsample, spatialsample, recipes, parsnip, workflows,
#   dials, tune, yardstick, vip, DALEXtra, themis, tidysdm, CAST,
#   embarcadero, BART, dbarts, sf, terra, probably, car, dplyr
#
# Install (sub-packages of tidymodels/tidyverse excluded):
# install.packages(c(
#   "tidymodels",    # meta-package: rsample, recipes, parsnip, workflows,
#                    #   dials, tune, yardstick, dplyr, ggplot2, etc.
#   "tidyverse",     # meta-package: dplyr, ggplot2, tidyr, purrr, readr, etc.
#   "spatialsample", # spatial CV folds
#   "sf",            # spatial features
#   "terra",         # raster data
#   "vip",           # variable importance plots
#   "DALEXtra",      # partial dependence (installs DALEX)
#   "themis",        # SMOTE / class imbalance
#   "CAST",          # area of applicability
#   "car",           # VIF calculation
#   "probably",      # conformal prediction intervals
#   "tidysdm",       # SDM-specific tools
#   "embarcadero",   # BART wrapper (installs dbarts)
#   "BART"           # Bayesian Additive Regression Trees
# ))

# Data Loading and Column Renaming ####
# Rename columns to match code examples throughout.
# LONGITUDE → lon, LATITUDE → lat, Depth → depth, Temperature → sst
# No 'chlorophyll' analog — SLIDE 24 PDP uses Salinity instead.
# Response: presence/absence from Cuckoo ray CPUE (12.3% prevalence).

setwd(
  "/home/simon/Documents/Si Work/PostDoc Work/FIU/2026-01 Al Harborne Stats Class/2026-03-04 BRTs TidyModel Workflow/Code/"
)

samples <- readRDS("samples.rds") |>
  dplyr::rename(
    lon = LONGITUDE,
    lat = LATITUDE,
    depth = Depth,
    sst = Temperature
  ) |>
  dplyr::mutate(
    presence = factor(
      ifelse(Cuckoo > 0, "present", "absent"),
      levels = c("absent", "present")
    )
  ) |>
  # Keep only coordinates, predictors, and response.
  # Drops: Survey_StNo_HaulNo_Year (ID), Cuckoo (source of presence),
  #        Thornback, Blonde, Spotted (other species responses).
  dplyr::select(lon, lat, depth, sst, Salinity, Current_Speed,
                Grain_Size, Distance_to_Shore, F_LPUE, Effort, presence)

grids <- readRDS("grids.rds") |>
  dplyr::rename(
    lon = LONGITUDE,
    lat = LATITUDE,
    depth = Depth,
    sst = Temperature
  )


# SEDI (Symmetric Extremal Dependence Index) — Custom yardstick metric ####
# Not in yardstick; preferred over MCC/TSS when prevalence < 2.5%
# (Wunderlich et al. 2019; Ferro & Stephenson 2011).
# Defined here so it can be used in metric_set(), tune_grid(), select_best().

sedi_vec <- function(truth, estimate, na_rm = TRUE, case_weights = NULL,
                     event_level = "first", ...) {
  yardstick:::abort_if_class_pred(truth)
  estimate <- yardstick:::as_factor_from_class_pred(estimate)
  estimator <- yardstick:::finalize_estimator(truth, estimator = "binary")
  yardstick:::check_class_metric(truth, estimate, case_weights, estimator)
  if (na_rm) {
    result <- yardstick:::yardstick_remove_missing(truth, estimate, case_weights)
    truth <- result$truth
    estimate <- result$estimate
    case_weights <- result$case_weights
  } else if (yardstick:::yardstick_any_missing(truth, estimate, case_weights)) {
    return(NA_real_)
  }
  data <- yardstick:::yardstick_table(truth, estimate, case_weights = case_weights)
  sens <- yardstick:::sens_binary(data, event_level)
  spec <- yardstick:::spec_binary(data, event_level)
  small <- 1e-9
  H <- max(min(sens, 1 - small), small)
  Fa <- max(min(1 - spec, 1 - small), small)
  (log(Fa) - log(H) - log(1 - Fa) + log(1 - H)) /
    (log(Fa) + log(H) + log(1 - Fa) + log(1 - H))
}

sedi <- function(data, ...) { UseMethod("sedi") }
sedi.data.frame <- function(data, truth, estimate, na_rm = TRUE,
                            case_weights = NULL, event_level = "first", ...) {
  yardstick:::class_metric_summarizer(
    name = "sedi", fn = sedi_vec,
    data = data, truth = !!rlang::enquo(truth), estimate = !!rlang::enquo(estimate),
    na_rm = na_rm, case_weights = !!rlang::enquo(case_weights),
    event_level = event_level
  )
}
sedi <- yardstick::new_class_metric(sedi, direction = "maximize")


# SLIDE 9 — rsample: Standard Data Splitting ####
# Hold out data the model hasn't seen, so we can test whether it generalises.

library(rsample)
library(dplyr)

# 80/20 split, stratified by outcome
split <- rsample::initial_split(samples, prop = 0.8, strata = presence)
train_data <- rsample::training(split)
test_data <- rsample::testing(split)

# Auto-select evaluation metric based on class count and prevalence
n_classes <- nlevels(train_data$presence)
if (n_classes > 2) {
  class_prevs <- table(train_data$presence) / nrow(train_data)
  min_prev <- min(class_prevs)
  message(sprintf("Multiclass classification (%d classes); min prevalence: %.1f%%",
                  n_classes, min_prev * 100))
} else {
  min_prev <- mean(train_data$presence == "present")
  message(sprintf("Binary classification; prevalence: %.1f%%", min_prev * 100))
}
if (min_prev < 0.025) {
  selection_metric <- "sedi"
  message("Prevalence < 2.5%: using SEDI for model selection (Wunderlich et al. 2019)")
} else {
  selection_metric <- "mcc"
  message("Using MCC for model selection (Chicco & Jurman 2020)")
}

# 10-fold cross-validation
folds <- rsample::vfold_cv(train_data, v = 10, strata = presence)




# SLIDE 10 — spatialsample: Critical for SDMs! ####
# Standard CV ignores spatial autocorrelation — nearby points are often more
# similar (Tobler's first law). Without spatial CV, models appear better than
# they truly are.

library(spatialsample)
library(sf)

# remove = FALSE keeps lon/lat as columns (recipe needs them for update_role).
# Project to UTM 29N: cellsize is in metres, not degrees (EPSG:4326 would fail).
train_sf <- sf::st_as_sf(train_data, coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
  sf::st_transform(crs = 32629)
spatial_folds <- spatialsample::spatial_block_cv(
  train_sf,
  v = 10,
  cellsize = 50000
)
ggplot2::autoplot(spatial_folds)




# SLIDE 12 — recipes: The Preprocessing Engine ####
# Define a reusable pipeline of variable selection and transformation steps.
# prep() and bake() are called automatically by workflows — you rarely need
# them directly.

library(recipes)

# Non-predictor columns (ID, other species) already removed in preamble select().
rec <- recipes::recipe(presence ~ ., data = train_data) |>
  recipes::update_role(lon, lat, new_role = "coordinates") |> # exclude coords from model
  recipes::step_impute_median(all_numeric_predictors()) |> # fill NAs with medians
  recipes::step_nzv(all_predictors()) |> # drop near-zero variance
  recipes::step_corr(all_numeric_predictors(), threshold = 0.9) |> # drop pairwise-correlated
  recipes::step_normalize(all_numeric_predictors()) |> # centre & scale
  recipes::step_dummy(all_nominal_predictors()) # one-hot encode factors




# SLIDE 14 — parsnip: Unified Model Interface ####
# A single interface to many model types — swap engines without rewriting code.
# Key concepts: Type (model category), Mode (task type), Engine (underlying
# package). tune() marks hyperparameters to be optimised later.

library(parsnip)

brt_spec <- parsnip::boost_tree(
  trees = 1000,
  tree_depth = tune::tune(),
  min_n = tune::tune(),
  learn_rate = tune::tune(),
  mtry = tune::tune(),
  sample_size = tune::tune(),
  loss_reduction = tune::tune()
) |>
  parsnip::set_engine("xgboost") |>
  parsnip::set_mode("classification")
# tune(): placeholder for tune_grid() [SLIDE 18]




# SLIDE 17 — workflows: Bundle Everything ####
# Bundle recipe + model into one object so fit() handles everything.
# Advantages: single fit() handles prep + train; consistent preprocessing at
# prediction; simplified tuning integration.
# NOTE: fit() below would fail because brt_spec has tune() placeholders.
# This is a concept demo; see Complete Workflow (SLIDE 40) for the full sequence.

library(workflows)

brt_wf <- workflows::workflow() |>
  workflows::add_recipe(rec) |>
  workflows::add_model(brt_spec)

# brt_fit <- fit(brt_wf, data = train_data)  # requires finalized hyperparameters
# preds <- predict(brt_fit, new_data = test_data)




# SLIDE 18 — dials & tune: Hyperparameter Tuning ####
# Create a hyperparameter search grid and tune across spatial CV folds.
# Recommended strategy for SDMs:
#   1. Start: Space-filling grid (30-50 points)
#   2. If needed: Racing methods for larger grids
#   3. Refine: Bayesian optimisation around best regions
#   4. Select: select_by_one_std_err() — picks the simplest model whose
#      performance is within one standard error of the best (parsimony over
#      raw score).

library(dials)
library(tune)

brt_params <- brt_wf |>
  workflows::extract_parameter_set_dials() |>
  dials::finalize(train_data)

brt_grid <- dials::grid_space_filling(brt_params, size = 30)

tune_results <- tune::tune_grid(
  brt_wf,
  resamples = spatial_folds,
  grid = brt_grid,
  metrics = yardstick::metric_set(yardstick::roc_auc, yardstick::accuracy),
  control = tune::control_grid(save_pred = TRUE)
)

# Pick simplest model within 1 SE of best
best <- tune::select_by_one_std_err(
  tune_results,
  metric = "roc_auc",
  dplyr::desc(learn_rate)
)
## Finalize & Predict ####
final_wf <- tune::finalize_workflow(brt_wf, best)
brt_fit <- fit(final_wf, data = train_data)
predict(brt_fit, new_data = test_data)




# SLIDE 20 — yardstick: Performance Metrics (Create Metric Set) ####
# Quantify model skill using metrics suited to presence/absence data.

library(yardstick)

sdm_metrics <- yardstick::metric_set(
  # Probability-based
  yardstick::roc_auc,             # discrimination (threshold-free)
  yardstick::pr_auc,              # precision-recall AUC (better for rare species)
  yardstick::mn_log_loss,         # calibration (lower is better)
  yardstick::brier_class,         # calibration (lower is better)
  yardstick::gain_capture,        # cumulative gains
  # Hard-class
  yardstick::accuracy,            # overall correct rate
  yardstick::j_index,             # TSS: sens + spec - 1 (Allouche et al. 2006)
                                   # bal_accuracy omitted: linear equivalent ((TSS+1)/2)
  yardstick::mcc,                 # Matthews correlation coefficient
  yardstick::kap,                 # Cohen's Kappa
  yardstick::f_meas,              # F1 score
  yardstick::sens,                # sensitivity / recall / TPR
  yardstick::spec,                # specificity / TNR
  yardstick::ppv,                 # positive predictive value / precision
  yardstick::npv,                 # negative predictive value
  yardstick::detection_prevalence,# rate of positive predictions
  sedi                            # SEDI (Wunderlich et al. 2019; defined above)
)




# SLIDE 20 — yardstick: Finalize and Evaluate ####
# Select best parameters, finalize workflow, and evaluate on held-out test set.
# NOTE: references brt_wf and split from earlier concept-demo sections.

# Select best parameters
best_params <- tune::select_best(tune_results,
                                 metric = selection_metric)

# Finalize workflow
final_workflow <- tune::finalize_workflow(brt_wf, best_params)

# Evaluate on test set
final_fit <- tune::last_fit(
  final_workflow,
  split = split,
  metrics = sdm_metrics
)

tune::collect_metrics(final_fit)




# SLIDE 23 — Variable Importance with vip ####
# Identify which predictors drive the model's predictions.
# Native importance is based on split frequency/gain; permutation importance is
# model-agnostic and measures true impact — more reliable for ecological
# interpretation. Always report which method was used.

library(vip)

fitted_model <- workflows::extract_fit_parsnip(final_fit)

# Native XGBoost importance
vip::vip(fitted_model, num_features = 15)

# Permutation importance — needs pred_wrapper (returns numeric vector)
# and a metric function with signature metric(truth, estimate).
pfun <- function(object, newdata) {
  predict(object, new_data = newdata, type = "prob")$.pred_present
}
auc_metric <- function(truth, estimate) {
  yardstick::roc_auc_vec(truth, estimate, event_level = "second")
}
vi_perm <- vip::vi(
  fitted_model,
  method = "permute",
  train = train_data,
  target = "presence",
  metric = auc_metric,
  pred_wrapper = pfun,
  smaller_is_better = FALSE,
  nsim = 30
)

vip::vip(vi_perm, num_features = 15)




# SLIDE 24 — Partial Dependence with DALEX ####
# Visualise how each predictor affects the predicted outcome.

library(DALEXtra)

# Create explainer — pass the fitted workflow, not the last_fit result.
explainer <- DALEXtra::explain_tidymodels(
  final_fit$.workflow[[1]],
  data = train_data |> dplyr::select(-presence),
  y = as.numeric(train_data$presence == "present"),
  label = "BRT SDM",
  verbose = FALSE
)

# Partial dependence plots
pdp_depth <- DALEX::model_profile(explainer, variables = "depth", N = 500)
plot(pdp_depth)

# Multiple variables (no chlorophyll analog in this dataset; using Salinity)
pdp_multi <- DALEX::model_profile(
  explainer,
  variables = c("depth", "sst", "Salinity"),
  N = 500
)
plot(pdp_multi)




# SLIDE 26 — Handling Class Imbalance in SDMs (themis: Recipe with SMOTE) ####
# SDM data often have far more absences than presences.
# SMOTE interpolates minority neighbours to create synthetic points.
# Note: SMOTE requires dummy encoding before application.

library(themis)

rec_balanced <- recipes::recipe(presence ~ ., data = train_data) |>
  recipes::step_impute_median(all_numeric_predictors()) |>
  themis::step_smote(presence, over_ratio = 0.5) |>
  recipes::step_normalize(all_numeric_predictors())




# SLIDE 26 — Handling Class Imbalance in SDMs (Alternative: Class Weights) ####
# For rare species with <1% prevalence, class weights in the engine are an
# alternative to resampling methods.

# Calculate ratio
neg_pos_ratio <-
  sum(train_data$presence == "absent") /
  sum(train_data$presence == "present")

# Add to engine
brt_spec_weighted <- parsnip::boost_tree(trees = 1000) |>
  parsnip::set_engine("xgboost", scale_pos_weight = neg_pos_ratio) |>
  parsnip::set_mode("classification")




# SLIDE 28 — Spatial Data Packages ####
# Build a raster environmental stack from grids.rds, convert occurrences
# to sf, and extract environmental values at occurrence points.

library(terra)
library(sf)

# -- Build env_stack from grids.rds --
# grids has irregular spacing, so we rasterize onto a regular 0.01° template.
env_cols <- c("depth", "sst", "Salinity", "Current_Speed",
              "Grain_Size", "Distance_to_Shore", "F_LPUE", "Effort")
grids_vect <- terra::vect(grids, geom = c("lon", "lat"), crs = "EPSG:4326")
# Regular template covering the grids extent
template <- terra::rast(terra::ext(grids_vect), resolution = 0.01,
                        crs = "EPSG:4326")
env_stack <- terra::rasterize(grids_vect, template, field = env_cols,
                              fun = mean)

# -- Build occurrences_sf (presence-only points) --
occurrences_sf <- samples |>
  dplyr::filter(presence == "present") |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# -- Extract environmental values at occurrence points --
env_values <- terra::extract(env_stack, terra::vect(occurrences_sf),
                             ID = FALSE)




# SLIDE 29 — SDM-Specific Packages: tidysdm Example ####
# Occurrence thinning, pseudo-absence generation, and threshold calibration.
# Random pseudoabsence generation within the activity space performs better
# than random walk (Hazen et al. 2021).
# Uses occurrences_sf and env_stack built in SLIDE 28.

library(tidysdm)

# Thin occurrences to one per raster cell
thinned <- tidysdm::thin_by_cell(occurrences_sf, raster = env_stack,
                                 coords = c("lon", "lat"))

# Generate pseudo-absences (~4:1 ratio for demo)
# thin_by_cell renames coords to X/Y, which sample_pseudoabs auto-detects.
pseudoabs <- tidysdm::sample_pseudoabs(
  thinned,
  n = 1000,
  raster = env_stack,
  method = "random"
)

# Threshold calibration (requires a fitted tidysdm ensemble — not built here)
# calibrated <- tidysdm::calib_class_thresh(sdm_ensemble, class_thresh = "tss")




# SLIDE 31 — Area of Applicability (CAST) ####
# Compare model predictions to training data coverage — flag extrapolation.
# Identifies regions where your model is extrapolating beyond training
# conditions. gbm.auto::gbm.rsb() provides a related approach.

library(CAST)

# Calculate area of applicability
aoa_result <- CAST::aoa(
  newdata = env_df,
  model = final_fit,
  trainDI = CAST::trainDI(model = final_fit, train = train_data)
)

# Visualize where predictions are reliable
plot(aoa_result)

# Mask predictions outside AOA
pred_masked <- pred_rast
pred_masked[aoa_result$AOA == 0] <- NA




# SLIDE 33 — Variable Selection: VIF ####
# Use VIF < 10 threshold for collinearity (Yavas et al. 2024).
# car::vif() is more thorough than step_corr() — it accounts for
# multicollinearity from all predictors simultaneously, not just pairwise.
# Run BEFORE building the recipe to decide which variables to keep.

library(car)

# Subset to environmental predictors only (exclude coordinates, ID, species)
vif_data <- train_data |>
  dplyr::select(
    presence,
    depth,
    sst,
    Salinity,
    Current_Speed,
    Grain_Size,
    Distance_to_Shore
  )
vif_model <- glm(presence ~ ., data = vif_data, family = binomial)
car::vif(vif_model)
# Drop variables with VIF > 10, then build recipe with remaining predictors




# SLIDE 37 — BART: Bayesian Alternative to BRT ####
# BART provides full posterior uncertainty quantification, minimal tuning, and
# natural causal inference. Slots directly into the tidymodels workflow via
# parsnip::bart() with the dbarts engine.

library(embarcadero) # Carlson 2020
# Or use the base BART package:
library(BART)
bart_fit <- BART::pbart(x.train, y.train)

# BART in TidyModels via parsnip
bart_spec <- parsnip::bart(trees = 200) |>
  parsnip::set_engine("dbarts") |>
  parsnip::set_mode("classification")




# SLIDE 37 — Uncertainty in xgboost/BRTs ####
# Bootstrapping, quantile regression, and conformal prediction.

# Bootstrapping: fit multiple models on resampled data, get prediction intervals.
# For PDPs, VIP barplots, spatial predictions.
# See gbm.auto::gbm.loop (Thesis S5.2.4.4 & 6.6.4)

# Quantile regression: fit xgboost models on 5th/50th/95th percentiles

# Conformal prediction: post-hoc prediction intervals via probably
# library(probably)
# See https://probably.tidymodels.org/




# SLIDE 40 — Complete Workflow (Part 1) ####
# Full reproducible SDM workflow from data preparation through spatial CV folds.

library(tidymodels)
library(spatialsample)
library(sf)
library(terra)
library(vip)
library(themis)
tidymodels::tidymodels_prefer()
# Resolves function name conflicts (e.g. dplyr::filter vs stats::filter) in
# favour of tidymodels packages. Redundant when using package:: prefixes but
# good practice for interactive use.

## 1. Prepare data ####
# samples already cleaned in preamble: only coords, predictors, and presence remain.

## 2. Split data ####
set.seed(42)
data_split <- rsample::initial_split(samples, prop = 0.8, strata = presence)
train_data <- rsample::training(data_split)
test_data <- rsample::testing(data_split)

# Auto-select evaluation metric based on class count and prevalence
n_classes <- nlevels(train_data$presence)
if (n_classes > 2) {
  class_prevs <- table(train_data$presence) / nrow(train_data)
  min_prev <- min(class_prevs)
  message(sprintf("Multiclass classification (%d classes); min prevalence: %.1f%%",
                  n_classes, min_prev * 100))
} else {
  min_prev <- mean(train_data$presence == "present")
  message(sprintf("Binary classification; prevalence: %.1f%%", min_prev * 100))
}
if (min_prev < 0.025) {
  selection_metric <- "sedi"
  message("Prevalence < 2.5%: using SEDI for model selection (Wunderlich et al. 2019)")
} else {
  selection_metric <- "mcc"
  message("Using MCC for model selection (Chicco & Jurman 2020)")
}

## 3. Spatial CV folds ####
# remove = FALSE keeps lon/lat as columns (recipe needs them for update_role).
# Project to UTM 29N: cellsize is in metres, not degrees (EPSG:4326 would fail).
train_sf <- sf::st_as_sf(train_data, coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
  sf::st_transform(crs = 32629)
spatial_folds <- spatialsample::spatial_block_cv(
  train_sf,
  v = 10,
  cellsize = 50000
)




# SLIDE 41 — Complete Workflow (Part 2) ####
# Recipe, model specification, workflow bundling, and tuning grid.

## 4. Recipe ####
sdm_recipe <- recipes::recipe(presence ~ ., data = train_data) |>
  recipes::update_role(lon, lat, new_role = "coordinates") |>
  recipes::step_impute_median(all_numeric_predictors()) |>
  recipes::step_nzv(all_predictors()) |>
  recipes::step_corr(all_numeric_predictors(), threshold = 0.9) |>
  recipes::step_normalize(all_numeric_predictors())

## 5. Model specification ####
brt_spec <- parsnip::boost_tree(
  trees = 1000,
  tree_depth = tune::tune(),
  min_n = tune::tune(),
  learn_rate = tune::tune(),
  mtry = tune::tune()
) |>
  parsnip::set_engine("xgboost") |>
  parsnip::set_mode("classification")

## 6. Workflow ####
brt_workflow <- workflows::workflow() |>
  workflows::add_recipe(sdm_recipe) |>
  workflows::add_model(brt_spec)

## 7. Tuning grid ####
brt_params <- brt_workflow |>
  workflows::extract_parameter_set_dials() |>
  dials::finalize(train_data)
brt_grid <- dials::grid_space_filling(brt_params, size = 30)




# SLIDE 42 — Complete Workflow (Part 3) ####
# Tuning, best-parameter selection, final fit, evaluation, and raster
# prediction — the final steps to a deployable SDM.

## 8. Tune with spatial CV ####
sdm_metrics <- yardstick::metric_set(
  # Probability-based
  yardstick::roc_auc,             # discrimination (threshold-free)
  yardstick::pr_auc,              # precision-recall AUC (better for rare species)
  yardstick::mn_log_loss,         # calibration (lower is better)
  yardstick::brier_class,         # calibration (lower is better)
  yardstick::gain_capture,        # cumulative gains
  # Hard-class
  yardstick::accuracy,            # overall correct rate
  yardstick::j_index,             # TSS: sens + spec - 1 (Allouche et al. 2006)
                                   # bal_accuracy omitted: linear equivalent ((TSS+1)/2)
  yardstick::mcc,                 # Matthews correlation coefficient
  yardstick::kap,                 # Cohen's Kappa
  yardstick::f_meas,              # F1 score
  yardstick::sens,                # sensitivity / recall / TPR
  yardstick::spec,                # specificity / TNR
  yardstick::ppv,                 # positive predictive value / precision
  yardstick::npv,                 # negative predictive value
  yardstick::detection_prevalence,# rate of positive predictions
  sedi                            # SEDI (Wunderlich et al. 2019; defined above)
)

tune_results <- tune::tune_grid(
  brt_workflow,
  resamples = spatial_folds,
  grid = brt_grid,
  metrics = sdm_metrics
)

## 9. Select best and finalize ####
best_params <- tune::select_best(tune_results,
                                 metric = selection_metric)
final_workflow <- tune::finalize_workflow(brt_workflow, best_params)

## 10. Final fit and evaluation ####
final_fit <- tune::last_fit(
  final_workflow,
  split = data_split,
  metrics = sdm_metrics
)
final_metrics <- tune::collect_metrics(final_fit)
final_metrics

## 11. Predict to grids ####
# grids.rds serves as the prediction surface (environmental values at each cell).
# Select only the columns the model was trained on.
env_df <- grids |>
  dplyr::select(
    lon,
    lat,
    depth,
    sst,
    Salinity,
    Current_Speed,
    Grain_Size,
    Distance_to_Shore,
    F_LPUE,
    Effort
  )
predictions <- predict(
  final_fit$.workflow[[1]],
  new_data = env_df,
  type = "prob"
)

# Convert to sf points then rasterize to a regular grid.
# terra::rast(type = "xyz") requires perfectly regular spacing;
# grids.rds coordinates are irregular, so we rasterize instead.
pred_sf <- sf::st_as_sf(
  dplyr::bind_cols(env_df[, c("lon", "lat")], predictions),
  coords = c("lon", "lat"),
  crs = 4326
)
# Create a template raster at ~0.01° resolution (~1 km)
template <- terra::rast(terra::ext(terra::vect(pred_sf)),
                        resolution = 0.01, crs = "EPSG:4326")
pred_raster <- terra::rasterize(terra::vect(pred_sf), template,
                                field = ".pred_present", fun = mean)
terra::writeRaster(pred_raster, "species_prediction.tif", overwrite = TRUE)




# SLIDE 43 — Parameter Mapping: dismo/gbm to TidyModels ####
# Reference for translating classic BRT parameters to parsnip equivalents.
#
# Classic dismo:
#   gbm.step(data, gbm.x, gbm.y, family = "bernoulli",
#            tree.complexity = 5, learning.rate = 0.01, bag.fraction = 0.75)
#
# gbm.auto:
#   gbm.auto(samples = data, expvar = 2:10, resvar = 1,
#            tc = 5, lr = 0.01, bf = 0.75)
#
# TidyModels equivalent:

parsnip::boost_tree(
  trees = 1000, # n.trees / n.trees
  tree_depth = 5, # interaction.depth
  min_n = 10, # n.minobsinnode
  learn_rate = 0.01, # shrinkage / learning.rate
  sample_size = 0.75 # bag.fraction
) |>
  parsnip::set_engine("xgboost") |>
  parsnip::set_mode("classification")

# Parallelisation (equivalent to n.cores in gbm.auto):
# library(doParallel)
# cl <- parallel::makeCluster(parallel::detectCores() - 1)
# doParallel::registerDoParallel(cl)
# ... run tune::tune_grid() ...
# parallel::stopCluster(cl)
