
# statscourse

Teaching materials for the FIU graduate statistics course (Spring 2026), covering data wrangling, exploration, and species distribution modelling with R.

## Lectures

| # | Topic | Slides | Script |
|---|-------|--------|--------|
| 01 | Tidy Data | | [R/01_tidy-data.R](R/01_tidy-data.R) |
| 02 | Data Transformation | | [R/02_transform.R](R/02_transform.R) |
| 08 | TidyModels, BRTs & SDMs | [HTML slides](https://simondedman.github.io/statscourse/lectures/tidymodels_sdm_workflow.html) | [R/tidymodels_sdm_workflow.R](R/tidymodels_sdm_workflow.R) |

## Lecture 08: TidyModels for Species Distribution Modelling

A complete workflow for building Boosted Regression Tree (BRT/xgboost) species distribution models using the tidymodels framework. Covers:

- **Data splitting** with rsample and spatialsample (spatial block CV)
- **Preprocessing** with recipes (imputation, normalisation, VIF)
- **Model specification** with parsnip (boost_tree/xgboost)
- **Hyperparameter tuning** with dials and tune
- **Evaluation** with yardstick (MCC, TSS/j_index, AUC, SEDI, and 12 other metrics)
- **Variable importance** with vip and partial dependence with DALEX
- **Class imbalance** handling with themis (SMOTE, class weights)
- **Spatial packages**: terra, sf, tidyterra, tidysdm
- **Prediction** to raster grids
- **SEDI metric**: custom yardstick implementation for low-prevalence species (< 2.5%)

### Key metric choices

- **Model selection**: MCC (Matthews correlation coefficient) — uses all four confusion matrix quadrants
- **Low prevalence (< 2.5%)**: switch to SEDI (Wunderlich et al. 2019) — prevalence-independent via log transform
- **Reporting**: AUC + TSS + MCC (standard); add SEDI for rare species

## Data

Example datasets use Irish Sea survey trawl data:

- `samples.rds` (2,244 records, training) and `grids.rds` (378,570 cells, prediction surface) are required for Lecture 08 but not included in the repo due to size. Available from the course instructor.
- `sharkdata.rda` and associated files are used in Lectures 01-02.

## Installation

```r
# install.packages("pak")
pak::pak("SimonDedman/statscourse")
```

## References

- Elith et al. (2008). [A working guide to boosted regression trees](https://doi.org/10.1111/j.1365-2656.2008.01390.x). *Journal of Animal Ecology*.
- Dedman et al. (2017). [gbm.auto: A software tool for simplifying spatial modelling and MPA planning](https://doi.org/10.1371/journal.pone.0188955). *PLOS ONE*.
- Allouche et al. (2006). [Assessing the accuracy of SDMs: TSS](https://doi.org/10.1111/j.1365-2664.2006.01214.x). *Journal of Applied Ecology*.
- Wunderlich et al. (2019). [Two alternative evaluation metrics to replace TSS](https://doi.org/10.3897/natureconservation.35.33918). *Nature Conservation*.
- Chicco & Jurman (2020). [MCC more reliable than balanced accuracy and F1](https://doi.org/10.1186/s13040-021-00244-z). *BioData Mining*.
