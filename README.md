
# statscourse

Teaching materials for the FIU graduate statistics course (Spring 2026), covering data wrangling, exploration, and species distribution modelling with R.

## Lectures

| # | Topic | Slides | Script |
|---|-------|--------|--------|
| 04 | Tidy Data | [PDF slides](https://simondedman.github.io/statscourse/lectures/lecture04_data_wrangling.pdf) | [R/01_tidy-data.R](R/01_tidy-data.R) |
| 04 | Data Transformation | | [R/02_transform.R](R/02_transform.R) |
| 08 | TidyModels, BRTs & SDMs | [HTML slides](https://simondedman.github.io/statscourse/lectures/tidymodels_sdm_workflow.html) | [R/tidymodels_sdm_workflow.R](R/tidymodels_sdm_workflow.R) |
| 11 | Causal Modelling, DAGs & SEMs | [PDF slides](https://simondedman.github.io/statscourse/lectures/lecture11_causal_dags_sems.pdf) | |

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

## Lecture 11: Causal Modelling, DAGs & SEMs for Marine Scientists

Causal inference from observational data using directed acyclic graphs (DAGs) and structural equation models (SEMs). Covers the Schoolmaster et al. (2022) three-step framework and applies it to reef food-web ecology. Topics:

- **Why causal inference**: Pearl's Ladder of Causation (association, intervention, counterfactuals); correlation vs causation in marine ecology
- **DAGs**: encoding causal assumptions from domain knowledge; forks (confounders), pipes (mediators), colliders (selection bias)
- **The backdoor criterion**: identifying the minimal sufficient adjustment set with **dagitty** (`adjustmentSets()`, `impliedConditionalIndependencies()`, `localTests()`)
- **SEM frameworks**: classical covariance-based (lavaan), piecewise (piecewiseSEM), Bayesian (brms); DAG-informed regression with adjustment sets
- **Bayesian SEMs with brms**: informative priors, Student-t families, splines for non-linear ecology, `loo_R2()` for model comparison
- **Opposing DAGs as null models**: novel approach to testing top-down vs bottom-up control in food webs within the acyclic DAG framework
- **Case study**: French Polynesian coral reefs (24 reefs, 12 atolls/islands), testing the Exploitation Ecosystems Hypothesis with 11 trophic relationships fitted in each direction
- **Phylogenetic path analysis**: `phylopath` for accounting for shared evolutionary history (Aitchison et al. shark CFAR example)
- **Common pitfalls**: conditioning on colliders, controlling for mediators, treating DAGs as data-derived

### Key take-home messages

- **Draw a DAG before you run a regression** — causal assumptions should be explicit and testable
- **Use the backdoor criterion** to decide what to control for: not everything belongs in the model
- **Test your DAG** against data via implied conditional independencies
- **Bayesian SEMs** handle small marine datasets well and quantify uncertainty properly
- **Opposing DAGs** let you test directional hypotheses (top-down vs bottom-up) in food webs

## Data

Example datasets use Irish Sea survey trawl data:

- `samples.rds` (2,244 records, training) and `grids.rds` (378,570 cells, prediction surface) are required for Lecture 08 but not included in the repo due to size. Available from the course instructor.
- `sharkdata.rda` and associated files are used in Lecture 04.

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
