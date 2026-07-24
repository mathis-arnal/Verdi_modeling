# Notes — data_exploration.qmd

## Context

Working with a household COVID transmission dataset (`dataset_freezed2Nov25.csv`). The core goal is to estimate the **Secondary Attack Rate (SAR)** — i.e., the probability that a household contact gets infected from an index case — as accurately as possible.

**Important study condition:** Households were enrolled during the Omicron wave when Thailand had **no mandatory quarantine or isolation**. Household members continued normal daily activities (work, school, etc.) throughout follow-up. This has two direct consequences for modelling:

1. **Contact patterns are normal**, not inflated by home confinement — the Thai contact matrix (`TH_contact_data.txt`) is applicable as a prior on age-stratified within-household transmission probability.
2. **Community acquisition is a competing risk** — a contact who tests positive during follow-up may have been infected outside the household, not from the index case. The published phylogenetic analysis confirmed 12 community transmissions in 10 households. Any transmission model must include a background community hazard, and any SAR estimate must be interpreted with this in mind.

---

## What the script does

### 1. Data preparation

- Loads the dataset and creates derived variables:
  - `age_cat`: Child (≤18), Adult (18–59), Elderly (≥60)
  - `index_child`: label for whether the index case is a child or adult
- Propagates the household-level `exclude_mem` variable to all members of the same household (using `group_by(hhid)`).
- Defines variable lists for household-level and contact-level variables.

### 2. Missing data check on index cases

- Subsets to index cases only and inspects missingness with `gg_miss_var`.
- Finding: very little missing data among index cases.

---

## Main challenge: excluded household members

Some household members were not enrolled/observed in the study (`exclude_mem > 0`). Ignoring them biases the SAR because the denominator (number of exposed contacts) is underestimated.

**Two things need to be imputed for each missing member:**
1. Their **age category** (Child / Adult / Elderly)
2. Their **infection status** (Detected / Not detected)

---

## Imputation approach

### Building the imputation dataframe

A household-level summary (`imputation_df`) is constructed with:
- Observed infection counts and proportions
- Age composition of the household (counts of children, adults, elderly)
- Number of missing members

### Critical limitation: training and target populations differ

Complete households (the training set) are mostly **size 3**, while incomplete households (where imputation is needed) tend to be **larger**. This is a distribution shift between training and deployment — the OOB error rates below are measured on the training population and will likely underestimate true error on the incomplete households. This affects both models below and is the same problem flagged in `mics_data_notes.md` regarding validation.

### Imputation model 1 — Age category (Random Forest)

- Trained on **complete households only** (`include_all_members == "All members included"`).
- Uses leave-one-out (LOO) counts so the target member is excluded from its own predictors.
- Predictors: `child_nb_loo`, `adult_nb_loo`, `elderly_nb_loo`, `household_size`.
- With class weights: ~18% OOB error rate. Without: ~11%.
- **Issue noted**: the model over-predicts "Adult", even with class weights. A household structure prior could help.
- **Caveat**: OOB error measured on small complete HH — likely optimistic for larger incomplete HH.

### Imputation model 2 — Infection status (Random Forest)

- Same complete-household training set, contacts only (index excluded).
- LOO adjustments applied to infection counts and proportions.
- Full model predictors: `age_cat`, `age_index`, `n_infected_loo`, `n_observed_loo`, `prop_infected_loo` → ~10% OOB error.
- Parsimonious model (2 vars): `age_cat + prop_infected_loo` → ~9% OOB error.
- **Caveat**: same distribution shift applies — incomplete households are larger, so household exposure structure differs from training data.

### Sequential imputation function

`impute_age_sequential()` imputes missing members **one by one**, updating the household composition counts after each imputation so later members reflect earlier imputations. Applied to all households with missing members.

---

## Alternative approaches tested

| Approach | Notes |
|---|---|
| Poisson GLM (log link) | Mirrors Patumrat paper; tested for infection prediction |
| Logistic regression via `caret` | Cross-validated ROC, with `age_cat + prop_infected_loo` |
| **MICE** (`mice` package) | Predictive Mean Matching (`pmm`), 20 imputations, custom predictor matrix aligning with the RF models |

MICE is being explored as a more principled alternative to the custom sequential RF imputation.

---

## Current status / open questions

- The age imputation model tends to over-assign "Adult" — adding a prior on household composition (e.g., based on national age distributions) could improve calibration.
- The sequential imputation function is built but the full application to all households needs to be verified.
- MICE setup is in place but results haven't been evaluated yet.
- Next step: combine imputed members with observed members into `full_data` and re-estimate the SAR with the corrected denominators.
