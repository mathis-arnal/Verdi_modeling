# Notes — mics_data.qmd

## Context

This script uses **Thailand MICS6 2022** (UNICEF survey) as an external reference dataset to perform **hot-deck imputation** of missing household members in the main study (`dataset_freezed2Nov25.csv`). The idea: since we don't know the age/sex of excluded members, we can borrow plausible household compositions from a nationally representative survey.

---

## What the script does

### 1. Build the donor pool from MICS data

- Loads the MICS household-member file (`hl.sav`) and household file (`hh.sav`).
- Filters to **Chiang Mai only** (`HH7A == 50`).
- Computes age-group counts per household: `n_child` (0–17), `n_adult` (18–59), `n_elderly` (60+).
- Cross-checks household size against the official `HH48` variable — values match, confirming consistency.
- Filters to households with **at least 2 members and at least 1 child** (to match the study population).
- Builds `comp_table`: frequency table of all observed household compositions `(household_size, n_child, n_adult, n_elderly)`.

---

### 2. Hot-deck imputation logic

Three functions work together:

#### `find_donors(obs_child, obs_adult, obs_elderly, hh_size, donor_pool)`
Finds all MICS households that:
- Have the same total `household_size` as the target
- Have **at least as many** children/adults/elderly as already observed in the target household

Assigns sampling probability proportional to frequency in the MICS data.

#### `draw_donor(donors)`
Draws one donor at random, weighted by `freq`. Returns `NULL` if no matching donor exists.

#### `multiple_hot_deck(incomplete_HH, donor_pool)`
Iterates over all incomplete households from the study:
- Calls `find_donors` + `draw_donor` for each
- Computes imputed counts: `imp_child = donor$n_child - obs_child`, etc.
- Handles the no-donor edge case by returning NAs for that household.

---

### 3. Apply to the study data

- Reloads `dataset_freezed2Nov25.csv` and computes observed age-group counts per household.
- Isolates `incomplete_HH` (households with `n_missing_contacts > 0`).
- Runs `multiple_hot_deck(incomplete_HH, comp_table)`.

---

### 4. Validation attempt

- Plots household size distribution for complete vs incomplete households in the study.
- **Finding:** the 34 complete households are mostly size 3 — not a representative validation set for the incomplete ones (which tend to be larger). Validation is therefore limited.

---

## Infection imputation (open section)

The script ends with an unfinished section on imputing **infection status** for missing members. The key insight noted:

> This is not a standard imputation problem — it is a probabilistic inference problem on a transmission network.

One proposed approach: build an **exposure score** based on the infection state of observed household members and age group, possibly using a contact matrix by age.

---

## Relationship to `data_exploration.qmd`

| Aspect | `data_exploration.qmd` | `mics_data.qmd` |
|---|---|---|
| Age imputation method | Random Forest trained on complete study HH | Hot-deck from MICS donor pool |
| External data used | No | Yes (Thailand MICS6 2022, Chiang Mai) |
| Infection imputation | Random Forest + MICE | Not yet implemented |
| Validation | OOB error on RF | Weak (complete HH all size 3) |

Both scripts are exploring different strategies for the same problem. The hot-deck approach is more transparent and does not rely on model assumptions, but depends on MICS being representative of the study population.
