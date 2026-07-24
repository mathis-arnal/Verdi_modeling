# VERDI Modeling Project

## Background: the published study

This project is a methodological follow-up to a published manuscript:

> **"Household SARS-CoV-2 Transmission during Omicron wave in Chiang Mai, Thailand"**
> Khamduang et al., for the Thailand VERDI-RECOVER study team.
> Full text: `bib/verdi_cnx_2026.pdf`

For annotated summaries of all papers in `bib/`, see [bib.md](bib.md).

**Study design:** Prospective observational household transmission study. Enrolled households with ≥2 members including ≥1 child, with a confirmed COVID-19 index case. July 2022 – May 2024, Chiang Mai, Thailand.

**Recruitment and enrollment timing:** Index cases were identified through the Chiang Mai Provincial Public Health Office or primary healthcare hospitals, and referred to the study team. Eligibility required a positive SARS-CoV-2 RT-PCR or antigen result **within the previous 48 hours**. This directly constrains the index's pre-enrollment infectious period: the time from positive test to enrollment is at most 2 days (τ_enroll ≤ 2). Combined with the Omicron pre-test infectious period (~0–1 days before symptom onset/testing), the total pre-enrollment infectious period τ ∈ {1, 2, 3} days, with τ = 2 as the primary estimate. Values of τ = 0 (biologically impossible — index must have been infectious before testing positive) and τ ≥ 4 (ruled out by the 48 h criterion) are excluded.

**No quarantine:** Households were enrolled during the Omicron wave when Thailand had lifted mandatory isolation/quarantine. Household members continued normal daily activities (work, school, etc.) during follow-up. This has two important consequences:
- Contact patterns within the household reflect **normal behaviour**, not forced home confinement — the Thai contact matrix is therefore more applicable as a prior for transmission modelling.
- Community acquisition is a **real competing risk**: a contact testing positive during follow-up may have been infected outside the household, not from the index case. The phylogenetic analysis in the paper confirmed this (12 community transmissions in 10 households). Any transmission model must include a background community hazard term.

**Published results (main analysis):**
- 93 households, 93 index cases (48 children, 45 adults), 197 household contacts
- Individual-level SAR: **23.5%** (44/197 infected); Household-level SAR: **33.3%**
- Key protective factors: baseline IgG positivity (aRR 0.45, −55% risk), lower index viral load (N gene Ct, aRR 0.84 in phylogenetically confirmed sensitivity analysis)
- No significant SAR difference by index case age (child vs adult)
- Statistical method used: GEE with log link, exchangeable correlation structure

**Acknowledged limitation (paper p.10):** "About one third of household members did not consent to participate." 117 contacts declined (mainly fear of venipuncture) + 18 excluded for PCR positivity at enrollment. A sub-analysis comparing complete (35 HH) vs incomplete (58 HH) households found no significant SAR difference (28.6% vs 36.2%, p=0.45), but this is a crude check, not a principled correction.

**This modeling project's goals:**
1. Provide a more rigorous correction for the missing members, producing a revised SAR estimate with proper uncertainty quantification.
2. Quantify the bias introduced by partial household observation — ideally via a simulation study comparing SAR estimates under full vs. partial observation (observed network vs. partially observed), following a suggestion from the Belgium collaborators.
3. Assess the additional bias introduced by potential **multiple index cases / co-primary infections** — contacts who tested positive during follow-up may have been independently exposed from the community, not from the enrolled index case.

---

## What this project is

Epidemiological modeling to correct for **missing household members** in the VERDI-RECOVER cohort. The primary deliverable is a corrected **Secondary Attack Rate (SAR)** that accounts for contacts who refused to enroll, expressed with a posterior credible interval.

**Exclusion mechanism:** members refused consent (fear of venipuncture / not present) — likely **MCAR**, not related to infection status. This justifies treating missing infection status as missing at random conditional on age and household exposure.

The goal of infection imputation is to **correct the SAR denominator and numerator** — estimate how many excluded members were infected — not to produce individual-level predictions.

**Multiple index case bias (raised by Belgium collaborators):** The model assumes the enrolled index case is the sole primary infector. However, some contacts may have been co-primary cases infected from the community simultaneously and independently. The phylogenetic analysis confirmed 12 such community acquisitions. The background hazard term P_BG partially accounts for this, but if a contact was a true co-primary case, attributing their infection to household transmission inflates q and consequently overestimates P(Y_hk = 1) for missing members. This should be addressed as a sensitivity analysis.

**Community hazard calibration:** P_BG = 0.003/day is derived from the study's own phylogenetic data (12/40 infected contacts, ~21-day window) — more reliable than Thai provincial surveillance data for this period, since Thailand moved to endemic management in July 2022 (exactly when enrollment started), making official case counts severe undercounts with no stable correction factor. Provincial situation reports can inform the *shape* of the incidence curve (peak timing) for a time-varying sensitivity analysis, but not absolute levels.

---

## Data

| File | Description |
|---|---|
| `data/dataset_freezed2Nov25.csv` | Main study dataset. One row per household member. Frozen 2 Nov 2025. |
| `Thailand MICS6 2022 SPSS Datasets/hl.sav` | UNICEF MICS6 household-member file (Thailand 2022). Used as donor pool for hot-deck imputation of age composition. |
| `Thailand MICS6 2022 SPSS Datasets/hh.sav` | UNICEF MICS6 household-level file. |
| `TH_contact_data.txt` | Thai social contact survey (source unknown — 257 participants, ages 14–52 only, no elderly). **Not used in the model** — superseded by Prem et al. 2021. |
| `data/contact_home.rdata` | Prem et al. 2021 synthetic home contact matrix for Thailand (auto-downloaded by `chain_binomial_mcmc.qmd` on first run). Named list of 177 16×16 matrices keyed by ISO3; Thailand = `contact_home[["THA"]]`. **Not used by the primary model** (`1-chain_binomial_simple.qmd`), which is age-homogeneous; only used by the age-stratified extension `chain_binomial_mcmc.qmd`. |

All MICS analyses are filtered to **Chiang Mai** (`HH7A == 50`) to match the study population.

---

## Key variables in the study dataset

### Identifiers and structure

| Variable | Meaning |
|---|---|
| `hhid` | Household identifier |
| `pid` / `index` | Person identifier; `index == 1` marks the index case |
| `exclude_mem` | Number of household members NOT enrolled (propagated to all rows in the HH) |
| `household_size` | Total household size including excluded members |
| `include_all_members` | `"All members included"` for complete households |
| `relation_cat` | Contact's relationship to index: Child/Grandchild (n=40), Father/Mother (65), Grandfather/Grandmother (21), Husband/Wife/Lover (22), Relative (49). Proxy for within-household contact intensity. |
| `sex` | Male / Female — no missing values. Sex-stratified home contact matrix for Thailand is unavailable, so not currently used as a transmission modifier. |

### Infection outcomes

| Variable | Meaning |
|---|---|
| `result_infected` | `"Detected"` / `"Not detected"` — final binary infection outcome over the full follow-up |
| `result_day1` | PCR result on day 1 of study. For contacts: `"Detected"` identifies potential co-primary cases (infected before or at enrollment). Contacts already positive on day 1 were supposed to be excluded; residual positives here are edge cases. |

### Infection timing — critical for chain-binomial

| Variable | Meaning |
|---|---|
| `timetoinf` | Days from enrollment (day 1) to first positive PCR test for **contacts** only. NA for uninfected contacts and for all index cases. Provides interval-censored infection timing (PCR schedule is every other day → true infection day is approximately `timetoinf − 1` to `timetoinf`). Range: 1–12 days. |
| `day1_date` | Enrollment date (format `dd/m/yyyy`). Anchors `timetoinf` to calendar time. |
| `infected_date` | **Empty (NA) for all 93 index cases.** No record of when the index case was first detected or started being infectious before enrollment. |
| `infection_date` | **Empty (NA) for all 93 index cases.** Same as above. |

**Critical interpretation caveat for `timetoinf`:** `timetoinf` is measured from enrollment, not from the start of the index case's infectious period. The 48-hour enrollment criterion constrains the time from positive test to enrollment to ≤ 2 days. Combined with Omicron biology (~0–1 pre-test infectious days), the total pre-enrollment infectious period τ ∈ {1, 2, 3} days. In the model, `TAU = 2` is the primary estimate; sensitivity analysis at τ ∈ {1, 3} spans the plausible range given the enrollment criterion. τ = 0 is biologically implausible; τ ≥ 4 is ruled out by the 48 h criterion. A future enhancement is to treat τ_h as a per-household discrete latent variable with prior Uniform{0, 1, 2} (the test-to-enroll component), updated via an exact Gibbs step enumerating three values.

### Immunological and virological covariates

| Variable | Meaning |
|---|---|
| `result_IgG_day1_mod` | Baseline IgG positivity at enrollment — the strongest protective factor in the published analysis (aRR 0.45, −55% risk). Available for enrolled contacts only; missing for non-enrolled members. |
| `ngene_day1` / `sgene_day1` | N/S gene Ct values at day 1 (index viral load proxy). Lower Ct = higher viral load = higher infectiousness. `ngene_day1` available for 92/93 index cases (range 18.96–34.86, median 25.25); reserved for β_Ct extension, not used in base model. |
| `lineagecat` | Viral lineage category (BA.2, BA.4/5, XBB, etc.). |

### Household-level covariates

| Variable | Meaning |
|---|---|
| `bedroom_num` / `separate_bed` | Number of bedrooms and whether index slept separately — proxy for overnight contact intensity. |
| `index_wear_mask` / `member_wear_mask` | Mask use during follow-up — potential modifier of within-household transmission probability. |
| `age_cat` | Child (≤18) / Adult (18–59) / Elderly (≥60) — used to index the Prem contact matrix. |
| `index_child` | Whether the index case is a child — no significant SAR difference found in published analysis. |

---

## Scripts

| Script | Purpose |
|---|---|
| `scripts/1-chain_binomial_simple.qmd` | **Main model** — age-homogeneous timed chain-binomial with Metropolis-within-Gibbs MCMC; single scalar λ_h (a.k.a. `q` in code); latent Y_hk, T_hk (infection day) for missing members, updated via an **exact joint Gibbs step**; secondary-to-secondary chains; corrected SAR with posterior credible interval. Sources `R/model_functions_simple.R`. Fits on 80 households (93 minus 14 excluded for a contact PCR-positive at day-1 enrollment — ambiguous co-primary status). |
| `scripts/chain_binomial_mcmc.qmd` | **Age-stratified extension** (formerly the main model) — adds the Prem et al. (2021) 3×3 contact matrix `C_mat` and a latent `age_k` for missing members on top of the same timed chain-binomial likelihood. Sources `R/model_functions.R`. Kept as a planned refinement/sensitivity check on the primary model, not the current primary analysis. |
| `scripts/0-preprocess.qmd` | Stub — not yet implemented; intended as the shared preprocessing entry point for the numbered pipeline (`0-preprocess` → `1-chain_binomial_simple` → ...). |
| `scripts/model_presentation.qmd` | Beamer slide deck for collaborator presentations — covers study context, missing data problem, model structure, MCMC sampler, Omicron biology justification, Q&A backup slides |
| `scripts/sensitivity_analysis.qmd` | Sensitivity analyses — τ ∈ {1,3}, density-dependent contact correction γ ∈ {0,0.5,1}, time-varying P_BG. Planned addition: `q` prior sensitivity (`Beta(1,19)` vs vague `Uniform(0,1)`), motivated by comparison with Tsang & Cauchemez (2023)'s prior choices (see below). |
| `scripts/data_exploration.qmd` | Preprocessing, missing data exploration, RF-based imputation (exploratory, superseded by chain-binomial for primary analysis) |
| `scripts/mics_data.qmd` | Hot-deck imputation of excluded members' age composition using MICS6 Chiang Mai as donor pool |
| `scripts/contact_data.qmd` | Loading and exploration of the Thai contact matrix (feeds the age-stratified extension only) |
| `scripts/timed_chain_binomial_explainer.qmd` | Pedagogical worked example of the timed chain-binomial likelihood for collaborators — self-contained, no dependency on main model |

---

## Current modeling approach

**Primary method: age-homogeneous discrete-time timed chain-binomial with Metropolis-within-Gibbs MCMC** (`1-chain_binomial_simple.qmd`, sourcing `R/model_functions_simple.R`). No contact matrix, no age stratification — all contacts treated as homogeneous. The age-stratified version (contact matrix + `age_k` latent variable) that used to be primary now lives in `chain_binomial_mcmc.qmd` as a planned extension/sensitivity check, not the current primary analysis.

**Sample**: fit on **80 households** (of the published 93) — 14 households (15 contacts) with a contact already PCR-positive at day-1 enrollment are excluded outright, since it's then impossible to attribute their infection to the index case vs. the co-positive contact. This is a stricter exclusion than relying on `P_BG` alone to absorb that ambiguity.

### Methodological precedent chain

| Step | Our model | Citable precedent |
|---|---|---|
| Likelihood form | Additive daily hazards, exp(−λ) escape, E_i(t_i−1)−E_i(t_i) for timed infections | **Tsang & Cauchemez (2023, eLife)** `bib/tsang_cauchemez_2023.pdf` — shared likelihood structure (additive hazards, escape probability). Precedent for the *likelihood structure only* — see prior row below for why it is **not** precedent for our prior choice. |
| Infectiousness profile f(s) | Discretized Gamma serial interval, shape = 2.11, rate = 0.50 (mean ≈ 4.2 days) — replaces an earlier flat rectangular D=5-day window | **an der Heiden et al. (2022)** — Omicron-period household serial interval fitted from 11,512 German household clusters. Not yet added to `bib/`; cite directly pending PDF acquisition. Operational proxy: since `timetoinf` is the ATK detection day (not the transmission day), the distribution of time between consecutive positive tests is treated as an operational proxy for the serial interval, following **Cauchemez et al. (2009)**'s use of the serial interval in the hazard-of-infection model. |
| Discrete-time daily hazard (vs continuous) | Product of daily escape probabilities | **Cauchemez et al. (2004)** `bib/cauchemez2004.pdf` — fixed infectious period per infector |
| MCMC with infection time augmentation | T_hk, Y_hk as latent variables, updated **jointly** via exact (T_follow+1)-way categorical Gibbs (see §6 below) | **Cauchemez et al. (2004)** — augments ν_i, ψ_i (start/end of infectious period); **Ganyani et al. (2020)** — two-step MCMC confirming the augmentation pattern. The joint-categorical implementation (rather than updating Y_hk against a single candidate T) is our own correctness fix: Y_hk's true full conditional marginalizes over all T_hk values, so comparing "not infected" against one candidate day (instead of the sum over all days) understates the infected branch's posterior mass whenever it's spread across several plausible days. |
| Homogeneous mixing (β=1, δ_i=1; no household-size dilution, no individual infectiousness heterogeneity) | **Deliberate simplification** of Tsang & Cauchemez's fuller model (which has X_k^β dilution and δ_i ~ N(0,σ_var) heterogeneity) | Motivated by identifiability with only 44 secondary cases across 80 households; age stratification (contact matrix C[a_i,a_j]) and household-size dilution are reserved for the `chain_binomial_mcmc.qmd` extension and the γ ∈ {0,0.5,1} sensitivity analysis respectively |
| Missing member imputation | Y_hk, T_hk as discrete latent variables, exact joint Gibbs | **Unique contribution** — no prior COVID-19 household study has formally treated non-enrolled members as latent variables in a transmission model |
| Community hazard | Fixed P_BG = 0.003/day | **Khamduang et al. (2026)** — phylogenetic analysis: 12/40 community-acquired contacts over ~21 days |
| Prior on λ_h (a.k.a. `q`) | Informative `Beta(1,19)` (mean 0.05, SD ≈ 0.047) | **Our own choice**, not inherited from Tsang & Cauchemez — they used **vague** priors (`Uniform(0,10)` for their unbounded hazard-rate parameter λ_k; `Uniform(0,1)` for their bounded-[0,1] dilution exponent β). Our λ_h is bounded to (0,1) like their β, so `Uniform(0,1)` is the closer "vague" analogue if one were wanted — but with only 44 secondary cases plus a large latent state space (Y_hk, T_hk), a flat prior gives materially less regularization than `Beta(1,19)`, and `Uniform(0,1)` is not actually weakly-informative on a per-contact-per-day probability scale (it puts real mass on epidemiologically implausible values like 0.5–0.9). Decision: keep `Beta(1,19)` as primary, add `Uniform(0,1)` as a **planned sensitivity analysis** in `sensitivity_analysis.qmd` to show the result isn't prior-driven. |

### Architecture — script sections (`1-chain_binomial_simple.qmd` / `model_functions_simple.R`)

**§1 — Household data** (`build_hh_list_simple`)
- Reads `data/dataset_freezed2Nov25.csv` (**Khamduang et al. 2026**), filtered to the 80 households without a day-1-positive contact. Extracts per-household: `tau` (pre-enrollment infectious period, primary τ = 2), enrolled contacts with `infected_bin`, `timetoinf`, and `n_missing` (= `exclude_mem`). No `age_cat` — age-homogeneous model.
- τ ∈ {1,2,3} is constrained by the 48-hour enrollment criterion (see CLAUDE.md Background).

**§2 — Likelihood** (`compute_escape_seq_simple`, `hh_loglik_timed_simple`)
- Daily force of infection on susceptible i at day t: `λ_i(t) = −log(1−P_BG) + q·f(τ+t) + Σ[j active on t: q·f(t−t_j+1)]`, where `f` is the discretized Gamma serial-interval profile.
- Cumulative escape: `E_i(t) = exp(−Σ_{s=1}^{t} λ_i(s))`.
- Likelihood: `E_i(t_i−1) − E_i(t_i)` if infected on day t_i (from `timetoinf`); `E_i(T_follow)` if uninfected; `1 − E_i(T_follow)` fallback if infected but `timetoinf` missing.
- Secondary-to-secondary chains fully modelled: observed infected contacts and infected missing members both enter the infector list from their respective infection days (**Cauchemez et al. 2009** — tertiary infections allowed).
- Single scalar `q` (λ_h) — no per-household modifier. β_Ct extension planned.

**§3 — Missing member prior** (`prior_miss_simple`)
- Uses index + community hazard + already-enrolled infectors (fixed/observed, so no circular dependence) to avoid circularity in the Gibbs update.
- `P(Y_hk=0) = E_k^prior(T_follow)`; `P(Y_hk=1, T_hk=t) ∝ E_k^prior(t−1) − E_k^prior(t)`.
- Also used to initialise Y_hk^(0)/T_hk^(0) at the start of the chain (not a flat `Uniform{1..14}` or `Bernoulli(p)` — draws from this same model-implied distribution).
- Data augmentation pattern: **Cauchemez et al. (2004)**, **Ganyani et al. (2020)**.

**§4 — MCMC sampler** (`run_mcmc_simple`): two update components per iteration
1. **MH for `q`** — log-scale random walk (`Q_SD` proposal SD); `Beta(1,19)` prior. `q_init`/`q_sd` in the qmd are MCMC mechanics (starting value, proposal step size), **not** prior parameters — the prior itself (`alpha_q=1, beta_q=19`) is set inside `run_mcmc_simple`.
2. **Joint Gibbs for (Y_hk, T_hk)** — exact `(T_follow+1)`-way categorical per missing member: state 1 = "not infected", state t+1 = "infected on day t" for t = 1..T_follow. Each state scored by `hh_loglik_timed_simple` + `prior_miss_simple` log-weight, normalised via log-sum-exp, and **sampled** (not argmax — taking the mode would collapse the chain to ICM/coordinate-ascent MAP estimation and destroy the posterior uncertainty that the corrected-SAR credible interval depends on).

**§5 — MCMC run**: single chain, `N_ITER=20000`, `N_BURNIN=5000`, `N_THIN=10` → 1,500 posterior draws. Reports MH acceptance rate for `q` (target 20–40%).

**§6 — Results**: posterior of λ_h/`q` (median + 95% CrI); corrected SAR = (n_obs_infected + E[Σ Y_hk]) / (n_obs_contacts + n_total_missing) with 95% CrI; expected infections among missing members.

### Age-stratified extension (`chain_binomial_mcmc.qmd` / `model_functions.R`)

Same likelihood/augmentation framework as above, plus:
- **Contact matrix** `C_mat` (3×3 Child/Adult/Elderly) from **Prem et al. (2021)**, home-setting matrix for Thailand (THA), aggregated from 16×16 to 3 groups; downloaded on first run, cached to `data/contact_home.rdata`. Used as a fixed point estimate (Prem's 200 bootstrap draws per country not propagated — **Zhang et al. (2026)** §3.3).
- Fixing `C` and estimating only scalar `q` is a **necessary identifiability constraint**: `q·C[a_i,a_j]` is structurally non-identifiable (any (q,C) pair with the same product is observationally equivalent) — **Zhang et al. (2026)** §4.2. Precedent for age-stratified rates: **Cauchemez et al. (2009)** `bib/cauchemez2009.pdf`.
- Adds a fourth Gibbs step: **age_k** — exact three-way (Child/Adult/Elderly), joint posterior with likelihood + `prior_miss()` + MICS6 marginal (**UNICEF MICS6 2022**, `compute_mics_priors()`).
- Ct infectiousness modifier: **not implemented**, planned extension motivated by **Khamduang et al. (2026)** aRR 0.84 per Ct unit and **Li et al. (2026)** `bib/hhbayes_2026.pdf` (HHBayes).

**Exploratory methods (superseded for primary analysis):**
- Random Forest (`data_exploration.qmd`): predicts age and infection for excluded members; over-predicts "Adult"; OOB error evaluated on complete households only — not representative.
- MICE (`data_exploration.qmd`): PMM is wrong for binary outcomes; `logreg` should be used instead.

---

## Key model assumptions and known limitations

| Assumption | Status | Impact if violated |
|---|---|---|
| TAU = 2 days fixed across all households | Primary; sensitivity at τ ∈ {1,3} | True τ_h varies per household (0–2 days test-to-enroll); augmenting τ_h as discrete latent Uniform{0,1,2} would propagate this uncertainty |
| β_Ct viral load modifier | **Not in base model — planned extension** | Single q absorbs infectiousness heterogeneity; adding β_Ct requires identifiability check given n=44 secondary cases |
| Single transmission parameter q (no IgG modifier) | **Deliberate simplification** | No reason to believe IgG prevalence differs between enrolled and missing members (MCAR via venipuncture refusal); under MCAR the correction is unnecessary. IgG not observed for missing members anyway. |
| Missing members assumed IgG-naive | Unavoidable (no blood draw) | Conservative upward bias on corrected SAR; magnitude depends on true IgG prevalence in missing members |
| Enrolled index case is sole primary infector | **Known gap — sensitivity analysis needed** | Co-primary community infections (n=12 phylogenetically confirmed) inflate q if misattributed; partially addressed by P_BG |
| Secondary-to-secondary transmission | **Implemented** via timetoinf daily hazard model | HH 42 (timetoinf=7,10) now correctly attributed to transmission chains beyond index infectious window |
| Fixed P_BG = 0.003/day (not time-varying) | Reasonable as primary; sensitivity warranted | Community hazard varied across July 2022 – May 2024 sub-waves; provincial surveillance data unreliable post-endemic declaration |
| Age of missing members from MICS hot-deck | Well-motivated for the age-stratified extension; **not used in the primary model** (age-homogeneous) | MICS is Thailand 2022, filtered to Chiang Mai; minor if composition differs |
| Infectiousness profile f(s) is a Gamma-shaped serial interval, not a flat window | Not a rectangular D=5-day window in either model (that description was stale) — the profile itself is still not augmented per-infector | Cauchemez (2004) augments ψ_i−ν_i (duration) as a continuous latent variable per infector; we instead fix one population-level Gamma shape for all infectors, absorbing duration/infectiousness heterogeneity into q. Primary model uses a serial-interval parameterization (an der Heiden et al. 2022, shape≈2.11); the age-stratified extension (`chain_binomial_mcmc.qmd`) still uses a generation-time parameterization (shape≈4.84) — the two scripts are not yet reconciled to the same f(s), worth aligning |
| Homogeneous mixing — no age stratification, no contact matrix (primary model) | **Deliberate primary-model simplification**, not a limitation of the framework itself | Age-stratified transmission (contact matrix C[a_i,a_j], latent age_k) is implemented in `chain_binomial_mcmc.qmd`; dropped from the primary model to reduce latent-state complexity relative to 44 secondary cases across 80 households |
| Excluding 14 households with a contact PCR-positive at day-1 enrollment | **Primary-model choice** — stricter than the previous approach of retaining all 93 households and relying on P_BG to absorb the ambiguity | Reduces the analysed sample from 93 to 80 households; cleaner causal attribution (no household with a plausible alternative infector already positive at enrollment) at the cost of a smaller, and slightly different, sample than the published GEE analysis |
| Prior on λ_h/`q`: informative `Beta(1,19)` rather than a vague prior | Deliberate — see precedent-chain table above | With sparse data (44 secondary cases) a flat `Uniform(0,1)` prior gives materially less regularization and isn't truly weakly-informative on a probability scale; `Beta(1,19)` encodes real epidemiological plausibility. Planned `Uniform(0,1)` sensitivity run will quantify how much the informative prior actually drives the posterior |
| MCAR missingness | Assumed; crude sub-analysis in paper supports it | If IgG-positive members more likely to refuse, corrected SAR further overestimated |
| Prem contact matrix C[a_i, a_j] is a population average; model is density-dependent (no household-size correction) | **Known structural limitation** | The Prem matrix gives mean daily home contacts between age groups across the Thai population; the hazard q_h · C[a_i, a_j] does not decrease with household size (density-dependent). Goeyvaerts et al. (2018) `bib/hhmem_nmar.pdf` showed empirically (318 Belgian households) that within-household contact density decreases with size (≈1.00 for n=2, ≈0.85 for n≥6) and that ignoring this biases transmission rate estimates. They do not provide a plug-in correction for Thailand, but frame the problem as mean contact degree ∝ n^w with w ∈ (0,1); w=0 is our current model, w=1 is Cauchemez's frequency-dependent β/n. Planned sensitivity: re-run at effective hazard q_h · C[a_i, a_j] / n_h^γ for γ ∈ {0, 0.5, 1}; practical impact expected modest for n=3–5 households (density difference ~6 pp per Goeyvaerts Table 1). |
| Prem contact matrix used as a single point estimate (sampling uncertainty not propagated) | **Known gap** — Zhang et al. (2026) §3.3 | Prem provides 200 bootstrap draws per country stored separately on GitHub (not in `contact_home.rdata`). Running the MCMC on a subsample of draws would propagate matrix sampling uncertainty. Practical impact expected modest: Thai demographic projection is well-constrained; Prem's posterior over C[THA] is narrow. Mitigating factor: no-quarantine study design makes the pre-pandemic Prem matrix more applicable here than for lockdown studies (Zhang §4.1). |

---

## R environment

R project file: `VERDI_modeling.Rproj`. Key packages: `dplyr`, `stringr`, `naniar`, `randomForest`, `mice`, `ggplot2`, `caret`, `haven`.
