# Annotated Bibliography — VERDI Modeling Project

Papers are organized thematically by their primary contribution to the project. Each entry includes full citation, core finding, and explicit notes on how the paper relates to our SARS-CoV-2 household transmission model correcting for missing members.

---

## 1. Primary Data Source

### Khamduang et al. (2026) — `verdi_cnx_2026.pdf`
**"Household SARS-CoV-2 Transmission during Omicron Wave in Chiang Mai, Thailand: a prospective observational study"**
*Lancet Regional Health — Southeast Asia*, 44: 100711.

**Core findings:** 93 households (48 child index, 45 adult index), 197 enrolled contacts, Chiang Mai, July 2022–May 2024. Individual-level SAR 24% (95% CI 17–32%); household-level SAR 33% (95% CI 24–44%). Protective factors: baseline anti-NCP IgG positivity (aRR 0.42, 95% CI 0.22–0.83) and lower index N-gene Ct (aRR 0.82, 95% CI 0.74–0.92). No significant SAR difference by index age. Phylogenetic analysis confirmed 28 intra-household transmissions in 23 households and **12 community acquisitions** in 10 households. 117 contacts declined participation (venipuncture/swabbing fear); 18 excluded for PCR positivity at enrollment.

**Relevance:** This is the dataset being modeled. Key quantities used directly:
- `n_missing`, `exclude_mem`, `household_size` → drives the imputation problem
- `result_infected`, `timetoinf` → likelihood inputs
- `result_IgG_day1_mod` → observed covariate for enrolled members (not modeled for missing)
- `ngene_day1` → Ct values for β_Ct extension
- The 12 community acquisitions calibrate `P_BG = 0.003/day` (~12/40 infected × 1/21 days)
- The crude sub-analysis (complete vs. incomplete households: 28.6% vs. 36.2%, p = 0.45) motivates but does not replace our rigorous correction

---

## 2. Methodological Backbone — MCMC for Household Transmission

### Cauchemez et al. (2004) — `cauchemez2004.pdf`
**"A Bayesian MCMC approach to study transmission of influenza: application to household longitudinal data"**
*Statistics in Medicine*, 23: 3469–3487.

**Core findings:** 334 households, France, 1999–2000 influenza season. Three-level hierarchical model: observation level (compatibility of augmented and observed data), transmission level (instantaneous hazard), and prior level. Augments unobserved start (ν_i) and end (ψ_i) of each infector's infectious period as continuous latent variables. Community hazard α_s; household hazard ε_s β_i / n (frequency-dependent: β/n). Infectious duration ~ Gamma(μ_i, σ_i). Mean infectious period 3.8 days (95% CrI 3.1–4.6). Children more likely to transmit than adults (posterior probability >99%).

**Relevance — structural backbone of our model:**
- The three-level hierarchical likelihood structure (observation ↔ transmission ↔ prior) is the direct ancestor of our `hh_loglik_timed`
- The data augmentation pattern — treating unobserved infection timing as latent variables updated by MCMC — is exactly what we do for Y_hk, T_hk (missing members) and T_hk (observed infected contacts)
- The β/n frequency-dependent household hazard informs our density-dependent sensitivity analysis (γ = 1 case)
- Their MH scheme for β + Gibbs for ν, ψ parallels our MH for q + Gibbs for Y_hk, T_hk, age_k
- **Key difference:** they augment the continuous start/end of each infector's infectious period; we fix D = 5 days and instead augment discrete infection days for missing members

### O'Neill & Roberts (1999) — `oneill1999.pdf`
**"Bayesian inference for partially observed stochastic epidemics"**
*Journal of the Royal Statistical Society Series A*, 162(1): 121–129.

**Core findings:** Foundational paper for Bayesian MCMC on partially observed epidemics. Two models: (i) Reed–Frost (discrete-time, avoidance probability q) — Gibbs sampler, augments the latent epidemic pattern vector n_{21}; (ii) General stochastic SIR (continuous-time) — MH with three moves (move/add/remove infection times). Applied to Providence measles data (households of size 3). Posterior for q robust to prior specification.

**Relevance:**
- Establishes the conceptual framework underlying all subsequent household MCMC work in our bibliography
- The Reed–Frost formulation (avoidance probability q, escape structure q^{I_t}) is the discrete-time ancestor of our per-day escape probability
- Demonstrates that augmenting latent infection counts/times is computationally feasible and yields valid posteriors
- The three-move MH algorithm for infection times (move/add/remove) is the continuous-time precursor to our discrete Gibbs step for T_hk

### Cauchemez et al. (2009) — `cauchemez2009.pdf`
**"Household Transmission of 2009 Pandemic Influenza A (H1N1) Virus in the United States"**
*New England Journal of Medicine*, 361: 2619–2627.

**Core findings:** 216 complete households, 600 contacts, US H1N1 pandemic, spring 2009. Bayesian MCMC with data augmentation. Age-stratified susceptibility: children ≤18 twice as susceptible as contacts 19–50 (relative susceptibility 1.96, 95% CrI 1.05–3.78); contacts >50 substantially less susceptible (0.17). Infectivity did not vary with age. Transmission decreases with household size. Mean serial interval 2.6 days. Tertiary infections allowed.

**Relevance:**
- Directly justifies our age-stratified hazard via C[a_i, a_j]
- Confirms that tertiary infections should be modelled (secondary-to-secondary chains in our §4)
- The household-size effect supports our density-dependent sensitivity analysis
- The Bayesian MCMC + augmentation approach is the same class as ours; their supplement details the likelihood form closest to ours

### Walker et al. (2017) — `walker2017.pdf`
**"Inference of epidemiological parameters from household stratified data"**
*PLoS ONE*, 12(10): e0185910.

**Core findings:** Continuous-time SIR Markov chain with two mixing levels (within-household rate β, between-household rate α, recovery rate γ). Two Bayesian MCMC methods for FF100-type data: (i) exact DA-MCMC augmenting infection and recovery times with five move types; (ii) branching process approximation that is more scalable with more households but introduces slight positive bias in α. Both methods recover true parameters well in simulations.

**Relevance:**
- Confirms the general class of DA-MCMC for household models with missing infection times
- The branching process approximation is potentially relevant if we need to scale up to larger datasets or add a between-household mixing term
- The five-move DA-MCMC structure (including adding/removing recovery events) contrasts with our simpler discrete-time Gibbs, useful for justifying our design choices in the methods section

---

## 3. Primary Likelihood Reference

### Tsang & Cauchemez (2023) — `tsang_cauchemez_2023.pdf`
**"The effect of variation of individual infectiousness on SARS-CoV-2 transmission in households"**
*eLife*, 12: e82611.

**Core findings:** Meta-analysis of 17 SARS-CoV-2 household studies (13,098 index cases, 31,359 contacts). Individual-based model with additive hazard: λ_{i→j}(t) = (λ_k / X_k^β) × exp(δ_i) × f(t − t_i), where δ_i ~ N(0, σ_var) is individual infectiousness random effect, X_k^β is household-size dilution, and f(·) is an infectiousness profile. Bayesian DA-MCMC augmenting missing infection times. Pooled σ_var = 1.33; top 20% of cases are 3.1× more infectious than average. Community infections and tertiary infections explicitly modelled.

**Relevance — this is the primary likelihood precedent for our model:**
- The additive hazard structure λ_{i→j}(t) = λ_community + Σ_infectors λ_{i→j}(t) is directly what we implement in `compute_escape_seq`
- The escape probability product exp(−Σ_t λ(t)) and the likelihood contribution E(t−1) − E(t) for timed infections is identical to our formulation
- The data augmentation of missing infection times via MCMC is the same pattern as our T_hk Gibbs step
- The household-size dilution parameter β (their X_k^β) motivates our γ ∈ {0, 0.5, 1} sensitivity analysis
- σ_var (infectiousness heterogeneity) is the parameter we deliberately excluded from the base model for identifiability; this paper shows how to add it as a future extension

---

## 4. Contact Matrix

### Prem et al. (2021) — `prem2021.pdf`
**"Projecting contact matrices in 177 geographical regions: An update and comparison with empirical data for the COVID-19 era"**
*PLoS Computational Biology*, 17(7): e1009098.

**Core findings:** Updated synthetic contact matrices (home, work, school, other) for 177 geographical regions covering 97.2% of world population. Bayesian hierarchical model using POLYMOD European empirical data combined with country-specific household age structure (DHS), labour force participation (ILO), and school enrollment (UNESCO). Rural/urban stratification for 43 countries. Synthetic matrices show good qualitative agreement with out-of-sample empirical matrices. Corrected in 2024 (see below).

**Relevance:**
- Provides `contact_home[["THA"]]` — the 16×16 Thailand home contact matrix aggregated to our 3×3 Child/Adult/Elderly matrix C[a_i, a_j]
- The matrix directly encodes age-assortative contact patterns that modulate the hazard λ_i(t) in `hh_loglik_timed`
- Thailand has no empirical household contact survey for this age range, making the synthetic matrix the only available option
- The 2024 correction (`correc_prem2024.pdf`) fixes coding errors in HAM construction and visitor extrapolation; qualitative conclusions unchanged but updated matrices should be used

### Prem et al. (2024) — `correc_prem2024.pdf`
**"Correction: Projecting contact matrices in 177 geographical regions..."**
*PLoS Computational Biology*, 20(9): e1012454.

**Core findings:** Correction notice. Fixed errors: (i) HAM subscript error for countries without household age data; (ii) missing visitor parameter δ^H in household contact extrapolation; (iii) age-specific population ratio applied to wrong matrix dimension; (iv) frequency-dependent adjustment removed for workplace/school. Corrected code at GitHub (kieshaprem/synthetic-contact-matrices). Results quantitatively slightly improved; main conclusions unchanged.

**Relevance:** The matrices downloaded by `chain_binomial_mcmc.qmd` should come from the corrected repository. If `contact_home[["THA"]]` was cached before September 2024, verify against the corrected version.

### Goeyvaerts et al. (2018) — `hhmem_nmar2018.pdf`
**"Household members do not contact each other at random: implications for infectious disease modelling"**
*Proceedings of the Royal Society B*, 285: 20182201.

**Core findings:** First social contact survey specifically designed for within-household contact networks. 318 Belgian households (2–7 members), 1266 participants, Flanders and Brussels 2010–2011. Exponential-family random graph models (ERGMs). Key results: (i) within-household contact networks are almost always complete (density ≈ 0.93); (ii) contact density **decreases with household size** on weekdays (proportion of complete networks: 1.00 for n=2, 0.77 for n=4, 0.46 for n≥6); (iii) homogeneous random mixing is an adequate structural approximation for epidemic simulations; (iv) however, ignoring contact density biases upward estimates of the within-household transmission rate. The scaling exponent w ∈ (0, 1) where mean contact degree ∝ n^w (w=0 = density-dependent/our current model; w=1 = frequency-dependent β/n).

**Relevance:**
- Directly cited in our CLAUDE.md as the empirical basis for the density-dependent sensitivity analysis
- Motivates re-running our model with effective hazard q × C[a_i, a_j] / n_h^γ for γ ∈ {0, 0.5, 1}
- Practical impact expected modest for n = 3–5 households (density difference ~6 pp per their Table 1)
- Confirms that using Prem et al.'s population-average C[a_i, a_j] without household-size correction is a known bias direction (upward on q)

### Zhang et al. (2026) — `zhang2026.pdf`
**"Constructing Contact and Connectivity Matrices for Infectious Disease Modelling"**
*arXiv:2605.30034* (preprint, not peer-reviewed)

**Core findings:** Comprehensive review of data types, uncertainty methods, and identifiability issues for contact/connectivity matrices in epidemic modelling. Distinguishes *aleatory* uncertainty (irreducible stochasticity in individual contact behaviour) from *epistemic* uncertainty (incomplete knowledge, finite survey samples, projection assumptions). Organises methods along two axes (Figure 4): model stratification (national → local + age) and parameter inference complexity (fixed matrix → joint inference of matrix and epidemic state). Calls out that structural identifiability is widely ignored: the standard parameterisation β_ij = p · c_ij means p and c_ij are jointly non-identifiable without external constraints. Recommends propagating sampling uncertainty by running the epidemic model on multiple draws from the matrix posterior, though notes this is rarely done in practice. For historical (pre-epidemic) matrices applied to epidemic data, warns that mismatches in behaviour can bias age-specific attack rates.

**Relevance:**
- **Validates our fixed-matrix + scalar q approach (§3.2):** using Prem et al. as a fixed input and estimating a scalar transmission probability on top is the standard practice, classified by Zhang as "Fixed/Scaled connectivity matrix." Cite §3.2 in the methods section to justify the design
- **Formalises the structural non-identifiability justifying single-q (§4.2):** β_ij = q · c_{ij} is non-identifiable — any (q, c_{ij}) pair with the same product is observationally equivalent. Fixing C from Prem and estimating only q is therefore not just convenient but *necessary* for identifiability with n = 44 secondary cases. Zhang §4.2 provides the theoretical language to state this explicitly in the manuscript
- **Identifies a gap: matrix sampling uncertainty is not propagated (§3.3):** our model uses a single point-estimate `contact_home[["THA"]]`. Zhang recommends running the model on multiple draws from the matrix posterior. Prem's bootstrap samples (200 per country) are available on the GitHub repository in separate files but are not in `contact_home.rdata`. Practical impact is expected to be modest because Prem's demographic projection for Thailand is well-constrained and matrix entries do not vary greatly across draws. Should be acknowledged as a limitation
- **Historical matrix vs. epidemic-period behaviour (§4.1):** Prem is a pre-pandemic projection; Zhang warns this can bias age-specific attack rate estimates. For VERDI this concern is substantially mitigated by study design — no mandatory quarantine, normal household behaviour during the Omicron wave — so the Prem home matrix, which reflects normal (non-epidemic) household contact patterns, is more applicable here than for lockdown studies. This mitigating factor should be stated explicitly
- **Figure 4 taxonomy:** our model sits in the "Fixed pre-epidemic age-contact structure, national/local non-spatial" cell — confirms we are in good company with Prem et al. [74] themselves and Cauchemez et al.
- **Not peer-reviewed** — cite as preprint

---

## 5. SAR Estimation — Biases and Alternatives

### Sharker & Kenah (2021) — `sharker2021.pdf`
**"Estimating and interpreting secondary attack risk: Binomial considered biased"**
*PLoS Computational Biology*, 17(1): e1008601.

**Core findings:** Mathematical proof using probability generating functions that binomial models (logistic regression, GEE) produce **upward-biased SAR estimates** even for small true SARs, because multiple generations of transmission inflate the final attack rate (FAR) above the SAR. The FAR/SAR divergence increases with household size and SAR. Cluster-adjusted variances (GEE) correct coverage probabilities only slightly and do not fix the point estimate bias. Longitudinal chain binomial models and pairwise survival analysis produce unbiased estimates. Applied to LA County H1N1 2009 data: binomial SAR ≈ 17–20% vs. chain binomial SAR ≈ 11–13%.

**Relevance — critical methodological justification:**
- The VERDI paper used GEE with log link and exchangeable correlation to estimate SAR = 23.5% — this is precisely the estimator Sharker & Kenah prove is biased upward
- Our chain-binomial MCMC model for q is methodologically superior and should yield a lower estimate than 23.5%
- Provides language to explain why our corrected estimate may differ from the published figure
- The FAR vs. SAR distinction maps directly onto our secondary-to-secondary chain modelling (§4)

### Lindstrøm et al. (2024) — `hh_incomplete_chain_bin_2024.pdf`
**"Estimating the household secondary attack rate with the Incomplete Chain Binomial model"**
*arXiv:2403.03948*

**Core findings:** Derives the Incomplete Chain Binomial (ICB) distribution for outbreaks not yet concluded at observation time. Shows that using the final-size distribution when the outbreak is still ongoing systematically underestimates the SAR (bias increases with larger households and shorter follow-up). Provides an R package `chainbinomial` implementing MLE for the SAR (and SAR-as-GLM) using both final and incomplete distributions. Simulation study shows that confidence intervals based on Wilks' theorem outperform normal-approximation CIs.

**Relevance:**
- Complements Sharker & Kenah by addressing incomplete follow-up bias (in addition to multiple-generation bias)
- Our 14-day follow-up for contacts and 21-day symptom monitoring is sufficient for Omicron (generation time ~3 days), so incomplete-observation bias is likely small, but this paper provides the formal justification
- The GLM-like framework (SAR as function of predictors) is an alternative to our Bayesian approach worth citing when discussing inferential choices
- The `chainbinomial` R package could be used for a quick sensitivity check of our posterior q against a frequentist chain-binomial MLE

---

## 6. Model Extensions and Related Methods

### Cauchemez et al. (2014) — `cauchemez2014.pdf`
**"Determinants of Influenza Transmission in South East Asia: Insights from a Household Cohort Study in Vietnam"**
*PLoS Pathogens*, 10(8): e1004310.

**Core findings:** 940 participants, 270 households, Ha Nam province Vietnam, 2007–2010. Seasonal (H1N1, H3N2, B) and pandemic H1N1. Bayesian digraph data-augmentation MCMC (Demiris & O'Neill framework) reconstructing unobserved transmission chains. HI titer (pre-season antibody) as immunity marker: intermediate titer (1:20–40) reduces hazard by 59%; high titer (≥1:80) by 87%. Even after correcting for HI titers, adults have half the susceptibility of children. Household-size effect on β: inversely proportional to n_k (frequency-dependent). Average household transmission probability 8% (95% CI 6–10%).

**Relevance:**
- Closest geographic/epidemiological precedent to our VERDI study (South East Asia, household cohort, children included)
- The HI titer protective effect is the serological analogue of our `result_IgG_day1_mod` (anti-NCP IgG); aRR 0.42 in VERDI is consistent with their 59–87% reduction
- Their β(n_k) = β/n_k (frequency-dependent) is the γ = 1 end of our sensitivity analysis
- The Bayesian digraph approach for reconstructing transmission chains is a more general version of what we implement — useful reference when discussing our secondary-to-secondary chain modelling

### Li et al. (2026) — `hhbayes_2026.pdf`
**"HHBayes: A Flexible Bayesian Framework for Simulating and Analyzing Household Transmission Dynamics"**
*medRxiv* (preprint, not peer-reviewed)

**Core findings:** Open-source R package (`HHBayes`) with four core functionalities: (1) flexible simulation of household transmission data with customizable demographics, testing schedules, and viral kinetics; (2) Ct-value or viral-load-dependent time-varying infectiousness modelled via a biphasic viral load curve and logistic link; (3) Bayesian HMC parameter estimation via Stan (Hamiltonian Monte Carlo, 4 chains, 1000 warm-up + 1000 sampling); (4) visualization tools. Force of infection: λ_{ih}(t) = φ_i [α_comm S(t) + (1/max(n_h, 1))^δ Σ_{j≠i} κ_j c_{ij} βI(t,j)], where δ is household-size scaling. Covariates modify susceptibility φ_i and infectivity κ_j.

**Relevance:**
- The most direct methodological parallel to our model; key similarities and differences:
  - **Similar:** discrete-time hazard, community + household pathways, age contact matrix c_{ij}, density-dependent household term (δ ≡ our γ), Bayesian inference with data augmentation for infection intervals
  - **Different:** HMC (Stan) vs. our Metropolis-within-Gibbs; does not treat non-enrolled members as latent variables (our unique contribution); more complex viral kinetics model
- **Identifiability justification:** our `q · C[a_i, a_j]` is a restricted special case of their `β · φ_i · κ_j · c_{ij}` with φ ≡ 1, κ ≡ 1, δ = 0. With only 44 secondary cases in 93 households, separately estimating age-specific susceptibility (φ) and infectivity (κ) on top of the Prem matrix is practically non-identifiable. The single-q restriction is therefore a deliberate identifiability choice, not a simplification; this paper provides the general framework within which our model is a named, principled special case
- The δ parameter (household-size scaling) is the HHBayes formalisation of our planned γ sensitivity analysis (γ ∈ {0, 0.5, 1}); unlike ours, they treat δ as estimable — cite both as precedent for the sensitivity approach
- The Ct-value infectiousness function is the implementation reference for our planned β_Ct extension
- Their `prepare_stan_data()` imputation of infection intervals (uniform between last negative and first positive test) informs our interval-censoring treatment of `timetoinf`
- **Not peer-reviewed** — cite as preprint

### Tsang et al. (2023, PNAS) — `tsang_2023.pdf`
**"Reconstructing household transmission dynamics to estimate the infectiousness of asymptomatic influenza virus infections"**
*PNAS*, 120(33): e2304750120.

**Core findings:** 727 households, 2515 members, Hong Kong, 2009 pandemic H1N1. Bayesian model with additive hazard (community + household), allowing tertiary infections and asymptomatic cases. Children 3.2× more susceptible than adults. Symptomatic transmission probability 15% (children) and 5% (adults). Asymptomatic cases are 0.82× as infectious as symptomatic (relative infectiousness 0.57, 95% CrI 0.11–1.54; posterior probability of being less infectious = 82%). 12% of infections attributed to household transmission, 68% to symptomatic community contacts.

**Relevance:**
- Another application of the Cauchemez/Tsang household model framework to a large cohort — confirms our likelihood structure is appropriate at scale
- The asymptomatic infectiousness estimation is a methodological template for future extensions handling the heterogeneity between IgG-positive (partially immune, potentially asymptomatic) contacts in VERDI
- The community transmission probability (~12% attributed) contextualizes our P_BG calibration (12 community acquisitions out of 44 infected contacts = 27%, but our P_BG is a per-day rate not an attributable fraction)

### Ganyani et al. (2020) — `generation_interval_torneri2020.pdf`
**"Estimating the generation interval for coronavirus disease (COVID-19) based on symptom onset data, March 2020"**
*Eurosurveillance*, 25(17): pii=2000257.

**Core findings:** Bayesian MCMC (2-step: update missing transmission links v(i), then generation interval parameters Θ) to estimate generation interval from symptom onset data in Singapore (GI mean 5.20 days, SD 1.72) and Tianjin (GI mean 3.95 days, SD 1.51). Pre-symptomatic transmission fraction: 48% (Singapore) and 62% (Tianjin). Missing transmission links modelled as latent variables with equal prior probabilities; updated via independence sampler.

**Relevance:**
- The MCMC two-step pattern — (i) impute missing transmission links, (ii) update parameters — is the same architecture as our Gibbs updates for Y_hk (impute infection status) and q (transmission parameter)
- The generation interval estimates (~4–5 days) are consistent with Omicron serial intervals and inform the choice of infectious window D = 5 days in our model
- Pre-symptomatic transmission fractions motivate including the pre-enrollment infectious period τ ∈ {1, 2, 3} days in our model

### Walker et al. (2017) — see Section 2 above

---

## 7. Methodological Critique of Standard SAR Estimation

*(See also Section 5)*

### Cauchemez et al. (2009, NEJM) — see Section 2 above

Both Sharker & Kenah (2021) and Lindstrøm et al. (2024) demonstrate that the VERDI paper's statistical method (GEE log-binomial) is likely to produce an upward-biased SAR estimate. Our chain-binomial MCMC provides both a methodological correction and a principled imputation of missing members.

---

## Summary Table

| File | Authors | Year | Primary Use |
|---|---|---|---|
| `verdi_cnx_2026.pdf` | Khamduang et al. | 2026 | **Primary dataset** — source of all observations |
| `cauchemez2004.pdf` | Cauchemez et al. | 2004 | **MCMC backbone** — data augmentation of infection timing, 3-level hierarchy |
| `oneill1999.pdf` | O'Neill & Roberts | 1999 | **Foundational MCMC** — Bayesian inference for partially observed epidemics |
| `tsang_cauchemez_2023.pdf` | Tsang & Cauchemez | 2023 | **Primary likelihood** — additive hazard, escape probability, tertiary infections |
| `cauchemez2009.pdf` | Cauchemez et al. | 2009 | **Age stratification** — age-stratified susceptibility, tertiary chains |
| `cauchemez2014.pdf` | Cauchemez et al. | 2014 | **SEA context + IgG analogue** — HI titer protection, β/n, Vietnam |
| `prem2021.pdf` | Prem et al. | 2021 | **Contact matrix** — Thailand home matrix C[a_i, a_j] |
| `correc_prem2024.pdf` | Prem et al. | 2024 | **Contact matrix correction** — use corrected repository |
| `hhmem_nmar2018.pdf` | Goeyvaerts et al. | 2018 | **Density-dependence** — household-size contact density correction |
| `sharker2021.pdf` | Sharker & Kenah | 2021 | **Bias critique** — GEE/binomial SAR is upward biased; justifies our approach |
| `hh_incomplete_chain_bin_2024.pdf` | Lindstrøm et al. | 2024 | **Incomplete observation** — bias from short follow-up; ICB model |
| `hhbayes_2026.pdf` | Li et al. | 2026 | **Parallel approach** — HHBayes R package; β_Ct extension template |
| `tsang_2023.pdf` | Tsang et al. | 2023 (PNAS) | **Asymptomatic infectiousness** — model at scale; community attribution |
| `generation_interval_torneri2020.pdf` | Ganyani et al. | 2020 | **Generation interval + MCMC pattern** — two-step augmentation |
| `walker2017.pdf` | Walker et al. | 2017 | **Scalable MCMC** — DA-MCMC + branching process approximation |
| `zhang2026.pdf` | Zhang et al. | 2026 | **Contact matrix review** — validates fixed-C + scalar q; formalises non-identifiability of q and C jointly |
