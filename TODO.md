# TODO — VERDI missing-member correction

---

## 1. Check Patumrat's SAR computation

- [ ] Re-read Patumrat's analysis script / report and confirm exactly which members were included in the numerator and denominator.
- [ ] Check specifically: are missing (non-enrolled) members included in the denominator, or was the SAR computed on enrolled contacts only?
- [ ] Compare her denominator definition to the one used in the published paper (197 contacts, 93 index cases).
- [ ] If she used a different denominator, document the discrepancy and decide which definition the modeling project should align to.

---

## 2. Literature

### Read (done)

- [x] **Li et al. (2026)** — *"HHBayes: A Flexible Bayesian Framework for Simulating and Analyzing Household Transmission Dynamics"* — medRxiv preprint. `bib/hhbayes_2026.pdf`
  Key takeaways: (1) force of infection factorises as φ_i · κ_j · c_{ij} · βI(t,j) — separates susceptibility (φ_i) from infectivity (κ_j), the right framing for IgG (susceptibility) vs Ct (infectivity); (2) time-varying Ct trajectory via piecewise linear model Ct(t) plus logistic infectiousness function — our β_Ct uses only day-1 Ct; (3) HMC/Stan inference with R̂ and ESS diagnostics — no discrete latent variables, so cannot handle missing members; (4) gamma-distributed infectious duration per infector — same gap as Cauchemez 2004; (5) `simulate_multiple_households_comm()` is a ready-made forward simulator usable for our validation study; (6) no missing member imputation — our unique contribution.

- [x] **Cauchemez et al. (2004)** — *"A Bayesian MCMC approach to study transmission of influenza: application to household longitudinal data"* — Statistics in Medicine 23(22):3469–3487.
  Key takeaways: (1) continuous-time hazard λ_s(t) = α_s + ε_s Σ β_i/n — our model is the discrete-time analog; (2) augments start/end of infectious period (ν_i, ψ_i) within MCMC — we augment T_hk for missing members only; (3) no contact matrix, uses β_i/n normalization — we use Prem 2021 C_mat instead; (4) augments infectious duration — we fix D=5, which is a known gap.

- [x] **Ganyani et al. (2020)** — *"Estimating the generation interval for COVID-19"* — Eurosurveillance.
  Two-step MCMC: step 1 updates missing links, step 2 updates parameters. Confirms data augmentation pattern.

- [x] **Goeyvaerts et al. (2018)** — `bib/hhmem_nmar.pdf`
  Contact density decreases with household size; β/n normalization is approximately valid. Random mixing OK for simulation but biased for transmission rate inference.

### Still to read (high priority)

- [ ] **O'Neill & Roberts (1999)** — *"Bayesian inference for partially observed stochastic epidemics"* — JRSS-A 162(1):121–129.
  Theoretical basis for Gibbs-within-Metropolis with augmented latent states in epidemic models.

- [ ] **Cauchemez et al. (2009)** — *"Household transmission of 2009 pandemic influenza A (H1N1) virus in the United States"* — NEJM 361:2619–2627. Available at `bib/cauchemez2009.pdf`.
  Age-stratified household MCMC model; how they handle age-dependent transmission rates and multiple index cases.

### To skim (context / benchmarks)

- [ ] **Mossong et al. (2008)** — POLYMOD — *"Social contacts and mixing patterns relevant to the spread of infectious diseases"* — PLOS Medicine.
  Empirical foundation for contact matrices; useful for justifying Prem et al. prior.

- [ ] **Madewell et al. (2022)** — *"Household secondary attack rates of SARS-CoV-2 by variant and vaccination status"* — Open Forum Infectious Diseases.
  Meta-analysis; benchmark our corrected SAR against Omicron estimates elsewhere.

- [ ] **Bi et al. (2020)** — *"Epidemiology and transmission of COVID-19 in 391 cases and 1286 close contacts in Shenzhen"* — Lancet Infectious Diseases 20(8):911–919.
  COVID-19 SAR with partial observation; check methodology section.

### Targeted search

- [ ] PubMed / Google Scholar: `household transmission` + `missing data` + `data augmentation` + `MCMC` + `SARS-CoV-2`, years 2020–2025.
- [ ] Check whether any COVID-19 household study formally treated non-enrolled members as latent variables rather than excluding them.

---

## 3. Model implementation

### Done

- [x] Discrete-time timed chain-binomial likelihood (`compute_escape_seq`) — daily hazard from index, observed secondary infectors, missing infected members, community
- [x] Infection timing for observed contacts (`timetoinf`) — likelihood E_i(t_i−1) − E_i(t_i) for infected; E_i(T_follow) for uninfected
- [x] Secondary-to-secondary chains — HH 42 (timetoinf=7,10) now correctly attributed
- [x] β_Ct Ct infectiousness modifier — q_h = q × exp(β_Ct × (CT_MEAN − Ct_h)), MH update with N(0,1) prior
- [x] T_hk augmentation for missing members — exact Gibbs step enumerating days 1..14 (Y_hk=1 only)
- [x] MICS6 Chiang Mai age prior (`compute_mics_priors`) — marginal probability over Child/Adult/Elderly
- [x] age_k Gibbs update for missing members
- [x] Single MCMC chain with 5 update components (q, β_Ct, Y_hk, T_hk, age_k)
- [x] Pedagogical explainer document (`scripts/timed_chain_binomial_explainer.qmd`)

### Pending — model extensions

- [ ] **Run the model and check convergence** — verify MH acceptance rates (target 20–40% for q and β_Ct), trace plots, ESS for both q and β_Ct
- [ ] **Separate susceptibility and infectivity; add IgG modifier** — following Li et al. (2026) `bib/hhbayes_2026.pdf` φ_i · κ_j factorisation: IgG acts on susceptibility as φ_i^eff = exp(ψ_IgG · IgG_i) with ψ_IgG ≈ log(0.45) ≈ −0.80 (published aRR); Ct acts on infectivity as κ_j^eff = exp(β_Ct · (CT_MEAN − Ct_j)) — already implemented but currently conflated with q·C_mat. Missing members remain IgG-naive (φ_k = 1, conservative). This separates the two effects cleanly and is the highest-priority model extension.
- [ ] **Time-varying Ct infectiousness** — following Li et al. (2026) `bib/hhbayes_2026.pdf`: model Ct(t) as piecewise linear trajectory (rise to peak, then decay) and compute infectiousness f(Ct(t)) = [1 + exp((Ct(t) − vl_mid)/vl_slope)]^{−1} per day within the infectious window; we have only day-1 Ct so treat it as Ct_peak and use published Omicron kinetics for rise/decay slopes. Makes q_h(t) time-varying within D instead of constant. Lower priority than susceptibility/infectivity separation but strengthens the Ct component.
- [ ] **τ_h per-household augmentation** — treat test-to-enroll delay as discrete latent Uniform{0,1,2} per household; exact Gibbs step with 3 evaluations per iteration (future enhancement, not urgent)
- [ ] **Augment infectious duration D** — follow Cauchemez 2004 and Li et al. (2026) `bib/hhbayes_2026.pdf`: sample D_i from gamma prior per infector within MCMC instead of fixing D=5; moderate priority given fixed D absorbs heterogeneity into q

---

## 4. Sensitivity analyses (planned)

- [ ] **τ sensitivity** — run model at τ ∈ {1, 3} in addition to primary τ=2; report corrected SAR and β_Ct for each
- [ ] **P_BG sensitivity** — vary background hazard (0.001, 0.003, 0.006/day); quantify effect on corrected SAR and q
- [ ] **Multiple index case / co-primary sensitivity** — increase P_BG for households with early timetoinf contacts (potential co-primaries); assess impact on q and corrected SAR
- [ ] **Simulation validation study** — generate synthetic households with known ground truth, apply partial observation (remove fraction of members), check whether MCMC recovers true SAR; compare full-observation vs partial-observation estimates. Consider using Li et al. (2026) `bib/hhbayes_2026.pdf` `simulate_multiple_households_comm()` as the forward simulator rather than writing one from scratch.
- [ ] **Contact matrix sensitivity — density-dependent vs frequency-dependent** — the Prem C_mat gives population-average contact rates and implicitly assumes density-dependent transmission (hazard q_h · C[a_i, a_j] does not decrease with household size). Goeyvaerts et al. (2018) `bib/hhmem_nmar.pdf` showed empirically (318 Belgian households) that within-household contact density decreases with size (density ≈ 1.00 for n=2 down to ≈ 0.85 for n≥6) and that ignoring this biases transmission rate estimates. They do NOT provide a plug-in formula transferable to Thailand, but they frame the problem as: mean contact degree ∝ n^w with w ∈ (0,1). The two extreme models are w=0 (our current model, density-dependent) and w=1 (Cauchemez's β/n, frequency-dependent). Sensitivity analysis: re-run the MCMC with the effective hazard q_h · C[a_i, a_j] / n_h^γ for γ ∈ {0, 0.5, 1}, report corrected SAR and q under each. γ is not a free parameter to be estimated given n=93 households — fix it as a sensitivity knob. Note: the practical impact is likely modest for our households (most n=3–5, density difference ~6 percentage points per Goeyvaerts Table 1).

---

## 5. Email to the Belgium team

- [ ] Send email to Belgium team (cc study PIs).
- [ ] Follow up if no reply within 2 weeks.

### Draft email

---

**To:** [Belgium team]
**Cc:** [Study PIs]
**Subject:** VERDI-RECOVER Thailand — Update on the missing-member correction model

Dear [Name],

Thank you for your detailed and thoughtful reply. Your comments have been very helpful in shaping the methodological approach, and I wanted to share where things currently stand and get your further input on a few points.

**Where we are**

Following your suggestion, we have implemented a discrete-time chain-binomial model with Bayesian inference, which avoids the need for a full dynamical model while still rigorously accounting for the probability of escaping infection from the index case — as you correctly anticipated would be necessary. The model is fitted via Metropolis-within-Gibbs MCMC, where:

- The **transmission parameter** q (per-contact-per-day probability) is updated via random-walk Metropolis-Hastings;
- The **Ct infectiousness modifier** β_Ct scales q per household as q_h = q × exp(β_Ct × (CT_mean − Ct_h)), capturing the well-documented association between index viral load and transmission (aRR 0.84 per Ct unit in the published analysis);
- **Infection timing** is explicitly modelled via the `timetoinf` variable (days from enrollment to first positive PCR): the likelihood of contact i being infected on day t_i is E_i(t_i−1) − E_i(t_i), where E_i(t) is the cumulative escape probability computed day by day. This naturally accommodates secondary-to-secondary transmission chains;
- The **infection status Y_hk** and **infection day T_hk** of each missing member are treated as latent variables and updated via exact Gibbs steps at each iteration;
- The **age of missing members** is drawn from the MICS6 2022 Chiang Mai donor pool (hot-deck), updated as a Gibbs variable within the chain;
- The **corrected SAR** is computed from the joint posterior of (q, β_Ct, Y_hk, T_hk, age_k), expressed with a 95% credible interval.

This formulation is the discrete-time analog of the Cauchemez et al. (2004) continuous-time hazard model, with the contact matrix from Prem et al. (2021) substituting for their β_i/n parameterization.

**On age and contact rates**

You asked whether data on the age distribution by household size is available to infer the age of non-participating individuals. We addressed this using the UNICEF MICS6 2022 survey, filtered to Chiang Mai province, which provides the empirical joint distribution of household compositions (number of children, adults, and elderly) by total household size. We use this as a hot-deck donor pool: for each incomplete study household, we match to MICS6 donor households with the same total size and at least as many members in each age category as observed, then sample the full composition. This gives us a principled, Chiang Mai-specific prior for the age of missing members.

For contact rates by age, we use the Prem et al. 2021 synthetic home contact matrix for Thailand (16×16, 5-year age bands, aggregated to Child/Adult/Elderly). You mentioned that time-use data by age and sex could help reconstruct contact rates — we investigated this. The Prem matrix is age-stratified only, with no sex dimension. We did find a small Thai social contact survey (n=257, ages 14–52) that records contacts at home by sex of the respondent, but it is too small, covers too narrow an age range, and does not record the age or sex of contacts in a way that would allow construction of a sex-stratified matrix. We therefore retain Prem 2021 as the best available option and note the absence of sex-stratified home contact rates as a limitation.

**On the index pre-enrollment infectious period**

One nuance specific to our study design: index cases were enrolled within 48 hours of their positive test. This constrains the pre-enrollment infectious period τ to at most 2 days from test to enrollment, plus ~0–1 days of pre-symptomatic infectiousness before testing — so τ ∈ {1, 2, 3} days in total. We use τ = 2 as the primary estimate, with sensitivity analysis at τ ∈ {1, 3}.

**On community acquisition**

You suggested that population-level prevalence data could help estimate the probability of becoming infected outside the household. We investigated provincial surveillance data from Thailand's DDC, but found that Thailand moved to endemic management in July 2022 — exactly when enrolment started — after which official case counts became severe undercounts with no stable correction factor. We therefore retain our current estimate of P_BG = 0.003/day, derived from the study's own phylogenetic analysis (12 confirmed community acquisitions in 40 infected contacts over approximately 21 days), which we consider more reliable and specific to this population than external surveillance data.

**On multiple index cases**

Your point about multiple index cases as an additional source of bias is well-taken and not yet fully addressed in the model. The current framework treats the enrolled index case as the sole primary infector. However, some contacts who tested positive during follow-up may have been co-primary cases independently exposed from the community — the phylogenetic analysis confirmed 12 such cases. The background hazard term P_BG partially accounts for this, but if a contact was truly a co-primary case, attributing their infection to household transmission inflates q and consequently overestimates P(Y_hk = 1) for missing members. We plan to address this in a sensitivity analysis.

**On a simulation validation study**

Your suggestion of comparing estimated SARs when the contact network is fully vs. partially observed resonates with us as a valuable methodological contribution. We are planning a simulation study using the fitted model in forward mode: generate synthetic households with known ground truth, apply partial observation (removing a fraction of members at random), and check whether the MCMC recovers the true SAR.

**Key remaining questions for you**

1. The largest unresolved bias in our model is the absence of IgG data for missing members. IgG positivity was the strongest protective factor in the published analysis (aRR 0.45). Missing members are assumed IgG-naive, which conservatively overestimates their probability of being infected. We are uncertain whether to (a) treat this as a fixed limitation and report a conservative upper bound on the corrected SAR, or (b) impute a population IgG prevalence for missing members based on the enrolled contacts and incorporate it as an additional sensitivity analysis. Do you have a preference, or have you encountered a similar problem in your studies?

2. Do you see any concern with the MCAR assumption for the missingness mechanism? The main barrier to participation was fear of venipuncture (blood draw for IgG). If IgG-positive individuals were systematically less likely to fear the blood draw — for example, because prior infection made them more comfortable with medical procedures — then MCAR would be violated in a way that further overestimates infections among missing members.

3. Are you aware of any COVID-19 household transmission study that formally treated non-enrolled members as latent variables in a transmission model, rather than simply excluding them from the denominator? A PubMed search has not yet turned up a close precedent.

We would very much welcome your thoughts, and are happy to share the current R/Quarto implementation if that would be useful for your review.

Best regards,
Mathis

---
