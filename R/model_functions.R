# Model functions for the VERDI household SARS-CoV-2 transmission model.
# Sourced by scripts/chain_binomial_mcmc.qmd and scripts/sensitivity_analysis.qmd.
#
# Default argument values match the primary analysis (Omicron biology, 14-day follow-up,
# p_bg from phylogenetic data). Override them explicitly in sensitivity runs.

library(dplyr)

AGE_CATS <- c("Child", "Adult", "Elderly")

# ── Data preparation ──────────────────────────────────────────────────────────

build_hh_list <- function(data, tau = 2) {
  lapply(unique(data$hhid), function(h) {
    hh  <- data %>% filter(hhid == h)
    idx <- hh %>% filter(index == 1)
    list(
      hhid      = h,
      index_age = as.character(idx$age_cat[1]),
      tau       = tau,
      contacts  = hh %>%
                    filter(index == 0, !is.na(infected_bin)) %>%
                    select(age_cat, infected_bin, timetoinf),
      n_missing = first(hh$exclude_mem)
    )
  })
}

# ── MICS hot-deck age prior ───────────────────────────────────────────────────

find_donors <- function(obs_c, obs_a, obs_e, hh_size, pool) {
  pool %>%
    filter(household_size == hh_size,
           n_child >= obs_c, n_adult >= obs_a, n_elderly >= obs_e) %>%
    mutate(freq = count / sum(count))
}

# Returns a named list keyed by hhid (character).
# Each element is a probability vector over (Child, Adult, Elderly) for missing members.
compute_mics_priors <- function(hh_list, comp_table, data) {
  priors <- list()
  for (h in hh_list) {
    if (h$n_missing == 0) next
    hid   <- as.character(h$hhid)
    obs   <- data %>% filter(hhid == h$hhid, index == 0)
    obs_c <- sum(obs$age_cat == "Child",   na.rm = TRUE)
    obs_a <- sum(obs$age_cat == "Adult",   na.rm = TRUE)
    obs_e <- sum(obs$age_cat == "Elderly", na.rm = TRUE)
    hh_sz <- nrow(obs) + h$n_missing + 1L
    donors <- find_donors(obs_c, obs_a, obs_e, hh_sz, comp_table)
    if (nrow(donors) == 0) {
      priors[[hid]] <- c(Child = 0.15, Adult = 0.75, Elderly = 0.10)
    } else {
      extra_c <- sum(donors$freq * pmax(0, donors$n_child   - obs_c))
      extra_a <- sum(donors$freq * pmax(0, donors$n_adult   - obs_a))
      extra_e <- sum(donors$freq * pmax(0, donors$n_elderly - obs_e))
      total   <- extra_c + extra_a + extra_e
      if (total < 1e-9)
        priors[[hid]] <- c(Child = 0.15, Adult = 0.75, Elderly = 0.10)
      else
        priors[[hid]] <- c(Child   = extra_c / total,
                           Adult   = extra_a / total,
                           Elderly = extra_e / total)
    }
  }
  priors
}

# ── Likelihood ────────────────────────────────────────────────────────────────

# Returns escape probability sequence of length T_follow + 1.
# esc[t+1] = P(susceptible through day t);  esc[1] = 1 (day 0, before any exposure).
# infectors: list of list(age = "Child"/"Adult"/"Elderly", start = first_infectious_day).
# f_prof: discretized Gamma infectiousness profile (normalized, length D_max).
#   f_prof[s] = relative infectiousness on day s of the infectious period.
#   Index: at follow-up day t, index is at infectious day tau + t  →  f_prof[tau + t].
#   Secondary infector infected at day t_j: at day t they are at day t - t_j + 1  →  f_prof[t - t_j + 1].
compute_escape_seq <- function(q_h, sus_age, index_age, tau,
                                infectors, T_follow, C_mat,
                                p_bg = 0.003, f_prof) {
  esc   <- numeric(T_follow + 1L)
  esc[1L] <- 1.0
  D_max <- length(f_prof)
  lam_c <- -log(1 - p_bg)

  for (t in seq_len(T_follow)) {
    lam <- lam_c

    s_idx <- tau + t
    if (s_idx >= 1L && s_idx <= D_max)
      lam <- lam + q_h * C_mat[sus_age, index_age] * f_prof[s_idx]

    for (inf in infectors) {
      s_inf <- t - inf$start + 1L
      if (s_inf >= 1L && s_inf <= D_max)
        lam <- lam + q_h * C_mat[sus_age, inf$age] * f_prof[s_inf]
    }

    esc[t + 1L] <- esc[t] * exp(-lam)
  }
  esc
}

# miss_Y  : integer vector 0/1, length h$n_missing
# miss_T  : integer vector 1..T_follow (NA when miss_Y = 0)
# miss_age: character vector
hh_loglik_timed <- function(q, h, miss_Y, miss_T, miss_age, C_mat,
                             f_prof, p_bg = 0.003, T_follow = 14) {
  n_c <- nrow(h$contacts)
  if (n_c == 0L) return(0)

  q_h <- max(q, 1e-9)

  obs_inf <- lapply(seq_len(n_c), function(i) {
    if (h$contacts$infected_bin[i] == 1L && !is.na(h$contacts$timetoinf[i]))
      list(age   = as.character(h$contacts$age_cat[i]),
           start = as.integer(h$contacts$timetoinf[i]))
    else NULL
  })

  miss_inf <- list()
  for (k in seq_along(miss_Y)) {
    if (!is.na(miss_Y[k]) && miss_Y[k] == 1L && !is.na(miss_T[k]))
      miss_inf[[length(miss_inf) + 1L]] <- list(age   = miss_age[k],
                                                  start = as.integer(miss_T[k]))
  }

  loglik <- 0

  for (i in seq_len(n_c)) {
    sus_age  <- as.character(h$contacts$age_cat[i])
    infected <- h$contacts$infected_bin[i]
    t_inf    <- h$contacts$timetoinf[i]

    infectors_i <- c(Filter(Negate(is.null), obs_inf[-i]), miss_inf)

    esc <- compute_escape_seq(q_h, sus_age, h$index_age, h$tau,
                               infectors_i, T_follow, C_mat, p_bg, f_prof)

    loglik <- loglik + if (infected == 0L) {
      log(max(esc[T_follow + 1L], 1e-12))
    } else if (!is.na(t_inf)) {
      log(max(esc[t_inf] - esc[t_inf + 1L], 1e-12))
    } else {
      log(max(1 - esc[T_follow + 1L], 1e-12))
    }
  }

  loglik
}

total_loglik_timed <- function(q, hh_list, Y_state, T_state, age_state, C_mat,
                                f_prof, p_bg = 0.003, T_follow = 14) {
  ll <- 0
  for (h in hh_list) {
    hid <- as.character(h$hhid)
    if (h$n_missing > 0)
      ll <- ll + hh_loglik_timed(q, h,
                                  Y_state[[hid]], T_state[[hid]], age_state[[hid]],
                                  C_mat, f_prof, p_bg, T_follow)
    else
      ll <- ll + hh_loglik_timed(q, h,
                                  integer(0L), integer(0L), character(0L),
                                  C_mat, f_prof, p_bg, T_follow)
  }
  ll
}

# ── Missing member prior ──────────────────────────────────────────────────────

# Returns list(Y0 = P(Y=0),
#              T_probs = length-T_follow vector where T_probs[t] = P(Y=1, T=t)).
# Uses index + community hazard only — no secondary infectors — to avoid
# circular dependence in the Gibbs update.
prior_miss <- function(q_h, miss_age, index_age, tau, C_mat,
                       f_prof, p_bg = 0.003, T_follow = 14) {
  esc <- compute_escape_seq(q_h, miss_age, index_age, tau, list(),
                             T_follow, C_mat, p_bg, f_prof)
  T_probs <- pmax(esc[1L:T_follow] - esc[2L:(T_follow + 1L)], 0)
  list(Y0 = esc[T_follow + 1L], T_probs = T_probs)
}

# ── MCMC sampler ─────────────────────────────────────────────────────────────

run_mcmc <- function(hh_list, mics_priors, C_mat,
                     n_iter   = 10000,
                     n_burnin = 2000,
                     n_thin   = 5,
                     q_init   = 0.02,
                     q_sd     = 0.015,
                     alpha_q  = 1,
                     beta_q   = 19,
                     T_follow = 14,
                     p_bg     = 0.003,
                     f_prof) {

  # ── Initialise latent states ─────────────────────────────────────────────
  Y_state   <- list()
  T_state   <- list()
  age_state <- list()

  for (h in hh_list) {
    if (h$n_missing > 0) {
      hid  <- as.character(h$hhid)
      Y0   <- rbinom(h$n_missing, 1L, 0.3)
      age0 <- sample(AGE_CATS, h$n_missing, replace = TRUE, prob = mics_priors[[hid]])
      T0   <- ifelse(Y0 == 1L,
                     sample.int(T_follow, h$n_missing, replace = TRUE),
                     NA_integer_)
      Y_state[[hid]]   <- Y0
      T_state[[hid]]   <- T0
      age_state[[hid]] <- age0
    }
  }

  q       <- q_init
  ll      <- total_loglik_timed(q, hh_list, Y_state, T_state, age_state, C_mat,
                                 f_prof = f_prof, p_bg = p_bg, T_follow = T_follow)
  n_acc_q <- 0L

  n_keep          <- floor((n_iter - n_burnin) / n_thin)
  q_chain         <- numeric(n_keep)
  E_miss_infected <- numeric(n_keep)
  keep_idx        <- 0L

  for (iter in seq_len(n_iter)) {

    # ── 1. MH update for q ────────────────────────────────────────────────
    q_prop <- exp(log(q) + rnorm(1L, 0, q_sd))
    if (q_prop > 0 && q_prop < 1) {
      ll_prop <- total_loglik_timed(q_prop, hh_list,
                                    Y_state, T_state, age_state, C_mat,
                                    f_prof = f_prof, p_bg = p_bg, T_follow = T_follow)
      log_pr  <- (alpha_q - 1) * (log(q_prop) - log(q)) +
                 (beta_q  - 1) * (log(1 - q_prop) - log(1 - q))
      if (log(runif(1L)) < (ll_prop - ll) + log_pr) {
        q <- q_prop; ll <- ll_prop; n_acc_q <- n_acc_q + 1L
      }
    }

    # ── 2. Gibbs for Y_hk ─────────────────────────────────────────────────
    for (h in hh_list) {
      if (h$n_missing == 0L) next
      hid <- as.character(h$hhid)
      q_h <- max(q, 1e-9)

      for (k in seq_len(h$n_missing)) {
        age_k <- age_state[[hid]][k]
        pm    <- prior_miss(q_h, age_k, h$index_age, h$tau, C_mat,
                            f_prof = f_prof, p_bg = p_bg, T_follow = T_follow)

        T_cur  <- T_state[[hid]][k]
        T_prop <- if (!is.na(T_cur)) T_cur else {
          tot <- sum(pm$T_probs)
          if (tot > 1e-12) sample.int(T_follow, 1L, prob = pm$T_probs / tot)
          else             sample.int(T_follow, 1L)
        }

        mY0 <- Y_state[[hid]]; mY0[k] <- 0L
        mT0 <- T_state[[hid]]; mT0[k] <- NA_integer_
        mY1 <- Y_state[[hid]]; mY1[k] <- 1L
        mT1 <- T_state[[hid]]; mT1[k] <- T_prop

        ll0 <- hh_loglik_timed(q, h, mY0, mT0, age_state[[hid]], C_mat,
                                f_prof = f_prof, p_bg = p_bg, T_follow = T_follow)
        ll1 <- hh_loglik_timed(q, h, mY1, mT1, age_state[[hid]], C_mat,
                                f_prof = f_prof, p_bg = p_bg, T_follow = T_follow)

        lp0 <- ll0 + log(max(pm$Y0, 1e-12))
        lp1 <- ll1 + log(max(pm$T_probs[T_prop], 1e-12))

        mx    <- max(lp0, lp1)
        p1    <- exp(lp1 - mx) / (exp(lp0 - mx) + exp(lp1 - mx))
        new_Y <- rbinom(1L, 1L, p1)

        Y_state[[hid]][k] <- new_Y
        T_state[[hid]][k] <- if (new_Y == 1L) T_prop else NA_integer_
      }
    }

    # ── 3. Gibbs for T_hk ─────────────────────────────────────────────────
    for (h in hh_list) {
      if (h$n_missing == 0L) next
      hid <- as.character(h$hhid)
      q_h <- max(q, 1e-9)

      for (k in seq_len(h$n_missing)) {
        if (Y_state[[hid]][k] != 1L) next
        age_k <- age_state[[hid]][k]
        pm    <- prior_miss(q_h, age_k, h$index_age, h$tau, C_mat,
                            f_prof = f_prof, p_bg = p_bg, T_follow = T_follow)

        log_probs_t <- numeric(T_follow)
        for (t_try in seq_len(T_follow)) {
          mT_tmp             <- T_state[[hid]]; mT_tmp[k] <- t_try
          ll_t               <- hh_loglik_timed(q, h,
                                                 Y_state[[hid]], mT_tmp,
                                                 age_state[[hid]], C_mat,
                                                 f_prof = f_prof, p_bg = p_bg, T_follow = T_follow)
          log_probs_t[t_try] <- ll_t + log(max(pm$T_probs[t_try], 1e-12))
        }

        lp_mx   <- max(log_probs_t)
        probs_t <- exp(log_probs_t - lp_mx)
        probs_t <- probs_t / sum(probs_t)
        T_state[[hid]][k] <- sample.int(T_follow, 1L, prob = probs_t)
      }
    }

    # ── 4. Gibbs for age_k ────────────────────────────────────────────────
    for (h in hh_list) {
      if (h$n_missing == 0L) next
      hid <- as.character(h$hhid)
      q_h <- max(q, 1e-9)

      for (k in seq_len(h$n_missing)) {
        y_k <- Y_state[[hid]][k]
        t_k <- T_state[[hid]][k]

        log_probs <- sapply(AGE_CATS, function(a) {
          ages_tmp    <- age_state[[hid]]; ages_tmp[k] <- a
          pm          <- prior_miss(q_h, a, h$index_age, h$tau, C_mat,
                                    f_prof = f_prof, p_bg = p_bg, T_follow = T_follow)
          log_p_prior <- if (y_k == 0L)
                           log(max(pm$Y0, 1e-12))
                         else
                           log(max(pm$T_probs[t_k], 1e-12))
          hh_loglik_timed(q, h,
                          Y_state[[hid]], T_state[[hid]], ages_tmp, C_mat,
                          f_prof = f_prof, p_bg = p_bg, T_follow = T_follow) +
            log_p_prior +
            log(max(mics_priors[[hid]][a], 1e-12))
        })

        lp_mx <- max(log_probs)
        probs  <- exp(log_probs - lp_mx)
        probs  <- probs / sum(probs)
        age_state[[hid]][k] <- sample(AGE_CATS, 1L, prob = probs)
      }
    }

    ll <- total_loglik_timed(q, hh_list, Y_state, T_state, age_state, C_mat,
                              f_prof = f_prof, p_bg = p_bg, T_follow = T_follow)

    # ── Store post-burnin thinned samples ─────────────────────────────────
    if (iter > n_burnin && (iter - n_burnin) %% n_thin == 0L) {
      keep_idx <- keep_idx + 1L
      q_chain[keep_idx] <- q
      E_miss_infected[keep_idx] <- sum(sapply(hh_list, function(hh) {
        if (hh$n_missing == 0L) return(0L)
        sum(Y_state[[as.character(hh$hhid)]], na.rm = TRUE)
      }))
    }
  }

  cat(sprintf("MH acceptance rate — q: %.1f%%\n", 100 * n_acc_q / n_iter))

  list(
    q_chain         = q_chain[seq_len(keep_idx)],
    E_miss_infected = E_miss_infected[seq_len(keep_idx)]
  )
}
