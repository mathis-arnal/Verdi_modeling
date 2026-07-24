# Simplified chain-binomial model — no age stratification, no contact matrix.
# Sourced by scripts/chain_binomial_simple.qmd.
#
# Hazard: λ_i(t) = λ_c + q·f(τ+t) + Σ_j q·f(t−t_j+1)
# Latent variables per missing member: (Y_hk, T_hk) only — no age_k.
# Contrast with model_functions.R which adds C_mat and age_k.

library(dplyr)

# ── Data preparation ──────────────────────────────────────────────────────────

build_hh_list_simple <- function(data, tau = 2) {
  lapply(unique(data$hhid), function(h) {
    hh  <- data %>% filter(hhid == h)
    list(
      hhid      = h,
      tau       = tau,
      contacts  = hh %>%
                    filter(index == 0, !is.na(infected_bin)) %>%
                    select(infected_bin, timetoinf),
      n_missing = first(hh$exclude_mem)
    )
  })
}

# ── Likelihood ────────────────────────────────────────────────────────────────

# Returns escape probability sequence of length T_follow + 1.
# esc[1] = 1 (before any exposure); esc[t+1] = P(susceptible through day t).
# infectors: list of list(start = first_infectious_day).
# f_prof: discretized Gamma infectiousness profile (normalized, length D_max).
compute_escape_seq_simple <- function(q_h, tau, infectors, T_follow,
                                       p_bg = 0.003, f_prof) {
  esc   <- numeric(T_follow + 1L)
  esc[1L] <- 1.0
  D_max <- length(f_prof)
  lam_c <- -log(1 - p_bg)

  for (t in seq_len(T_follow)) {
    lam <- lam_c
    s_idx <- tau + t
    if (s_idx >= 1L && s_idx <= D_max)
      lam <- lam + q_h * f_prof[s_idx]
    for (inf in infectors) {
      s_inf <- t - inf$start + 1L
      if (s_inf >= 1L && s_inf <= D_max)
        lam <- lam + q_h * f_prof[s_inf]
    }
    esc[t + 1L] <- esc[t] * exp(-lam)
  }
  esc
}

# miss_Y: integer vector 0/1, length h$n_missing
# miss_T: integer vector 1..T_follow (NA when miss_Y = 0)
hh_loglik_timed_simple <- function(q, h, miss_Y, miss_T,
                                    f_prof, p_bg = 0.003, T_follow = 14) {
  n_c <- nrow(h$contacts)
  if (n_c == 0L) return(0)
  q_h <- max(q, 1e-9)

  obs_inf <- lapply(seq_len(n_c), function(i) {
    if (h$contacts$infected_bin[i] == 1L && !is.na(h$contacts$timetoinf[i]))
      list(start = as.integer(h$contacts$timetoinf[i]))
    else NULL
  })

  miss_inf <- list()
  for (k in seq_along(miss_Y)) {
    if (!is.na(miss_Y[k]) && miss_Y[k] == 1L && !is.na(miss_T[k]))
      miss_inf[[length(miss_inf) + 1L]] <- list(start = as.integer(miss_T[k]))
  }

  loglik <- 0

  for (i in seq_len(n_c)) {
    infected    <- h$contacts$infected_bin[i]
    t_inf       <- h$contacts$timetoinf[i]
    infectors_i <- c(Filter(Negate(is.null), obs_inf[-i]), miss_inf)
    esc         <- compute_escape_seq_simple(q_h, h$tau, infectors_i,
                                             T_follow, p_bg, f_prof)

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

total_loglik_timed_simple <- function(q, hh_list, Y_state, T_state,
                                       f_prof, p_bg = 0.003, T_follow = 14) {
  ll <- 0
  for (h in hh_list) {
    hid <- as.character(h$hhid)
    if (h$n_missing > 0)
      ll <- ll + hh_loglik_timed_simple(q, h,
                                         Y_state[[hid]], T_state[[hid]],
                                         f_prof, p_bg, T_follow)
    else
      ll <- ll + hh_loglik_timed_simple(q, h,
                                         integer(0L), integer(0L),
                                         f_prof, p_bg, T_follow)
  }
  ll
}

# ── Missing member prior ──────────────────────────────────────────────────────

# Prior for a missing member: index + community + enrolled secondary contacts.
# Enrolled infectors are observed and fixed — including them adds no circular
# dependence and correctly captures secondary-chain exposure for missing members.
prior_miss_simple <- function(q_h, tau, enrolled_infectors = list(),
                               f_prof, p_bg = 0.003, T_follow = 14) {
  esc     <- compute_escape_seq_simple(q_h, tau, enrolled_infectors, T_follow, p_bg, f_prof)
  T_probs <- pmax(esc[1L:T_follow] - esc[2L:(T_follow + 1L)], 0)
  list(Y0 = esc[T_follow + 1L], T_probs = T_probs)
}

# ── MCMC sampler ─────────────────────────────────────────────────────────────

# Two update components per iteration (no age_k step):
#   1. MH for q (log-scale random walk, Uniform(0,1) prior — i.e. Beta(1,1),
#      matching the vague-prior convention of Cauchemez et al. and
#      Tsang & Cauchemez (2023) rather than an informative Beta(1,19))
#   2. Gibbs for (Y_hk, T_hk) jointly (exact (T_follow+1)-way enumeration:
#      "not infected" vs "infected on day t" for t = 1..T_follow). Updating
#      Y_hk against a single candidate T (rather than marginalizing over all
#      T) understates the true P(Y_hk=1), which is a sum over T, not a point
#      evaluation — so Y_hk and T_hk must be drawn jointly, not sequentially.
run_mcmc_simple <- function(hh_list,
                             n_iter      = 10000,
                             n_burnin    = 2000,
                             n_thin      = 5,
                             q_init      = 0.02,
                             q_sd        = 0.015,
                             alpha_q     = 1,
                             beta_q      = 1,
                             y_init_prob = 0.2,
                             T_follow    = 14,
                             p_bg        = 0.003,
                             f_prof) {

  # ── Initialise latent states ─────────────────────────────────────────────
  # Y_hk^(0) ~ Bernoulli(y_init_prob), fixed at ~observed SAR (0.2);
  # T_hk^(0) ~ empirical timetoinf distribution among infected enrolled
  # contacts. This is a starting value only — the model-implied
  # prior_miss_simple distribution is still what governs the Y_hk/T_hk
  # Gibbs updates every iteration; burn-in washes out the initialization
  # choice, so a simple, fixed observed-data-based start is used.
  timetoinf_obs <- unlist(lapply(hh_list, function(h) {
    h$contacts$timetoinf[h$contacts$infected_bin == 1L & !is.na(h$contacts$timetoinf)]
  }))

  Y_state <- list()
  T_state <- list()

  for (h in hh_list) {
    if (h$n_missing > 0) {
      hid <- as.character(h$hhid)
      Y0  <- rbinom(h$n_missing, 1L, y_init_prob)
      T0  <- ifelse(Y0 == 1L,
                    sample(timetoinf_obs, h$n_missing, replace = TRUE),
                    NA_integer_)
      Y_state[[hid]] <- Y0
      T_state[[hid]] <- T0
    }
  }

  q       <- q_init
  ll      <- total_loglik_timed_simple(q, hh_list, Y_state, T_state,
                                        f_prof = f_prof, p_bg = p_bg,
                                        T_follow = T_follow)
  n_acc_q <- 0L

  n_keep          <- floor((n_iter - n_burnin) / n_thin)
  q_chain         <- numeric(n_keep)
  E_miss_infected <- numeric(n_keep)
  keep_idx        <- 0L

  for (iter in seq_len(n_iter)) {

    # ── 1. MH update for q ────────────────────────────────────────────────
    # Proposal is a symmetric random walk on log(q); the change of variables
    # from q to log(q) contributes a Jacobian factor of q, so the prior-ratio
    # exponent on log(q_prop/q) is alpha_q (not alpha_q - 1, which would be
    # the correct exponent only for an untransformed/linear-scale proposal).
    q_prop <- exp(log(q) + rnorm(1L, 0, q_sd))
    if (q_prop > 0 && q_prop < 1) {
      ll_prop <- total_loglik_timed_simple(q_prop, hh_list, Y_state, T_state,
                                            f_prof = f_prof, p_bg = p_bg,
                                            T_follow = T_follow)
      log_pr  <- alpha_q * (log(q_prop) - log(q)) +
                 (beta_q  - 1) * (log(1 - q_prop) - log(1 - q))
      if (log(runif(1L)) < (ll_prop - ll) + log_pr) {
        q <- q_prop; ll <- ll_prop; n_acc_q <- n_acc_q + 1L
      }
    }

    # ── 2. Gibbs for (Y_hk, T_hk) jointly ────────────────────────────────
    # Exact (T_follow + 1)-way categorical: state 1 = "not infected",
    # state t+1 = "infected on day t". Enumerating and sampling all states
    # in one draw is the exact full conditional for the pair (Y_hk, T_hk);
    # updating Y_hk against only one candidate T would instead compare
    # against a single point of the Y=1 branch rather than its full mass.
    for (h in hh_list) {
      if (h$n_missing == 0L) next
      hid <- as.character(h$hhid)
      q_h <- max(q, 1e-9)

      enrolled_inf <- Filter(Negate(is.null),
        lapply(seq_len(nrow(h$contacts)), function(i) {
          if (h$contacts$infected_bin[i] == 1L && !is.na(h$contacts$timetoinf[i]))
            list(start = as.integer(h$contacts$timetoinf[i]))
          else NULL
        }))

      for (k in seq_len(h$n_missing)) {
        pm <- prior_miss_simple(q_h, h$tau, enrolled_inf, f_prof, p_bg, T_follow)

        log_probs <- numeric(T_follow + 1L)  # [1] = Y=0; [t+1] = Y=1, T=t

        mY0 <- Y_state[[hid]]; mY0[k] <- 0L
        mT0 <- T_state[[hid]]; mT0[k] <- NA_integer_
        ll0 <- hh_loglik_timed_simple(q, h, mY0, mT0, f_prof, p_bg, T_follow)
        log_probs[1L] <- ll0 + log(max(pm$Y0, 1e-12))

        for (t_try in seq_len(T_follow)) {
          mY1 <- Y_state[[hid]]; mY1[k] <- 1L
          mT1 <- T_state[[hid]]; mT1[k] <- t_try
          ll1 <- hh_loglik_timed_simple(q, h, mY1, mT1, f_prof, p_bg, T_follow)
          log_probs[t_try + 1L] <- ll1 + log(max(pm$T_probs[t_try], 1e-12))
        }

        lp_mx <- max(log_probs)
        probs <- exp(log_probs - lp_mx)
        probs <- probs / sum(probs)
        state <- sample.int(T_follow + 1L, 1L, prob = probs)

        if (state == 1L) {
          Y_state[[hid]][k] <- 0L
          T_state[[hid]][k] <- NA_integer_
        } else {
          Y_state[[hid]][k] <- 1L
          T_state[[hid]][k] <- state - 1L
        }
      }
    }

    ll <- total_loglik_timed_simple(q, hh_list, Y_state, T_state,
                                     f_prof = f_prof, p_bg = p_bg,
                                     T_follow = T_follow)

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
