# LOAD IN PACKAGES

library(data.table)
library(depmixS4)
library(stats)
library(parallel)

options(stringsAsFactors = FALSE) # keeps text as text 

# CREATE PROGRESS LOG MESSAGES

TS <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S") 

log_msg <- function(..., sep = "") {
  message(sprintf("[%s] %s", TS(), paste0(..., collapse = sep))) 
}

time_block <- function(label, expr) {
  log_msg(label, " ...")
  tm <- system.time({
    val <- eval.parent(substitute(expr))
  })
  log_msg(label, sprintf(" done (elapsed %.2f sec)", tm[["elapsed"]]))
  val
}

# DEFINE GLOBAL SETTINGS FOR MODELS

MAX_Z_BETA  <- 5  # exp(5) ~ 148; prevents beta exploding during optimisation
MAX_Z_ALPHA <- 10  # plogis(10) ~ 0.99995; prevents alpha too close to extreme edges
MAX_KAPPA   <- 8
MIN_SEG_LEN <- 5  # minimum number of trials to attempt fits

# OPTIONAL IF NEEDED FOR STABILISATION/DEBUGGING

USE_BETA_PENALTY <- FALSE   # keep FALSE for fair AIC/BIC comparisons across models

# SETTINGS FOR EMA-BASED RECENT REWARD RATE

EMA_LAMBDA <- 0.2 # how quickly estimate updates 
EMA_INIT   <- 0.5 # starting value
RR_CENTER  <- 0.5 # reference midpoint

ema_reward_rate <- function(reward, lambda = EMA_LAMBDA, init = EMA_INIT) {
  n <- length(reward)
  rr <- rep(init, n)
  if (n <= 1) return(rr)
  for (t in 2:n) {
    rr[t] <- (1 - lambda) * rr[t-1] + lambda * reward[t-1]
  }
  rr
}

# CREATE LOGGING TABLE FOR DEBUGGING AND FUNCTION FOR LOGGING FAILURES 

debug_failures <- data.table(
  seg_id     = character(),
  model      = character(),
  message    = character(),
  start_par  = list(),
  result_par = list()
)

log_failure <- function(segid, model, message, start_par = NULL, result_par = NULL) {
  newrow <- data.table(
    seg_id  = as.character(segid),
    model   = as.character(model),
    message = as.character(message)
  )
  newrow[, start_par := list(list(start_par))]
  newrow[, result_par := list(list(result_par))]
  debug_failures <<- rbindlist(list(debug_failures, newrow), use.names = TRUE, fill = TRUE)
}

# LEARNING ENGINES FOR REINFORCEMENT LEARNING MODELS WITH NO RESET AT REVERSAL 

# q-learning model with separate positive and negative learning rates 

dual_alpha_Q <- function(choice, reward, alpha_pos, alpha_neg, init_Q = c(0,0)){
  n <- length(choice) #number of trials
  Q <- matrix(0, n, 2) # table of trials (row = n) and choices (columns = 2 = left vs right)
  Q[1,] <- init_Q # start both choices at the defined initial values
  rpe <- rep(NA_real_, n) # create reward prediction error vector
  if (n < 1) return(list(Q=Q, rpe=rpe)) # safety if no trials 
  for (t in 1:n){ # for all trials
    a <- as.integer(choice[t]) + 1 # converts 0/1 in data to 1/2 for R 
    r <- reward[t] # reward
    rpe[t] <- r - Q[t, a] # compute rpe or difference between expected reward and actual reward
    lr <- ifelse(is.na(rpe[t]), NA, ifelse(rpe[t] >= 0, alpha_pos, alpha_neg)) # if rpe pos, then alpha_pos or if rpe neg, then alpha_neg
    if (t < n) {
      Q[t+1, ] <- Q[t, ]
      if (!is.null(lr) && !is.na(lr)) Q[t+1, a] <- Q[t, a] + lr * rpe[t]
    } # update the choice's Q-value using learning rate x prediction error 
  }
  return(list(Q=Q, rpe=rpe)) # assign full history of learned values and prediction errors
}

# q-learning model with single learning rate 

single_alpha_Q <- function(choice, reward, alpha, init_Q = c(0,0)){
  dual_alpha_Q(choice, reward, alpha, alpha, init_Q)
} # as before but uses the same alpha rate for both positive and negative outcomes

# SOFTMAX TO TURN INTERNAL VALUES INTO CHOICE PROBABILITIES

softmax_probs <- function(Qt_row, beta, prev_choice = NULL, kappa = 0, epsilon = 0){
  logits <- beta * Qt_row # scales q-values by betas before softmax 
  if (!is.null(prev_choice) && !is.na(prev_choice)) {
    logits[as.integer(prev_choice) + 1] <- logits[as.integer(prev_choice) + 1] + kappa
  } # if there was a prior choice, add kappa (+ kappa = repeat choice, - kappa = change choice)
  mx <- max(logits, na.rm = TRUE) 
  ex <- exp(logits - mx) # subtract max exponential keeps numbers manageable for computing but maintains relative differences
  p1 <- ex[2] / sum(ex) # probability of the second option after softmax
  p1 <- (1 - epsilon) * p1 + epsilon * 0.5 # mixes the model probabilities with pure randomness (allows for some randomness in choice selection)
  return(c(1 - p1, p1)) # return probabilities for each choice
}

# ONLY WHEN DEBUGGING

penalty_on_zbeta <- function(z_beta, sigma = 3) {
  if (!USE_BETA_PENALTY) return(0)
  if (!is.finite(z_beta)) return(1e10)
  0.5 * (z_beta^2) / (sigma^2)
}

# REJECTION OF EXTREME PARAMETER VALUE FUNCTION

reject_extreme_params <- function(z_pars) {
  if (any(!is.finite(z_pars))) return(TRUE)
  if (length(z_pars) >= 2) {
    z_beta <- z_pars[length(z_pars)]
    if (!is.finite(z_beta) || abs(z_beta) > MAX_Z_BETA) return(TRUE)
  }
  return(FALSE)
}

# DEFINE "RANDOM" SIDE-BIAS BERNOULLI MODEL (1 free parameter) FUNCTION

nll_random <- function(choice){
  p <- mean(choice, na.rm = TRUE) # estimate probability of right choice
  eps <- 1e-12 
  nll <- -sum(choice * log(p + eps) + (1 - choice) * log(1 - p + eps), na.rm = TRUE)
  return(list(nll = nll, p_right = p))
}



# DEFINE SINGLE ALPHA Q-LEARNING MODEL (2 free parameters) FUNCTION

nll_rw1 <- function(params, choice, reward, return_both = FALSE, epsilon = 0){
  if (reject_extreme_params(params)) return(.Machine$double.xmax/10)
  z_alpha <- params[1]; z_beta <- params[2]
  if (!is.finite(z_alpha) || !is.finite(z_beta)) return(.Machine$double.xmax/10)
  if (abs(z_alpha) > MAX_Z_ALPHA || abs(z_beta) > MAX_Z_BETA) return(.Machine$double.xmax/10)
  alpha <- plogis(z_alpha); beta <- exp(z_beta) # makes beta positive and alpha between 0 and 1
  out <- single_alpha_Q(choice, reward, alpha)
  Q <- out$Q; n <- length(choice); nll <- 0; eps <- 1e-12 # extract q-values trial-by-trial and start NLL accumulation
  for (t in 1:n){
    probs <- softmax_probs(Q[t,], beta, epsilon = epsilon)
    p <- probs[as.integer(choice[t]) + 1]; p <- pmin(pmax(p, eps), 1 - eps)
    nll <- nll - log(p)
  } # for each trial, compute the predicted choice probabilities, pick based on actual choice and add to total NLL
  nll_reg <- nll + penalty_on_zbeta(z_beta, sigma = 3)
  if (return_both) return(list(nll = nll, nll_reg = nll_reg))
  nll_reg
} # when penalty off, just nll 

# DFINE DUAL ALPHA Q-LEARNING MODEL (3 free parameters) FUNCTION

nll_rw2 <- function(params, choice, reward, return_both = FALSE, epsilon = 0){
  if (reject_extreme_params(params)) return(.Machine$double.xmax/10)
  if (any(!is.finite(params))) return(.Machine$double.xmax/10)
  if (abs(params[1]) > MAX_Z_ALPHA || abs(params[2]) > MAX_Z_ALPHA || abs(params[3]) > MAX_Z_BETA) return(.Machine$double.xmax/10)
  alpha_pos <- plogis(params[1]); alpha_neg <- plogis(params[2]); beta <- exp(params[3])
  out <- dual_alpha_Q(choice, reward, alpha_pos, alpha_neg)
  Q <- out$Q; n <- length(choice); nll <- 0; eps <- 1e-12
  for (t in 1:n){
    probs <- softmax_probs(Q[t,], beta, epsilon = epsilon)
    p <- probs[as.integer(choice[t]) + 1]; p <- pmin(pmax(p, eps), 1 - eps)
    nll <- nll - log(p)
  }
  nll_reg <- nll + penalty_on_zbeta(params[3], sigma = 3)
  if (return_both) return(list(nll = nll, nll_reg = nll_reg))
  nll_reg
}

# DEFINE SINGLE ALPHA Q-LEARNING MODEL WITH STICKINESS (4 free parameters) FUNCTION

nll_rw_stick <- function(params, choice, reward, return_both = FALSE, epsilon = 0){
  if (reject_extreme_params(params)) return(.Machine$double.xmax/10)
  if (any(!is.finite(params))) return(.Machine$double.xmax/10)
  if (abs(params[1]) > MAX_Z_ALPHA || abs(params[2]) > MAX_KAPPA || abs(params[3]) > MAX_KAPPA || abs(params[4]) > MAX_Z_BETA) return(.Machine$double.xmax/10)
  
  alpha  <- plogis(params[1]) # bounded between 0 and 1
  k0     <- params[2] # raw
  k1     <- params[3] # raw
  beta   <- exp(params[4]) # force positive
  
  rr <- ema_reward_rate(reward, lambda = EMA_LAMBDA, init = EMA_INIT) # compute recent reward rate estimate
  out <- single_alpha_Q(choice, reward, alpha)
  Q <- out$Q; n <- length(choice); nll <- 0; eps <- 1e-12
  
  for (t in 1:n){
    prev_choice <- if (t == 1) NULL else choice[t-1]
    kappa_t <- k0 + k1 * (rr[t] - RR_CENTER)
    probs <- softmax_probs(Q[t,], beta, prev_choice = prev_choice, kappa = kappa_t, epsilon = epsilon)
    p <- probs[as.integer(choice[t]) + 1]; p <- pmin(pmax(p, eps), 1 - eps)
    nll <- nll - log(p)
  }
  nll_reg <- nll + penalty_on_zbeta(params[4], sigma = 3)
  if (return_both) return(list(nll = nll, nll_reg = nll_reg))
  nll_reg
}

# DEFINE DUAL ALPHA Q-LEARNING MODEL WITH STICKINESS (5 free parameters) FUNCTION

nll_rw2_stick <- function(params, choice, reward, return_both = FALSE, epsilon = 0){
  if (reject_extreme_params(params)) return(.Machine$double.xmax/10)
  if (any(!is.finite(params))) return(.Machine$double.xmax/10)
  if (abs(params[1]) > MAX_Z_ALPHA || abs(params[2]) > MAX_Z_ALPHA ||
      abs(params[3]) > MAX_KAPPA || abs(params[4]) > MAX_KAPPA || abs(params[5]) > MAX_Z_BETA) {
    return(.Machine$double.xmax/10)
  }
  
  alpha_pos <- plogis(params[1])
  alpha_neg <- plogis(params[2])
  k0        <- params[3]
  k1        <- params[4]
  beta      <- exp(params[5])
  
  rr <- ema_reward_rate(reward, lambda = EMA_LAMBDA, init = EMA_INIT)
  out <- dual_alpha_Q(choice, reward, alpha_pos, alpha_neg)
  Q <- out$Q; n <- length(choice); nll <- 0; eps <- 1e-12
  
  for (t in 1:n){
    prev_choice <- if (t == 1) NULL else choice[t-1]
    kappa_t <- k0 + k1 * (rr[t] - RR_CENTER)
    probs <- softmax_probs(Q[t,], beta, prev_choice = prev_choice, kappa = kappa_t, epsilon = epsilon)
    p <- probs[as.integer(choice[t]) + 1]; p <- pmin(pmax(p, eps), 1 - eps)
    nll <- nll - log(p)
  }
  nll_reg <- nll + penalty_on_zbeta(params[5], sigma = 3)
  if (return_both) return(list(nll = nll, nll_reg = nll_reg))
  nll_reg
}

# DEFINE WIN-STAY/LOSE-SHIFT MODEL (2 free parameters) FUNCTION

nll_wsls <- function(params, choice, reward) {
  if (any(!is.finite(params))) return(.Machine$double.xmax/10)
  z_pw <- params[1]; z_pl <- params[2]
  if (abs(z_pw) > 10 || abs(z_pl) > 10) return(.Machine$double.xmax/10)
  p_win_stay <- plogis(z_pw) # bounded between 0 and 1
  p_loss_switch <- plogis(z_pl) # bounded between 0 and 1 
  
  n <- length(choice)
  if (n < 1) return(.Machine$double.xmax/10)
  eps <- 1e-12
  nll <- 0
  
  p_first <- 0.5 # handles first trial (with no trial history)
  nll <- nll - log(pmin(pmax(ifelse(choice[1] == 1, p_first, 1 - p_first), eps), 1 - eps))
  
  if (n >= 2) {
    for (t in 2:n) {
      prev_choice <- choice[t-1]
      prev_reward <- reward[t-1]
      if (!is.na(prev_choice) && !is.na(prev_reward)) {
        if (prev_reward > 0) {
          p_curr <- if (choice[t] == prev_choice) p_win_stay else (1 - p_win_stay)
        } else {
          p_curr <- if (choice[t] != prev_choice) p_loss_switch else (1 - p_loss_switch)
        }
      } else {
        p_curr <- 0.5
      }
      p_curr <- pmin(pmax(p_curr, eps), 1 - eps)
      nll <- nll - log(p_curr)
    }
  }
  nll
}

# GENERATE INFORMATION CRITERION FROM NLL

ic_from_nll <- function(nll, k, n){
  aic <- 2*k + 2*nll # AIC to punish per parameter
  bic <- k*log(n) + 2*nll # BIC to punish proportional to data size 
  return(list(aic = aic, bic = bic))
}

# DEFINE HIDDEN MARKOV MODEL FIT TO CHOICES ALONE (complex free parameters) FUNCTION

fit_hmm_segment <- function(choice, nstates = 2){ # assumes hidden two-states driving choices
  n <- length(choice)
  if (n < MIN_SEG_LEN) return(list(success = FALSE, message = "too_short"))
  dat <- data.frame(ch = factor(choice))
  mod <- depmix(ch ~ 1, data = dat, nstates = nstates, family = multinomial("identity"))
  fm <- try(suppressWarnings(fit(mod, verbose = FALSE)), silent = TRUE)
  if (inherits(fm, "try-error")) return(list(success = FALSE, message = as.character(fm))) # prevents code failing if fitting fails, returns failure message instead
  posterior_df <- posterior(fm) # approximate which state for each trial
  state_p <- tapply(choice, posterior_df$state, mean) # average choice within each state
  state_seq <- posterior_df$state
  p_switch_emp <- sum(diff(state_seq) != 0) / (length(state_seq) - 1) # approximate state switching frequency
  
  ll <- as.numeric(logLik(fm))
  k <- tryCatch(npar(fm), error = function(e) NA_integer_)
  ic <- if (is.finite(ll) && is.finite(k)) ic_from_nll(-ll, k, n) else list(aic = NA_real_, bic = NA_real_)
  nll <- -ll
  
  return(list(success = TRUE,
              p_switch = p_switch_emp,
              state_choice_prob = state_p,
              logLik = ll,
              nll = nll,
              npar = k,
              aic = ic$aic,
              bic = ic$bic,
              model = fm))
}

# CREATE GENERAL PURPOSE WRAPPER FOR THE OPTIMISER

run_safe_optim <- function(start, fn, lower = NULL, upper = NULL, control = list(maxit = 2000), ...){
  p <- length(start)
  if (is.null(lower)) lower <- rep(-50, p)
  if (is.null(upper)) upper <- rep(50, p)
  res <- tryCatch({
    optim(par = start, fn = function(pv) { # prevent errors from crashing the entire code, log instead 
      val <- tryCatch(fn(pv, ...), error = function(e) .Machine$double.xmax/10)
      if (!is.finite(val)) val <- .Machine$double.xmax/10
      val
    }, method = "L-BFGS-B", lower = lower, upper = upper, control = control) # L-BFGS-B optimisation algorithm allows me to constrain parameter searches for the optimiser
  }, error = function(e) NULL)
  if (is.null(res)) return(list(success = FALSE, message = "optim call failed", result = NULL)) # if outer call fails, log
  if (!is.null(res$convergence) && res$convergence != 0) { # if convergence failed or suspect, log as failure 
    return(list(success = FALSE, message = paste0("convergence=", res$convergence), result = res))
  }
  if (!is.finite(res$value)) return(list(success = FALSE, message = "non-finite objective", result = res)) # if final value invalid, reject fit
  return(list(success = TRUE, result = res)) # otherwise, accept and return output 
}

# CREATE GENERATORS TO GIVE OPTIMISER MULTIPLE START POINTS 

make_jittered_starts <- function(base_starts, n_jitter = 8, jitter_sd = 0.6, lower, upper) { # add noise to existing starting points
  stopifnot(is.matrix(base_starts))
  p <- ncol(base_starts)
  out <- base_starts
  for (k in seq_len(n_jitter)) {
    idx <- sample.int(nrow(base_starts), 1)
    s <- base_starts[idx, ] + rnorm(p, mean = 0, sd = jitter_sd)
    s <- pmin(pmax(s, lower), upper)
    out <- rbind(out, s)
  }
  out
}

make_random_starts <- function(n_rand, lower, upper) { # create random starting points within bounds
  p <- length(lower)
  mat <- matrix(NA_real_, nrow = n_rand, ncol = p)
  for (j in seq_len(p)) mat[, j] <- runif(n_rand, min = lower[j], max = upper[j])
  mat
}

# CLEAN UP DATA SEGMENTS (if choice or reward are missing)

clean_segment <- function(choice, reward, min_len = MIN_SEG_LEN){
  ok <- which(!is.na(choice) & !is.na(reward))
  if (length(ok) < min_len) return(list(success = FALSE, reason = "too_short_after_NA_filter"))
  choice <- as.numeric(choice[ok])
  reward <- as.numeric(reward[ok])
  return(list(success = TRUE, choice = choice, reward = reward))
}

# DEFINE FUNCTIONS TO FIT MODELS TO EACH DATA SEGMENT

fit_models_on_segment <- function(seg, segid = NULL){
  out <- list() # create an empty output list for that segment 
  cl <- clean_segment(seg$choice, seg$reward, min_len = MIN_SEG_LEN) # clean that data segment and ensure enough data remains
  if (!cl$success) {
    out$random <- list(nll = NA, p_right = NA, k = 1, aic = NA, bic = NA, fit_status = 'too_short_or_NA')
    out$rw1 <- list(status = 'skipped')
    out$rw2 <- list(status = 'skipped')
    out$rw_stick <- list(status = 'skipped')
    out$rw2_stick <- list(status = 'skipped')
    out$wsls <- list(status = 'skipped')
    out$hmm <- list(success = FALSE, message = 'too_short_or_NA')
    return(out)
  } # error message if data too short
  
  choice <- cl$choice; reward <- cl$reward; n <- length(choice) # extract variables from cleaned data
  
  log_msg("  Segment ", segid, " has n=", n, " valid trials") # print segment info after cleaning
  
  # RANDOM MODEL
  out$random <- time_block("  Fit RANDOM (side-bias coin)", {
    rr <- nll_random(choice)
    list(nll = rr$nll, p_right = rr$p_right, k = 1,
         aic = ic_from_nll(rr$nll, 1, n)$aic,
         bic = ic_from_nll(rr$nll, 1, n)$bic,
         fit_status = 'ok')
  })
  
  epsilon_lapse <- 0.0 # no extra allowance for random choice (beyond softmax)
  
  # SINGLE ALPHA Q-LEARNING (rw1)
  out$rw1 <- time_block("  Fit RW1 (q-learning)", {
    best_rw1 <- tryCatch({
      lower <- c(-MAX_Z_ALPHA, -MAX_Z_BETA)
      upper <- c( MAX_Z_ALPHA,  MAX_Z_BETA)
      starts0 <- rbind(c(qlogis(0.1), log(1)), # sensible guess 1, alpha = 0.1 & beta = 1
                       c(qlogis(0.3), log(2))) # sensible guess 2, alpha = 0.3 & beta = 2
      starts <- make_jittered_starts(starts0, n_jitter = 10, jitter_sd = 0.6, lower = lower, upper = upper) # creates nearby variations to sensible guesses for optimiser
      starts <- rbind(starts, make_random_starts(10, lower, upper)) # adds random starts 
      
      best_val <- Inf; best_res <- NULL
      for (i in 1:nrow(starts)){ # for each start, run the optimiser and keep only the best result
        attempt <- run_safe_optim(starts[i,], nll_rw1,
                                  lower = lower, upper = upper,
                                  control = list(maxit = 2000),
                                  choice = choice, reward = reward,
                                  return_both = FALSE, epsilon = epsilon_lapse)
        if (!attempt$success) { # if one attempt fails, move on
          log_failure(segid, "rw1", attempt$message, start_par = starts[i,], result_par = if(!is.null(attempt$result)) attempt$result$par else NULL)
          next
        }
        res <- attempt$result
        if (res$value < best_val) { best_val <- res$value; best_res <- res } # keep the best fit
      }
      if (is.null(best_res)) stop("optim failed RW1") # only fail if all starts fail
      nlls <- nll_rw1(best_res$par, choice, reward, return_both = TRUE, epsilon = epsilon_lapse)
      alpha <- plogis(best_res$par[1]); beta <- exp(best_res$par[2])
      nll <- nlls$nll; ic <- ic_from_nll(nll, 2, n) # 2 free parameters 
      list(alpha = alpha, beta = beta, nll = nll, aic = ic$aic, bic = ic$bic, status = 'ok', result = best_res)
    }, error = function(e) {
      log_failure(segid, "rw1", as.character(e))
      list(status = 'error', message = as.character(e))
    })
    best_rw1
  })
  
  # DUAL ALPHA Q-LEARNING (rw2)
  out$rw2 <- time_block("  Fit RW2 (asymmetric q-learning)", {
    best_rw2 <- tryCatch({
      lower <- c(-MAX_Z_ALPHA, -MAX_Z_ALPHA, -MAX_Z_BETA)
      upper <- c( MAX_Z_ALPHA,  MAX_Z_ALPHA,  MAX_Z_BETA)
      starts0 <- rbind(c(qlogis(0.1), qlogis(0.05), log(1)), # sensible guess 1, alpha_pos = 0.1, alpha_neg = 0.05, beta = 1
                       c(qlogis(0.2), qlogis(0.1),  log(2))) # sensible guess 2, alpha_pos = 0.2, alpha_neg = 0.1, beta = 2
      starts <- make_jittered_starts(starts0, n_jitter = 12, jitter_sd = 0.6, lower = lower, upper = upper)
      starts <- rbind(starts, make_random_starts(12, lower, upper))
      
      best_val <- Inf; best_res <- NULL
      for (i in 1:nrow(starts)){
        attempt <- run_safe_optim(starts[i,], nll_rw2,
                                  lower = lower, upper = upper,
                                  control = list(maxit = 2000),
                                  choice = choice, reward = reward,
                                  return_both = FALSE, epsilon = epsilon_lapse)
        if (!attempt$success) {
          log_failure(segid, "rw2", attempt$message, start_par = starts[i,], result_par = if(!is.null(attempt$result)) attempt$result$par else NULL)
          next
        }
        res <- attempt$result
        if (res$value < best_val) { best_val <- res$value; best_res <- res }
      }
      if (is.null(best_res)) stop("optim failed RW2")
      nlls <- nll_rw2(best_res$par, choice, reward, return_both = TRUE, epsilon = epsilon_lapse)
      alpha_pos <- plogis(best_res$par[1]); alpha_neg <- plogis(best_res$par[2]); beta <- exp(best_res$par[3])
      nll <- nlls$nll; ic <- ic_from_nll(nll, 3, n) # 3 free parameters 
      list(alpha_pos = alpha_pos, alpha_neg = alpha_neg, beta = beta, nll = nll, aic = ic$aic, bic = ic$bic, status = 'ok', result = best_res)
    }, error = function(e) {
      log_failure(segid, "rw2", as.character(e))
      list(status = 'error', message = as.character(e))
    })
    best_rw2
  })
  
  # SINGLE ALPHA Q-LEARNING WITH STICKINESS (rw_stick)
  out$rw_stick <- time_block("  Fit RW_stick (reward-rate–modulated stickiness)", {
    best_rw_stick <- tryCatch({
      lower <- c(-MAX_Z_ALPHA, -MAX_KAPPA, -MAX_KAPPA, -MAX_Z_BETA)
      upper <- c( MAX_Z_ALPHA,  MAX_KAPPA,  MAX_KAPPA,  MAX_Z_BETA)
      starts0 <- rbind(c(qlogis(0.1), 0.0, 0.0, log(1)), # sensible guess 1, alpha = 0.1, kappa0 = 0, kappa1 = 0, & beta = 1
                       c(qlogis(0.2), 0.5, 0.5, log(2))) # sensible guess 2, alpha = 0.2, kappa0 = 0.5, kappa1 = 0.5, & beta = 2
      starts <- make_jittered_starts(starts0, n_jitter = 14, jitter_sd = 0.7, lower = lower, upper = upper)
      starts <- rbind(starts, make_random_starts(14, lower, upper))
      
      best_val <- Inf; best_res <- NULL
      for (i in 1:nrow(starts)){
        attempt <- run_safe_optim(starts[i,], nll_rw_stick,
                                  lower = lower, upper = upper,
                                  control = list(maxit = 2000),
                                  choice = choice, reward = reward,
                                  return_both = FALSE, epsilon = epsilon_lapse)
        if (!attempt$success) {
          log_failure(segid, "rw_stick", attempt$message, start_par = starts[i,], result_par = if(!is.null(attempt$result)) attempt$result$par else NULL)
          next
        }
        res <- attempt$result
        if (res$value < best_val) { best_val <- res$value; best_res <- res }
      }
      if (is.null(best_res)) stop("optim failed RW_stick")
      nlls <- nll_rw_stick(best_res$par, choice, reward, return_both = TRUE, epsilon = epsilon_lapse)
      alpha <- plogis(best_res$par[1]); k0 <- best_res$par[2]; k1 <- best_res$par[3]; beta <- exp(best_res$par[4])
      nll <- nlls$nll; ic <- ic_from_nll(nll, 4, n) # 4 free parameters 
      list(alpha = alpha, kappa0 = k0, kappa1 = k1, beta = beta, nll = nll, aic = ic$aic, bic = ic$bic, status = 'ok', result = best_res)
    }, error = function(e) {
      log_failure(segid, "rw_stick", as.character(e))
      list(status = 'error', message = as.character(e))
    })
    best_rw_stick
  })
  
  # DUAL ALPHA Q-LEARNING WITH STICKINESS (rw2_stick)
  out$rw2_stick <- time_block("  Fit RW2_stick (dual-alpha + modulated stickiness)", {
    best_rw2_stick <- tryCatch({
      lower <- c(-MAX_Z_ALPHA, -MAX_Z_ALPHA, -MAX_KAPPA, -MAX_KAPPA, -MAX_Z_BETA)
      upper <- c( MAX_Z_ALPHA,  MAX_Z_ALPHA,  MAX_KAPPA,  MAX_KAPPA,  MAX_Z_BETA)
      starts0 <- rbind(c(qlogis(0.1), qlogis(0.05), 0.0, 0.0, log(1)), # sensible guess 1, alpha_pos = 0.1, alpha_neg = 0.05, kappa0 = 0, kappa1 = 0, & beta = 1
                       c(qlogis(0.2), qlogis(0.1),  0.5, 0.5, log(2))) # sensible guess 2, alpha_pos = 0.2, alpha_neg = 0.1, kappa0 = 0.5, kappa1 = 0.5, & beta = 2
      starts <- make_jittered_starts(starts0, n_jitter = 16, jitter_sd = 0.7, lower = lower, upper = upper)
      starts <- rbind(starts, make_random_starts(16, lower, upper))
      
      best_val <- Inf; best_res <- NULL
      for (i in 1:nrow(starts)){
        attempt <- run_safe_optim(starts[i,], nll_rw2_stick,
                                  lower = lower, upper = upper,
                                  control = list(maxit = 2000),
                                  choice = choice, reward = reward,
                                  return_both = FALSE, epsilon = epsilon_lapse)
        if (!attempt$success) {
          log_failure(segid, "rw2_stick", attempt$message, start_par = starts[i,], result_par = if(!is.null(attempt$result)) attempt$result$par else NULL)
          next
        }
        res <- attempt$result
        if (res$value < best_val) { best_val <- res$value; best_res <- res }
      }
      if (is.null(best_res)) stop("optim failed RW2_stick")
      nlls <- nll_rw2_stick(best_res$par, choice, reward, return_both = TRUE, epsilon = epsilon_lapse)
      alpha_pos <- plogis(best_res$par[1]); alpha_neg <- plogis(best_res$par[2])
      k0 <- best_res$par[3]; k1 <- best_res$par[4]; beta <- exp(best_res$par[5])
      nll <- nlls$nll; ic <- ic_from_nll(nll, 5, n) # 5 free parameters 
      list(alpha_pos = alpha_pos, alpha_neg = alpha_neg, kappa0 = k0, kappa1 = k1, beta = beta, nll = nll, aic = ic$aic, bic = ic$bic, status = 'ok', result = best_res)
    }, error = function(e) {
      log_failure(segid, "rw2_stick", as.character(e))
      list(status = 'error', message = as.character(e))
    })
    best_rw2_stick
  })
  
  out$q_learning <- out$rw1 # changed name to q_learning but didn't want to adjust whole code...sorry...
  
  # WSLS
  out$wsls <- time_block("  Fit WSLS", {
    best_wsls <- tryCatch({
      lower <- c(-10, -10) # maps near 0 after qlogis
      upper <- c( 10,  10) # maps near 1 after qlogis
      starts0 <- rbind(c(qlogis(0.7), qlogis(0.7)), # sensible guess 1, win_stay = 0.7, lose_shift = 0.7
                       c(qlogis(0.6), qlogis(0.8))) # sensible guess 2, win_stay = 0.6, lose_shift = 0.8
      starts <- make_jittered_starts(starts0, n_jitter = 10, jitter_sd = 0.8, lower = lower, upper = upper)
      starts <- rbind(starts, make_random_starts(10, lower, upper))
      
      best_val <- Inf; best_res <- NULL
      for (i in 1:nrow(starts)) {
        attempt <- run_safe_optim(starts[i,], nll_wsls,
                                  lower = lower, upper = upper,
                                  control = list(maxit = 2000),
                                  choice = choice, reward = reward)
        if (!attempt$success) {
          log_failure(segid, "wsls", attempt$message, start_par = starts[i,], result_par = if(!is.null(attempt$result)) attempt$result$par else NULL)
          next
        }
        res <- attempt$result
        if (res$value < best_val) { best_val <- res$value; best_res <- res }
      }
      if (is.null(best_res)) stop("optim failed WSLS")
      p_win_stay <- plogis(best_res$par[1]); p_loss_switch <- plogis(best_res$par[2])
      nll <- nll_wsls(best_res$par, choice, reward)
      ic <- ic_from_nll(nll, 2, n)
      list(p_win_stay = p_win_stay, p_loss_switch = p_loss_switch, nll = nll, aic = ic$aic, bic = ic$bic, status = 'ok', result = best_res)
    }, error = function(e) {
      log_failure(segid, "wsls", as.character(e))
      list(status = 'error', message = as.character(e))
    })
    best_wsls
  })
  
  # HIDDEN MARKOV MODEL FIT TO CHOICES ONLY
  out$hmm <- time_block("  Fit HMM (depmixS4)", {
    hmm_res <- tryCatch({
      fit_hmm_segment(choice)
    }, error = function(e) {
      log_failure(segid, "hmm", as.character(e))
      list(success = FALSE, message = as.character(e))
    })
    hmm_res
  })
  
  out
}

# MODELLING PIPELINE 

dat_fp <- "DRL-PRL-Model-Input.csv" # data file with this name must exist in working directory!!
if (!file.exists(dat_fp)) stop("Datafile not found: DRL-PRL-Model-Input.csv")
dt <- fread(dat_fp)

possible_keys <- c("mouse_id","session_type","session_id","group","block_prob") # some input files contained different keys but were cleaned - this is a relic
group_keys <- intersect(possible_keys, names(dt))
if (length(group_keys) == 0) {
  group_keys <- intersect(c("mouse_id","session_id","session_type"), names(dt))
}
if (length(group_keys) == 0) {
  stop("No grouping keys found. Please ensure columns like mouse_id/session_id exist.")
}

# Create seg_id for grouping
dt[, seg_id := do.call(paste, c(.SD, sep = "__")), .SDcols = group_keys]
segments <- unique(dt$seg_id)

log_msg("Detected ", length(segments), " segments (using keys: ", paste(group_keys, collapse = ", "), ")")

# PARALLELISE MODELLING (FOR MULTI-CORE PROCESSING ON UNIMELB RESEARCH COMPUTING)

n_cores <- parallel::detectCores()
n_workers <- min(28, n_cores - 2)  # leave headroom
if (!is.finite(n_workers) || n_workers < 1) n_workers <- 1
log_msg("parallel::detectCores() = ", n_cores, " | Using n_workers = ", n_workers)

# Create cluster
cl <- makeCluster(n_workers, type = "PSOCK")

# Load packages on each worker
clusterEvalQ(cl, {
  library(data.table)
  library(depmixS4)
  library(stats)
  options(stringsAsFactors = FALSE)
  NULL
})

# Export required objects/functions to workers
to_export <- c(
  # data / keys
  "dt","segments","group_keys",
  # logging helpers
  "TS","log_msg","time_block",
  # constants
  "MAX_Z_BETA","MAX_Z_ALPHA","MAX_KAPPA","MIN_SEG_LEN",
  "USE_BETA_PENALTY","EMA_LAMBDA","EMA_INIT","RR_CENTER",
  # debug + logger
  "debug_failures","log_failure",
  # functions
  "ema_reward_rate",
  "dual_alpha_Q","single_alpha_Q",
  "softmax_probs",
  "penalty_on_zbeta","reject_extreme_params",
  "nll_random","nll_rw1","nll_rw2","nll_rw_stick","nll_rw2_stick","nll_wsls",
  "ic_from_nll","fit_hmm_segment",
  "run_safe_optim",
  "make_jittered_starts","make_random_starts",
  "clean_segment",
  "fit_models_on_segment"
)
clusterExport(cl, varlist = to_export, envir = environment())

global_start <- Sys.time() # store start time

# Use load-balanced parallel apply to keep all workers busy
log_msg("Starting parallel fitting across segments...")

# Each worker identifies a unique segment and starts a segment specific timer
results <- parLapplyLB(cl, seq_along(segments), function(i) { 
  segid <- segments[i]
  seg_start <- Sys.time()
  
  # reset worker-local debug table for each new segment
  debug_failures <<- data.table(
    seg_id     = character(),
    model      = character(),
    message    = character(),
    start_par  = list(),
    result_par = list()
  )
  
  log_msg(sprintf("Starting segment %d / %d: %s", i, length(segments), segid))
  
  seg_dt <- dt[seg_id == segid]
  if ("trial" %in% names(seg_dt)) setorder(seg_dt, trial)
  if (!("choice" %in% names(seg_dt) && "reward" %in% names(seg_dt))) {
    warning("Segment ", segid, " missing choice/reward; skipping")
    return(list(i=i, segid=segid, out=NULL, group_keys=NULL, debug_failures=debug_failures, elapsed=NA_real_))
  }
  
  seg_dt[, choice := as.numeric(choice)]
  seg_dt[, reward := as.numeric(reward)]
  
  out <- fit_models_on_segment(seg_dt, segid = segid) # fit all models on this segment 
  
  gk <- seg_dt[1, ..group_keys]
  seg_elapsed <- as.numeric(difftime(Sys.time(), seg_start, units = "secs"))
  
  log_msg(sprintf("Finished segment %d / %d (elapsed %.2f sec): %s", i, length(segments), seg_elapsed, segid))
  
  list(i=i, segid=segid, out=out, group_keys=gk, debug_failures=debug_failures, elapsed=seg_elapsed) # return modelling outputs for that segment
})

stopCluster(cl)
log_msg("Parallel fitting finished.")

total_elapsed <- as.numeric(difftime(Sys.time(), global_start, units = "secs"))
log_msg(sprintf("Total elapsed time: %.2f minutes", total_elapsed/60))

# Reconstruct results in the correct order and merge debug failures across segments
res_list <- vector("list", length(segments))
names(res_list) <- segments

debug_failures_all <- data.table(
  seg_id     = character(),
  model      = character(),
  message    = character(),
  start_par  = list(),
  result_par = list()
)

for (r in results) {
  if (is.null(r$out)) next
  attr(r$out, "group_keys") <- r$group_keys
  res_list[[r$i]] <- r$out
  if (!is.null(r$debug_failures) && nrow(r$debug_failures) > 0) {
    debug_failures_all <- rbindlist(list(debug_failures_all, r$debug_failures), use.names = TRUE, fill = TRUE)
  }
}

# Replace global debug_failures with merged parallel output
debug_failures <- debug_failures_all

# COMBINE RESULTS TABLES INTO ONE

rows <- list()
for (i in seq_along(res_list)){
  segid <- names(res_list)[i]; out <- res_list[[i]]
  if (is.null(out)) next
  
  gvals <- as.list(attr(out, "group_keys"))
  base <- data.table(seg_id = segid, i = i)
  if (length(gvals) > 0) for (nm in names(gvals)) base[[nm]] <- gvals[[nm]]
  
  # random
  if (!is.null(out$random)) {
    base_rand <- copy(base); base_rand[, model := "random"]
    base_rand[, nll := out$random$nll]; base_rand[, aic := out$random$aic]; base_rand[, bic := out$random$bic]
    base_rand[, p_right := out$random$p_right]; base_rand[, fit_status := out$random$fit_status]
    rows[[length(rows) + 1]] <- base_rand
  }
  
  # single alpha q-learning (rw1)
  if (!is.null(out$rw1)) {
    base_rw1 <- copy(base); base_rw1[, model := "q-learning"]
    if (!is.null(out$rw1$status) && out$rw1$status == "ok") {
      base_rw1[, alpha := out$rw1$alpha]; base_rw1[, beta := out$rw1$beta]
      base_rw1[, nll := out$rw1$nll]; base_rw1[, aic := out$rw1$aic]; base_rw1[, bic := out$rw1$bic]
      base_rw1[, fit_status := out$rw1$status]
    } else {
      base_rw1[, fit_status := out$rw1$status]; base_rw1[, message := out$rw1$message]
    }
    rows[[length(rows) + 1]] <- base_rw1
  }
  
  # dual alpha q-learning (rw2)
  if (!is.null(out$rw2)) {
    base_rw2 <- copy(base); base_rw2[, model := "asymmetric q-learning"]
    if (!is.null(out$rw2$status) && out$rw2$status == "ok") {
      base_rw2[, alpha_pos := out$rw2$alpha_pos]; base_rw2[, alpha_neg := out$rw2$alpha_neg]; base_rw2[, beta := out$rw2$beta]
      base_rw2[, nll := out$rw2$nll]; base_rw2[, aic := out$rw2$aic]; base_rw2[, bic := out$rw2$bic]
      base_rw2[, fit_status := out$rw2$status]
    } else {
      base_rw2[, fit_status := out$rw2$status]; base_rw2[, message := out$rw2$message]
    }
    rows[[length(rows) + 1]] <- base_rw2
  }
  
  # single alpha q-learning with stickiness (rw_stick)
  if (!is.null(out$rw_stick)) {
    base_rwst <- copy(base); base_rwst[, model := "q-learning with stickiness"]
    if (!is.null(out$rw_stick$status) && out$rw_stick$status == "ok") {
      base_rwst[, alpha := out$rw_stick$alpha]
      base_rwst[, kappa0 := out$rw_stick$kappa0]; base_rwst[, kappa1 := out$rw_stick$kappa1]
      base_rwst[, beta := out$rw_stick$beta]
      base_rwst[, nll := out$rw_stick$nll]; base_rwst[, aic := out$rw_stick$aic]; base_rwst[, bic := out$rw_stick$bic]
      base_rwst[, fit_status := out$rw_stick$status]
    } else {
      base_rwst[, fit_status := out$rw_stick$status]; base_rwst[, message := out$rw_stick$message]
    }
    rows[[length(rows) + 1]] <- base_rwst
  }
  
  # dual alpha q-learning with stickiness (rw2_stick)
  if (!is.null(out$rw2_stick)) {
    base_rw2s <- copy(base); base_rw2s[, model := "q-learning dual-alpha with stickiness"]
    if (!is.null(out$rw2_stick$status) && out$rw2_stick$status == "ok") {
      base_rw2s[, alpha_pos := out$rw2_stick$alpha_pos]; base_rw2s[, alpha_neg := out$rw2_stick$alpha_neg]
      base_rw2s[, kappa0 := out$rw2_stick$kappa0]; base_rw2s[, kappa1 := out$rw2_stick$kappa1]
      base_rw2s[, beta := out$rw2_stick$beta]
      base_rw2s[, nll := out$rw2_stick$nll]; base_rw2s[, aic := out$rw2_stick$aic]; base_rw2s[, bic := out$rw2_stick$bic]
      base_rw2s[, fit_status := out$rw2_stick$status]
    } else {
      base_rw2s[, fit_status := out$rw2_stick$status]; base_rw2s[, message := out$rw2_stick$message]
    }
    rows[[length(rows) + 1]] <- base_rw2s
  }
  
  # win-stay/lose-shift
  if (!is.null(out$wsls)) {
    base_ws <- copy(base); base_ws[, model := "wsls"]
    if (!is.null(out$wsls$status) && out$wsls$status == "ok") {
      base_ws[, p_win_stay := out$wsls$p_win_stay]
      base_ws[, p_loss_switch := out$wsls$p_loss_switch]
      base_ws[, nll := out$wsls$nll]; base_ws[, aic := out$wsls$aic]; base_ws[, bic := out$wsls$bic]
      base_ws[, fit_status := out$wsls$status]
    } else {
      base_ws[, fit_status := out$wsls$status]; base_ws[, message := out$wsls$message]
    }
    rows[[length(rows) + 1]] <- base_ws
  }
  
  # hidden markov model
  if (!is.null(out$hmm)) {
    base_h <- copy(base); base_h[, model := "hmm"]
    if (!is.null(out$hmm$success) && out$hmm$success) {
      base_h[, p_switch := out$hmm$p_switch]
      base_h[, state_choice_prob := paste(names(out$hmm$state_choice_prob), round(unlist(out$hmm$state_choice_prob), 3), collapse = ";")]
      base_h[, logLik := as.numeric(out$hmm$logLik)]
      base_h[, nll := as.numeric(out$hmm$nll)]
      base_h[, npar := as.integer(out$hmm$npar)]
      base_h[, aic := out$hmm$aic]; base_h[, bic := out$hmm$bic]
      base_h[, fit_status := "ok"]
    } else {
      base_h[, fit_status := "error"]; base_h[, message := out$hmm$message]
    }
    rows[[length(rows) + 1]] <- base_h
  }
}

fits_dt <- if (length(rows) > 0) rbindlist(rows, fill = TRUE) else data.table()

# COMPUTE AIC/BIC WINNERS (LOWEST VALUE) PER SEGMENT

REQUIRED_FOR_WINNERS <- c("random", "wsls", "q-learning", "asymmetric q-learning", "hmm") # stickiness dropped to examine subset

if (nrow(fits_dt) > 0) {
  fits_dt[, ok_fit := (!is.na(aic) | !is.na(bic)) & (is.na(fit_status) | fit_status == "ok")]
  
  winners <- fits_dt[, {
    req_ok <- all(REQUIRED_FOR_WINNERS %in% model[ok_fit])
    if (!req_ok) {
      list(aic_winner = NA_character_, bic_winner = NA_character_, n_ok = sum(ok_fit), req_ok = FALSE)
    } else {
      aic_ok <- .SD[ok_fit & is.finite(aic)]
      bic_ok <- .SD[ok_fit & is.finite(bic)]
      aic_w <- if (nrow(aic_ok) > 0) aic_ok$model[which.min(aic_ok$aic)] else NA_character_
      bic_w <- if (nrow(bic_ok) > 0) bic_ok$model[which.min(bic_ok$bic)] else NA_character_
      list(aic_winner = aic_w, bic_winner = bic_w, n_ok = sum(ok_fit), req_ok = TRUE)
    }
  }, by = seg_id]
  
  fits_dt <- merge(fits_dt, winners, by = "seg_id", all.x = TRUE)
}

# EXPORT RESULTS

out_csv <- "model_fits_per_segment_R_exog_ema02_modkappa.csv" # main results stored here 
fwrite(fits_dt, out_csv)
log_msg("Saved model fits summary to ", out_csv)

# Save debug failures (if any)
if (nrow(debug_failures) > 0) {
  df_out <- copy(debug_failures)
  df_out[, start_par := sapply(start_par, function(x) {
    if (is.null(x)) return(NA_character_)
    paste(x, collapse = ",")
  })]
  df_out[, result_par := sapply(result_par, function(x) {
    if (is.null(x)) return(NA_character_)
    paste(x, collapse = ",")
  })]
  fwrite(df_out, "model_fit_debug_fails_exog_ema02_modkappa.csv") # save optimisation failures for debugging / troubleshooting
  log_msg("Saved optimisation debug failures to model_fit_debug_fails_exog_ema02_modkappa.csv")
} else {
  log_msg("No debug failures recorded.")
}

# EXPORT TRIAL-BY-TRIAL RECONSTRUCTIONS FOR FITTED Q-LEARNING MODELS 

recon_dir <- "reconstructions_R_exog_ema02_modkappa"
if (!dir.exists(recon_dir)) dir.create(recon_dir)

for (i in seq_along(res_list)) {
  segid <- names(res_list)[i]; out <- res_list[[i]]
  if (is.null(out)) next
  
  seg_dt_full <- dt[seg_id == segid]
  seg_dt <- seg_dt_full[!is.na(choice) & !is.na(reward)]
  if (nrow(seg_dt) == 0) next
  
  choice <- as.numeric(seg_dt$choice)
  reward <- as.numeric(seg_dt$reward)
  rr <- ema_reward_rate(reward, lambda = EMA_LAMBDA, init = EMA_INIT)
  
  # single alpha q-learning (rw1)
  if (!is.null(out$rw1) && !is.null(out$rw1$status) && out$rw1$status == 'ok') {
    alpha <- out$rw1$alpha; beta <- out$rw1$beta
    sim <- single_alpha_Q(choice, reward, alpha)
    probs <- sapply(1:nrow(sim$Q), function(t) softmax_probs(sim$Q[t,], beta)[2])
    recon <- data.table(seg_id = segid, trial = seq_len(nrow(sim$Q)),
                        q1 = sim$Q[,1], q2 = sim$Q[,2], p_choice1 = probs)
    fwrite(recon, file.path(recon_dir, paste0("recon_q_learning_",
                                              gsub("[^A-Za-z0-9]", "_", segid), ".csv")))
  }
  
  # dual alpha q-learning (rw2)
  if (!is.null(out$rw2) && !is.null(out$rw2$status) && out$rw2$status == 'ok') {
    alpha_pos <- out$rw2$alpha_pos; alpha_neg <- out$rw2$alpha_neg; beta <- out$rw2$beta
    sim <- dual_alpha_Q(choice, reward, alpha_pos, alpha_neg)
    probs <- sapply(1:nrow(sim$Q), function(t) softmax_probs(sim$Q[t,], beta)[2])
    recon <- data.table(seg_id = segid, trial = seq_len(nrow(sim$Q)),
                        q1 = sim$Q[,1], q2 = sim$Q[,2], p_choice1 = probs)
    fwrite(recon, file.path(recon_dir, paste0("recon_asymmetric_q_learning_",
                                              gsub("[^A-Za-z0-9]", "_", segid), ".csv")))
  }
  
  # single alpha q-learning with stickiness (rw_stick)
  if (!is.null(out$rw_stick) && !is.null(out$rw_stick$status) && out$rw_stick$status == 'ok') {
    alpha <- out$rw_stick$alpha
    k0 <- out$rw_stick$kappa0
    k1 <- out$rw_stick$kappa1
    beta <- out$rw_stick$beta
    sim <- single_alpha_Q(choice, reward, alpha)
    probs <- sapply(1:nrow(sim$Q), function(t) {
      prev_choice <- if (t == 1) NULL else choice[t-1]
      kappa_t <- k0 + k1 * (rr[t] - RR_CENTER)
      softmax_probs(sim$Q[t,], beta, prev_choice = prev_choice, kappa = kappa_t)[2]
    })
    recon <- data.table(seg_id = segid, trial = seq_len(nrow(sim$Q)),
                        q1 = sim$Q[,1], q2 = sim$Q[,2], p_choice1 = probs,
                        rr = rr, kappa_t = k0 + k1 * (rr - RR_CENTER))
    fwrite(recon, file.path(recon_dir, paste0("recon_q_learning_with_stickiness_",
                                              gsub("[^A-Za-z0-9]", "_", segid), ".csv")))
  }
  
  # dual alpha q-learning with stickiness (rw2_stick)
  if (!is.null(out$rw2_stick) && !is.null(out$rw2_stick$status) && out$rw2_stick$status == 'ok') {
    alpha_pos <- out$rw2_stick$alpha_pos
    alpha_neg <- out$rw2_stick$alpha_neg
    k0 <- out$rw2_stick$kappa0
    k1 <- out$rw2_stick$kappa1
    beta <- out$rw2_stick$beta
    sim <- dual_alpha_Q(choice, reward, alpha_pos, alpha_neg)
    probs <- sapply(1:nrow(sim$Q), function(t) {
      prev_choice <- if (t == 1) NULL else choice[t-1]
      kappa_t <- k0 + k1 * (rr[t] - RR_CENTER)
      softmax_probs(sim$Q[t,], beta, prev_choice = prev_choice, kappa = kappa_t)[2]
    })
    recon <- data.table(seg_id = segid, trial = seq_len(nrow(sim$Q)),
                        q1 = sim$Q[,1], q2 = sim$Q[,2], p_choice1 = probs,
                        rr = rr, kappa_t = k0 + k1 * (rr - RR_CENTER))
    fwrite(recon, file.path(recon_dir, paste0("recon_q_learning_dual_alpha_with_stickiness_",
                                              gsub("[^A-Za-z0-9]", "_", segid), ".csv")))
  }
}

log_msg("Done. Reconstructions saved to ", recon_dir)
