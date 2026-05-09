# INTSINT (INternational Transcription System for INTonation) — Pure R
# Port of momel_intsint.py lines 512-659 (itself a port of intsint.pl)
#
# Reference: Hirst (2005). Form and function in the representation of
# speech prosody. Speech Communication, 46(3-4), 334-347.

# Constants
MIN_PAUSE  <- 0.5   # seconds
MIN_RANGE  <- 0.5   # octaves
MAX_RANGE  <- 2.5
STEP_RANGE <- 0.1
MEAN_SHIFT <- 50    # Hz
STEP_SHIFT <- 1     # Hz
BIG_NUMBER <- 9999

# Tone estimation parameters
HIGHER <- 0.5; LOWER <- 0.5; UP <- 0.25; DOWN <- 0.25

octave <- function(hz) ifelse(hz > 0, log2(hz), 0)
linear <- function(oct) 2^oct

estimate_target <- function(tone, last_target, top, bottom, mid) {
  switch(tone,
    M = mid,
    S = last_target,
    T = top,
    B = bottom,
    H = last_target + (top - last_target) * HIGHER,
    L = last_target - (last_target - bottom) * LOWER,
    U = last_target + (top - last_target) * UP,
    D = last_target - (last_target - bottom) * DOWN,
    mid
  )
}

optimize_intsint <- function(targets, mid, range_oct) {
  top    <- mid + range_oct / 2
  bottom <- mid - range_oct / 2
  valid  <- targets$frequency > 0
  if (!any(valid)) return(list(tones = character(0), estimates = numeric(0), ss_error = BIG_NUMBER))

  tgt <- targets[valid, , drop = FALSE]
  all_tones <- c("M", "S", "T", "B", "H", "L", "U", "D")

  tones     <- character(nrow(tgt))
  estimates <- numeric(nrow(tgt))
  last_est  <- 0
  ss_error  <- 0
  last_time <- -BIG_NUMBER

  for (i in seq_len(nrow(tgt))) {
    target_oct <- octave(tgt$frequency[i])

    # After pause or first target: choose from {M, T, B}
    if (i == 1 || (tgt$time[i] - last_time > MIN_PAUSE * 1000 / PAS_TRAME)) {
      if (top - target_oct < abs(target_oct - mid)) {
        tone <- "T"
      } else if (target_oct - bottom < abs(target_oct - mid)) {
        tone <- "B"
      } else {
        tone <- "M"
      }
    } else {
      # Minimise distance over all tones except M
      min_diff <- BIG_NUMBER; best <- "S"
      for (t in all_tones) {
        if (t == "M") next
        est <- estimate_target(t, last_est, top, bottom, mid)
        d <- abs(target_oct - est)
        if (d < min_diff) { min_diff <- d; best <- t }
      }
      tone <- best
    }

    tones[i] <- tone
    est <- estimate_target(tone, last_est, top, bottom, mid)
    estimates[i] <- est
    ss_error <- ss_error + (est - target_oct)^2
    last_est  <- est
    last_time <- tgt$time[i]
  }
  list(tones = tones, estimates = estimates, ss_error = ss_error)
}

#' Run INTSINT algorithm on MOMEL targets
#'
#' @param targets data.frame with time (frames) and frequency (Hz) from momel()
#' @return list(targets_df, range, key) where targets_df has columns
#'   time_sec, tone, target_hz, estimate_hz
intsint <- function(targets) {
  valid <- targets$frequency > 0
  if (!any(valid)) return(list(targets_df = data.frame(), range = 0, key = 0))
  tgt <- targets[valid, , drop = FALSE]

  f0_oct   <- octave(tgt$frequency)
  mean_oct <- mean(f0_oct)
  linear_mean <- round(linear(mean_oct))

  min_mean <- linear_mean - MEAN_SHIFT
  max_mean <- linear_mean + MEAN_SHIFT

  best_ss    <- BIG_NUMBER
  best_tones <- character(0)
  best_est   <- numeric(0)
  best_mid   <- mean_oct
  best_range <- MIN_RANGE

  for (range_oct in seq(MIN_RANGE, MAX_RANGE, by = STEP_RANGE)) {
    for (lm in seq(as.integer(min_mean), as.integer(max_mean), by = STEP_SHIFT)) {
      mid <- octave(lm)
      res <- optimize_intsint(tgt, mid, range_oct)
      if (res$ss_error < best_ss) {
        best_ss    <- res$ss_error
        best_range <- range_oct
        best_mid   <- mid
        best_tones <- res$tones
        best_est   <- res$estimates
      }
    }
  }

  out <- data.frame(
    time_sec    = tgt$time * PAS_TRAME / 1000,
    tone        = best_tones,
    target_hz   = tgt$frequency,
    estimate_hz = linear(best_est),
    stringsAsFactors = FALSE
  )
  list(targets_df = out, range = best_range, key = linear(best_mid))
}
