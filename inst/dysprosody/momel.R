# MOMEL (MOdelling MELody) Algorithm — Pure R Implementation
# Port of momel_intsint.py lines 76-510 (itself a port of momel.c/momel.h)
#
# Reference: Hirst & Espesser (1993). Automatic modelling of fundamental
# frequency using a quadratic spline function.

# Constants
PAS_TRAME <- 10.0
SEUILV    <- 50.0
RAPP_GLITCH <- 0.05
HALO_BORNE_TRAME <- 4
FSIGMA <- 1.0

eliminate_glitches <- function(hz) {
  n <- length(hz)
  if (n < 3) return(hz)
  out <- hz
  for (i in 2:(n - 1)) {
    if (hz[i] > hz[i - 1] * (1 + RAPP_GLITCH) &&
        hz[i] > hz[i + 1] * (1 + RAPP_GLITCH)) {
      out[i] <- 0.0
    }
  }
  out
}

# Weighted quadratic regression (0-indexed frame coordinates)
# Returns list(a0, a1, a2, hzes) or NULL
calc_regression <- function(pond, dpx, fpx, hzptr) {
  n <- length(hzptr)
  pn <- 0; sx <- 0; sx2 <- 0; sx3 <- 0; sx4 <- 0
  sy <- 0; sxy <- 0; sx2y <- 0

  for (ix in dpx:fpx) {
    idx <- ix + 1L  # R 1-based
    if (idx < 1 || idx > n) next
    p <- pond[idx]
    if (p == 0) next
    val_ix <- as.double(ix)
    y <- hzptr[idx]
    x2 <- val_ix * val_ix
    pn   <- pn   + p
    sx   <- sx   + p * val_ix
    sx2  <- sx2  + p * x2
    sx3  <- sx3  + p * x2 * val_ix
    sx4  <- sx4  + p * x2 * x2
    sy   <- sy   + p * y
    sxy  <- sxy  + p * val_ix * y
    sx2y <- sx2y + p * x2 * y
  }
  if (pn < 3) return(NULL)

  spdxy  <- sxy  - (sx * sy) / pn
  spdx2  <- sx2  - (sx * sx) / pn
  spdx3  <- sx3  - (sx * sx2) / pn
  spdx4  <- sx4  - (sx2 * sx2) / pn
  spdx2y <- sx2y - (sx2 * sy) / pn

  muet <- spdx2 * spdx4 - spdx3 * spdx3
  if (spdx2 == 0 || muet == 0) return(NULL)

  a2 <- (spdx2y * spdx2 - spdxy * spdx3) / muet
  a1 <- (spdxy - a2 * spdx3) / spdx2
  a0 <- (sy - a1 * sx - a2 * sx2) / pn

  hzes <- numeric(n)
  for (ix in dpx:fpx) {
    idx <- ix + 1L
    if (idx < 1 || idx > n) next
    hzes[idx] <- a0 + (a1 + a2 * ix) * ix
  }
  list(a0 = a0, a1 = a1, a2 = a2, hzes = hzes)
}

cible <- function(nval, hzptr, lfen1, maxec, hzinf, hzsup) {
  pond <- ifelse(hzptr > SEUILV, 1.0, 0.0)
  lfens2 <- lfen1 %/% 2
  cib_time <- numeric(nval)
  cib_freq <- numeric(nval)

  for (ix_r in seq_len(nval)) {
    ix <- ix_r - 1L  # 0-based
    dpx <- ix - lfens2
    fpx <- dpx + lfen1
    nsup <- 0L; nsupr <- -1L
    pondloc <- pond

    while (nsup > nsupr) {
      nsupr <- nsup; nsup <- 0L
      reg <- calc_regression(pondloc, dpx, fpx, hzptr)
      if (is.null(reg)) break
      for (x in dpx:fpx) {
        idx <- x + 1L
        if (idx < 1 || idx > nval) next
        if (hzptr[idx] == 0 || reg$hzes[idx] / hzptr[idx] > maxec) {
          pondloc[idx] <- 0
          nsup <- nsup + 1L
        }
      }
    }

    xc <- 0; yc <- 0
    if (!is.null(reg) && reg$a2 != 0) {
      vxc <- -reg$a1 / (reg$a2 * 2)
      if ((ix - lfen1) < vxc && vxc < (ix + lfen1)) {
        vyc <- reg$a0 + (reg$a1 + reg$a2 * vxc) * vxc
        if (hzinf < vyc && vyc < hzsup) {
          xc <- vxc; yc <- vyc
        }
      }
    }
    cib_time[ix_r] <- xc
    cib_freq[ix_r] <- yc
  }
  data.frame(time = cib_time, frequency = cib_freq)
}

reduc <- function(nval, lfen2, seuildiff_x, seuilrapp_y, cib) {
  valid <- cib$frequency > 0
  if (!any(valid)) return(data.frame(time = numeric(0), frequency = numeric(0)))
  cib_s <- cib[valid, , drop = FALSE]
  cib_s <- cib_s[order(cib_s$time), , drop = FALSE]
  nc <- nrow(cib_s)
  if (nc < 2) return(cib_s)

  lf <- lfen2 %/% 2
  xdist <- numeric(nc - 1)
  ydist <- numeric(nc - 1)

  for (i in seq_len(nc - 1)) {
    j1 <- max(1, i - lf); j2 <- min(nc, i + lf + 1)
    left  <- cib_s[j1:i, , drop = FALSE]
    right <- cib_s[(i + 1):j2, , drop = FALSE]
    sxg <- mean(left$time);  syg <- mean(left$frequency)
    sxd <- mean(right$time); syd <- mean(right$frequency)
    xdist[i] <- abs(sxg - sxd)
    ydist[i] <- abs(syg - syd)
  }

  ok <- xdist > 0
  if (!any(ok)) return(cib_s[1, , drop = FALSE])
  xds <- sum(xdist[ok]); yds <- sum(ydist[ok])
  if (xds == 0 || yds == 0) return(cib_s[1, , drop = FALSE])
  nok <- sum(ok)
  px <- nok / xds; py <- nok / yds

  dist <- ifelse(ok, (xdist * px + ydist * py) / (px + py), -1)
  seuil <- 2.0 / (px + py)

  partitions <- 1L
  susseuil <- FALSE; xmax <- 1L
  for (i in seq_along(dist)) {
    if (!susseuil) {
      if (dist[i] > seuil) { susseuil <- TRUE; xmax <- i }
    } else {
      if (dist[i] > dist[xmax]) xmax <- i
      if (dist[i] < seuil) { partitions <- c(partitions, xmax); susseuil <- FALSE }
    }
  }
  if (susseuil) partitions <- c(partitions, xmax)
  partitions <- c(partitions, nc + 1L)

  cibred_t <- numeric(0); cibred_f <- numeric(0)
  for (ip in seq_len(length(partitions) - 1)) {
    parinf <- partitions[ip]; parsup <- partitions[ip + 1] - 1L
    grp <- cib_s[parinf:parsup, , drop = FALSE]
    n <- nrow(grp)
    if (n > 1) {
      xm <- mean(grp$time); ym <- mean(grp$frequency)
      varx <- max(0.1, mean(grp$time^2) - xm^2)
      vary <- max(0.1, mean(grp$frequency^2) - ym^2)
      et2x <- FSIGMA * sqrt(varx); et2y <- FSIGMA * sqrt(vary)
      keep <- grp$time >= xm - et2x & grp$time <= xm + et2x &
              grp$frequency >= ym - et2y & grp$frequency <= ym + et2y
      if (any(keep)) grp <- grp[keep, , drop = FALSE]
    }
    cibred_t <- c(cibred_t, mean(grp$time))
    cibred_f <- c(cibred_f, mean(grp$frequency))
  }
  if (length(cibred_t) == 0) return(data.frame(time = numeric(0), frequency = numeric(0)))

  # Merge close targets
  out_t <- cibred_t[1]; out_f <- cibred_f[1]
  for (i in seq_along(cibred_t)[-1]) {
    if (cibred_t[i] - out_t[length(out_t)] < seuildiff_x) {
      if (abs(cibred_f[i] - out_f[length(out_f)]) / out_f[length(out_f)] < seuilrapp_y) {
        out_t[length(out_t)] <- (out_t[length(out_t)] + cibred_t[i]) / 2
        out_f[length(out_f)] <- (out_f[length(out_f)] + cibred_f[i]) / 2
      }
    } else {
      out_t <- c(out_t, cibred_t[i])
      out_f <- c(out_f, cibred_f[i])
    }
  }
  data.frame(time = out_t, frequency = out_f)
}

borne <- function(nval, cibred, hzptr) {
  if (nrow(cibred) == 0) return(cibred)
  result_t <- numeric(0); result_f <- numeric(0)

  # First voiced frame (0-based)
  first_voiced <- 0L
  for (i in seq_len(nval)) {
    if (hzptr[i] >= SEUILV) { first_voiced <- i - 1L; break }
  }

  # Left boundary
  if (cibred$time[1] > first_voiced + HALO_BORNE_TRAME) {
    ancre_t <- cibred$time[1]; ancre_f <- cibred$frequency[1]
    sx2y <- 0; sx4 <- 0; j <- 0L
    for (i in as.integer(ancre_t):0) {
      idx <- i + 1L
      if (idx >= 1 && idx <= nval && hzptr[idx] > SEUILV) {
        x2 <- as.double(j * j)
        sx2y <- sx2y + x2 * (hzptr[idx] - ancre_f)
        sx4  <- sx4  + x2 * x2
      }
      j <- j + 1L
    }
    if (sx4 != 0) {
      a <- sx2y / sx4
      frontiere <- as.double(first_voiced)
      bx <- frontiere - (ancre_t - frontiere)
      by <- ancre_f + 2 * a * (ancre_t - frontiere)^2
      result_t <- bx; result_f <- by
    }
  }

  result_t <- c(result_t, cibred$time)
  result_f <- c(result_f, cibred$frequency)

  # Last voiced frame (0-based)
  last_voiced <- nval - 1L
  for (i in nval:1) {
    if (hzptr[i] >= SEUILV) { last_voiced <- i - 1L; break }
  }

  # Right boundary
  nr <- nrow(cibred)
  if (cibred$time[nr] < last_voiced - HALO_BORNE_TRAME) {
    ancre_t <- cibred$time[nr]; ancre_f <- cibred$frequency[nr]
    sx2y <- 0; sx4 <- 0; j <- 0L
    for (i in as.integer(ancre_t):(nval - 1L)) {
      idx <- i + 1L
      if (idx >= 1 && idx <= nval && hzptr[idx] > SEUILV) {
        x2 <- as.double(j * j)
        sx2y <- sx2y + x2 * (hzptr[idx] - ancre_f)
        sx4  <- sx4  + x2 * x2
      }
      j <- j + 1L
    }
    if (sx4 != 0) {
      a <- sx2y / sx4
      frontiere <- as.double(last_voiced)
      bx <- frontiere + (frontiere - ancre_t)
      by <- ancre_f + 2 * a * (ancre_t - frontiere)^2
      result_t <- c(result_t, bx)
      result_f <- c(result_f, by)
    }
  }
  data.frame(time = result_t, frequency = result_f)
}

#' Run MOMEL algorithm on F0 values
#'
#' @param f0_values numeric vector of F0 in Hz (10ms frames), 0 = unvoiced
#' @param window_length window in samples for cible (default 30ms / 10ms = 3)
#' @param min_f0 minimum F0 in Hz
#' @param max_f0 maximum F0 in Hz
#' @param max_error maximum error ratio
#' @param reduced_window_length window for reduction
#' @param minimal_distance minimal distance in frames
#' @param minimal_frequency_ratio minimal frequency ratio for merging
#' @return data.frame with columns time (frames, 0-based) and frequency (Hz)
momel <- function(f0_values,
                  window_length = 30L,
                  min_f0 = 60,
                  max_f0 = 750,
                  max_error = 1.04,
                  reduced_window_length = 20L,
                  minimal_distance = 5.0,
                  minimal_frequency_ratio = 0.05) {
  nval <- length(f0_values)
  hz <- pmin(pmax(f0_values, min_f0), max_f0)
  hz <- eliminate_glitches(hz)
  targets <- cible(nval, hz, window_length, max_error, min_f0, max_f0)
  targets <- reduc(nval, reduced_window_length, minimal_distance,
                   minimal_frequency_ratio, targets)
  borne(nval, targets, hz)
}
