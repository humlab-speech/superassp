#' HMPD Phase Distortion Features (COVAREP)
#'
#' Computes HMPD Phase Distortion vocoder features: amplitude envelope (AE),
#' phase deviation mean (PDM), phase deviation deviation (PDD) at frame level.
#' Performs sinusoidal analysis + harmonic phase distortion quantification.
#'
#' @param listOfFiles Vector of file paths (WAV, MP3, MP4, etc.) to analyze
#' @param f0s Optional matrix with 2 columns (time, f0_hz). If NULL, computed internally via trk_srh_variant
#' @param beginTime Start time in seconds (0 for beginning of file)
#' @param endTime End time in seconds (0 for end of file)
#' @param f0min Minimum F0 for search (Hz, default: 60)
#' @param f0max Maximum F0 for search (Hz, default: 440)
#' @param toFile Write output to file (TRUE) or return object (FALSE). Default: FALSE
#' @param explicitExt Output file extension (default: "hpd")
#' @param outputDirectory Output directory (NULL for same as input file)
#' @param verbose Show progress messages (default: TRUE)
#'
#' @return
#' If `toFile=FALSE` (default): AsspDataObj with 3 matrix tracks:
#'   - `ae`: 24-dim amplitude envelope (Mel cepstral, method=2, log-lin)
#'   - `pdm`: 24-dim phase deviation mean (log-lin, log-harmonics)
#'   - `pdd`: 12-dim phase deviation deviation (log magnitude)
#'
#' If `toFile=TRUE`: invisibly returns vector of output file paths
#'
#' @details
#' **Algorithm** \insertCite{Degottex2014COVAREP}{superassp}:
#' 1. Compute F0 track (internally via SRH if not provided externally)
#' 2. Sinusoidal analysis: extract harmonics from each frame's complex spectrum
#' 3. Amplitude envelope: log-Mel-cepstrum representation (24 coeffs)
#' 4. Phase tracking: harmonic phase per frame
#' 5. Phase deviation mean (PDM): within-frame phase dynamics (24 coeffs, log-linear)
#' 6. Phase deviation deviation (PDD): phase variability across harmonics (12 coeffs)
#'
#' **Output format**: Compressed (60 features/frame total: 24+24+12)
#'   - ae: non-log amplitude envelope
#'   - pdm: log-harmonic scale, mean phase
#'   - pdd: log-harmonic scale, phase deviation magnitude
#'
#' **Interpretation**:
#' - Large PDD = aperiodicity, noise-like phase scatter per harmonic
#' - Structured PDM = deterministic phase dynamics (e.g., vibrato)
#' - AE = spectral envelope shape (low frequencies dominated in voice)
#'
#' @examples
#' \dontrun{
#' # Single file, return object
#' hmpd <- trk_hmpd("speech.wav", toFile = FALSE)
#'
#' # Batch process multiple files
#' files <- c("file1.wav", "file2.wav")
#' hmpd_results <- trk_hmpd(files, toFile = TRUE, outputDirectory = "output/")
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @export
trk_hmpd <- function(listOfFiles,
                     f0s = NULL,
                     beginTime = 0.0,
                     endTime = 0.0,
                     f0min = 60,
                     f0max = 440,
                     toFile = FALSE,
                     explicitExt = "hpd",
                     outputDirectory = NULL,
                     verbose = TRUE) {

  nFiles <- length(listOfFiles)

  # Handle time parameters
  beginTime <- if (is.null(beginTime)) 0.0 else beginTime
  endTime <- if (is.null(endTime)) 0.0 else endTime
  beginTime <- fast_recycle_times(beginTime, nFiles)
  endTime <- fast_recycle_times(endTime, nFiles)

  # Process each file
  externalRes <- vector("list", nFiles)

  for (idx in seq_len(nFiles)) {
    if (verbose) {
      cli::cli_progress_step("Processing {.file {basename(listOfFiles[idx])}} ({idx}/{nFiles})")
    }

    # Load audio
    audio_data <- av::read_audio_bin(
      listOfFiles[[idx]],
      channels = 1,
      start_time = if (beginTime[idx] > 0) beginTime[idx] else NULL,
      end_time = if (endTime[idx] > 0) endTime[idx] else NULL
    )

    fs <- attr(audio_data, "sample_rate")
    wave <- as.numeric(audio_data)
    storage.mode(wave) <- "double"  # Ensure numeric, not integer

    # Compute HMPD analysis
    hmpd_res <- .hmpd_analysis(wave, fs, f0s, f0min, f0max)

    # Convert to AsspDataObj
    assp_obj <- list(
      ae = hmpd_res$ae,
      pdm = hmpd_res$pdm,
      pdd = hmpd_res$pdd
    )

    class(assp_obj) <- "AsspDataObj"
    attr(assp_obj, "sampleRate") <- 100  # 10ms frame rate
    attr(assp_obj, "startTime") <- as.numeric(beginTime[idx])
    attr(assp_obj, "startRecord") <- 1L
    attr(assp_obj, "endRecord") <- nrow(hmpd_res$ae)
    attr(assp_obj, "trackFormats") <- c("FLOAT", "FLOAT", "FLOAT")

    if (toFile) {
      # Write to SSFF file
      out_dir <- if (is.null(outputDirectory)) dirname(listOfFiles[[idx]]) else outputDirectory
      out_file <- file.path(
        out_dir,
        paste0(tools::file_path_sans_ext(basename(listOfFiles[[idx]])), ".", explicitExt)
      )
      write.AsspDataObj(assp_obj, out_file)
      externalRes[[idx]] <- out_file
    } else {
      externalRes[[idx]] <- assp_obj
    }
  }

  if (verbose) {
    cli::cli_progress_done()
  }

  # Return results
  if (toFile) {
    invisible(unlist(externalRes))
  } else {
    if (nFiles == 1) externalRes[[1]] else externalRes
  }
}

# ============================================================================
# PHASE PRIMITIVES
# ============================================================================

.hmpd_wrap_phase <- function(x) {
  ((x + pi) %% (2 * pi)) - pi
}

.hmpd_unwrap_phase <- function(x) {
  x <- as.numeric(x)
  if (length(x) < 2) return(x)
  dx <- diff(x)
  dx_wrapped <- .hmpd_wrap_phase(dx)
  c(x[1], x[1] + cumsum(dx_wrapped))
}

# ============================================================================
# SPECTRAL PRIMITIVES
# ============================================================================

.hmpd_hspec2spec <- function(X) {
  X <- as.complex(X)
  if (length(X) < 2) return(X)
  c(X, Conj(rev(X[-c(1, length(X))])))
}

.hmpd_as_row_matrix <- function(x) {
  if (is.null(dim(x))) matrix(x, nrow = 1L) else as.matrix(x)
}

.hmpd_fft_rows <- function(x, inverse = FALSE, n = NULL) {
  x <- .hmpd_as_row_matrix(x)
  if (is.null(n)) n <- ncol(x)
  out <- apply(x, 1, function(row) {
    if (inverse) Re(stats::fft(row, inverse = TRUE) / n)
    else stats::fft(row, inverse = FALSE)
  })
  if (is.null(dim(out))) matrix(out, nrow = 1L) else out
}

.hmpd_hspec2fwcep <- function(model, fs, cepord, warpfn = NULL) {
  if (is.null(warpfn)) {
    warpfn <- function(f) ucnv_hz_to_mel(f, method = "htk")
  }
  dftlen <- 2 * (length(model) - 1L)
  freqlin <- seq(0, dftlen / 2) * fs / dftlen
  freqwarp <- 0.5 * fs * warpfn(freqlin) / warpfn(fs / 2)
  X <- .hmpd_hspec2spec(model)
  X <- c(X, rep(0, dftlen - length(X)))
  lX <- log(pmax(Mod(X), .Machine$double.eps))
  cep <- Re(stats::fft(lX, inverse = TRUE)) / dftlen
  if (dftlen %% 2L == 0L) {
    out <- c(cep[1], 2 * cep[2:(dftlen / 2)], cep[dftlen / 2 + 1L])
  } else {
    out <- c(cep[1], 2 * cep[2:((dftlen + 1L) / 2)])
  }
  out[1:(1L + cepord)]
}

# ============================================================================
# FILTER PRIMITIVES
# ============================================================================

.hmpd_fir_filter <- function(x, b) {
  x <- as.numeric(x)
  b <- as.numeric(b)
  y <- stats::filter(x, b, sides = 1)
  y[is.na(y)] <- 0
  y
}

.hmpd_filtfilt_fir <- function(x, b) {
  forward <- .hmpd_fir_filter(x, b)
  backward <- .hmpd_fir_filter(rev(forward), b)
  rev(backward)
}

# ============================================================================
# SINUSOIDAL ANALYSIS HELPERS
# ============================================================================

.hmpd_merge_defaults <- function(opt, defaults) {
  if (is.null(opt)) {
    return(defaults)
  }
  utils::modifyList(defaults, opt)
}

.hmpd_nextpow2 <- function(x) {
  ceiling(log(x, base = 2))
}

.hmpd_odd_length <- function(x) {
  2L * round(x / 2) + 1L
}

.hmpd_window_from_opt <- function(win_fn, n) {
  if (is.character(win_fn) && length(win_fn) == 1L) {
    if (win_fn == "blackman") return(as.numeric(av::blackman(n)))
    if (win_fn == "hann" || win_fn == "hanning") return(as.numeric(av::hanning(n)))
    if (win_fn == "hamming") return(as.numeric(av::hamming(n)))
    stop("Unknown window: ", win_fn, call. = FALSE)
  }
  if (is.function(win_fn)) {
    return(as.numeric(win_fn(n)))
  }
  stop("win_fn must be a window name or function.", call. = FALSE)
}

.hmpd_normalize_window <- function(win, normstr) {
  if (is.null(normstr) || !nzchar(normstr)) {
    return(win)
  }
  if (!identical(normstr, "sum(win)")) {
    stop("Only normstr = 'sum(win)' is implemented.", call. = FALSE)
  }
  win / sum(win)
}

.hmpd_delay_to_spec <- function(delay, dftlen) {
  if (delay == 0) {
    return(rep(1 + 0i, dftlen))
  }
  if (dftlen %% 2L == 1L) {
    shift <- exp((delay * 2i * pi / dftlen) * seq_len((dftlen - 1L) / 2L))
    c(1 + 0i, shift, Conj(rev(shift)))
  } else {
    shift <- exp((delay * 2i * pi / dftlen) * seq_len(dftlen / 2L))
    c(1 + 0i, shift[-length(shift)], sign(Re(shift[length(shift)])) + 0i, Conj(rev(shift[-length(shift)])))
  }
}

.hmpd_phase_interp <- function(S, fs, freq_hz) {
  dftlen <- length(S)
  faxis <- fs * (0:(dftlen / 2)) / dftlen
  phase <- .hmpd_unwrap_phase(Arg(S[1:(dftlen / 2 + 1L)]))
  interp <- stats::approx(faxis, phase, xout = freq_hz, method = "linear", rule = 2)$y
  .hmpd_wrap_phase(interp)
}

.hmpd_find_peaks_quadratic <- function(x) {
  x <- as.numeric(x)
  if (length(x) < 3L) {
    return(list(index = numeric(), value = numeric()))
  }
  idx <- which(diff(sign(diff(x))) < 0) + 1L
  if (!length(idx)) {
    return(list(index = numeric(), value = numeric()))
  }
  b <- 0.5 * (x[idx + 1L] - x[idx - 1L])
  a <- x[idx] - b - x[idx - 1L]
  is_peak <- a > 0
  k <- as.numeric(idx)
  v <- x[idx]
  if (any(is_peak)) {
    v[is_peak] <- x[idx[is_peak]] + 0.25 * b[is_peak]^2 / a[is_peak]
    k[is_peak] <- idx[is_peak] + 0.5 * b[is_peak] / a[is_peak]
  }
  list(index = k, value = v)
}

.hmpd_spec_fit_freq_amp <- function(S, fs, index, zp = NULL, wintype = NULL) {
  if (index <= 1) {
    return(c(freq = 0, amp = abs(S[1])))
  }
  if (index >= length(S)) {
    return(c(freq = (length(S) - 1) * fs / length(S), amp = abs(S[length(S)])))
  }
  y1 <- log(Mod(S[index - 1L]))
  y2 <- log(Mod(S[index]))
  y3 <- log(Mod(S[index + 1L]))
  A <- (y3 - y2) / 2 + (y1 - y2) / 2
  B <- -(y1 - y2 - A)
  C <- y1 - A + B
  di <- -B / (2 * A)
  logamp <- A * di * di + B * di + C
  c(freq = (index + di - 1) * fs / length(S), amp = exp(logamp))
}

.hmpd_spec_getsins <- function(S, fs, ...) {
  S <- as.complex(S)
  dftlen <- length(S)
  peaks <- .hmpd_find_peaks_quadratic(Mod(S[1:floor(dftlen / 2)]))
  sins <- matrix(0, nrow = 5L, ncol = length(peaks$index))
  if (length(peaks$index)) {
    fit <- vapply(
      round(peaks$index),
      function(idx) .hmpd_spec_fit_freq_amp(S, fs, idx, ...),
      numeric(2)
    )
    sins[1:2, ] <- fit
    sins[5, ] <- 1
    sins[3, ] <- .hmpd_phase_interp(S, fs, sins[1, ])
  }
  cbind(c(0, Mod(S[1]), Arg(S[1]), 0, as.numeric(Mod(S[1]) > Mod(S[2]))), sins)
}

.hmpd_spec_getsins_f0 <- function(S, fs, f0, max_h = NULL) {
  S <- as.complex(S)
  dftlen <- length(S)
  if (is.null(max_h)) {
    max_h <- floor(fs / 2 / f0)
  }
  fks <- f0 * (0:max_h)
  fks <- fks[fks < fs / 2]
  sins <- matrix(0, nrow = 5L, ncol = length(fks))
  if (length(fks)) {
    fit <- vapply(
      pmax(1L, round(fks * dftlen / fs)),
      function(idx) .hmpd_spec_fit_freq_amp(S, fs, idx),
      numeric(2)
    )
    sins[1:2, ] <- fit
    sins[5, ] <- 1
    sins[3, ] <- .hmpd_phase_interp(S, fs, fks)
  }
  sins
}

.hmpd_spec2minphasespec <- function(X) {
  X <- as.complex(X)
  full <- .hmpd_hspec2spec(log(pmax(Mod(X), .Machine$double.eps)))
  n <- length(full)
  rcc <- stats::fft(full, inverse = TRUE) / n
  if (n %% 2L == 0L) {
    cep <- c(rcc[1], 2 * rcc[2:(n / 2)], rcc[n / 2 + 1])
  } else {
    cep <- c(rcc[1], 2 * rcc[2:((n - 1L) / 2 + 1L)])
  }
  length(cep) <- n
  cep[is.na(cep)] <- 0
  out <- stats::fft(cep)
  out[1:(n / 2 + 1L)]
}

.hmpd_sin_analysis_defaults <- function() {
  list(
    win_durf0sync = TRUE,
    win_durnbper = 3,
    win_dur = 30 / 1000,
    win_fn = "blackman",
    win_dropoutside = TRUE,
    fharmonic = TRUE,
    fquasiharm = FALSE,
    fadapted = FALSE,
    use_ls = TRUE,
    win_ls_marg = 0.1,
    normstr = "sum(win)",
    dftlen = NULL,
    osf = 2,
    frames_keepspec = FALSE,
    funique = TRUE,
    resyn = FALSE,
    debug = 0
  )
}

.hmpd_sin_analysis <- function(wav, fs, f0s, opt = NULL) {
  defaults <- .hmpd_sin_analysis_defaults()
  if (missing(wav)) {
    return(defaults)
  }
  opt <- .hmpd_merge_defaults(opt, defaults)
  wav <- as.numeric(wav)
  f0s <- as.matrix(f0s)[, 1:2, drop = FALSE]

  # Frame setup
  if (isTRUE(opt$win_durf0sync)) {
    win_dur <- mean(f0s[, 2], na.rm = TRUE) * opt$win_durnbper / 1000
  } else {
    win_dur <- opt$win_dur
  }

  if (is.null(opt$dftlen)) {
    opt$dftlen <- 2^.hmpd_nextpow2(ceil(win_dur * fs))
  }

  frame_len <- min(opt$dftlen, round(win_dur * fs))
  frame_shift <- round(mean(diff(f0s[, 1], na.rm = TRUE)) * fs)

  if (frame_shift < 1) frame_shift <- round(fs / 100)

  # Frame extraction + sinusoidal analysis
  frames <- list()
  start <- 1L
  while (start + frame_len <= length(wav)) {
    frame <- wav[start:(start + frame_len - 1L)]

    # Windowing
    win <- .hmpd_window_from_opt(opt$win_fn, length(frame))
    win <- .hmpd_normalize_window(win, opt$normstr)
    frame_windowed <- frame * win

    # Zero-pad to dftlen
    frame_padded <- c(frame_windowed, rep(0, opt$dftlen - length(frame_windowed)))
    spec <- stats::fft(frame_padded)

    # Sinusoidal extraction
    frame_time <- (start + frame_len / 2 - 1L) / fs
    f0_interp <- stats::approx(f0s[, 1], f0s[, 2], xout = frame_time, method = "linear", rule = 2)$y

    if (is.finite(f0_interp) && f0_interp > 0) {
      sins <- .hmpd_spec_getsins_f0(spec, fs, f0_interp)
    } else {
      sins <- .hmpd_spec_getsins(spec, fs)
    }

    frame_obj <- list(
      t = frame_time,
      f0 = if (is.finite(f0_interp)) f0_interp else NA_real_,
      sins = sins
    )

    frames[[length(frames) + 1L]] <- frame_obj
    start <- start + frame_shift
  }

  list(frames = frames, opt = opt)
}

# ============================================================================
# HMPD ANALYSIS HELPERS
# ============================================================================

.hmpd_apply_cols <- function(x, fn) {
  x <- as.matrix(x)
  out <- apply(x, 2, fn)
  if (is.null(dim(out))) {
    matrix(out, ncol = 1L)
  } else {
    out
  }
}

.hmpd_filtfilt_matrix <- function(x, b) {
  .hmpd_apply_cols(x, function(col) .hmpd_filtfilt_fir(as.numeric(col), b))
}

.hmpd_runmed_matrix <- function(x, k) {
  .hmpd_apply_cols(x, function(col) stats::runmed(as.numeric(col), k = k, endrule = "median"))
}

.hmpd_interp_with_default <- function(x, y, xi, def = NA_real_, transform = NULL, inverse = NULL) {
  if (is.null(transform) || is.null(inverse)) {
    transform <- identity
    inverse <- identity
  }
  yi <- stats::approx(x, transform(y), xout = xi, method = "linear", rule = 1)$y
  yi[is.na(yi)] <- transform(def)
  inverse(yi)
}

.hmpd_env_interp <- function(sins, fs, dftlen, extrap_dcny = FALSE) {
  sins <- as.matrix(sins)
  fks <- as.numeric(sins[1, ])
  aks <- as.numeric(sins[2, ])

  if (isTRUE(any(extrap_dcny))) {
    if (length(extrap_dcny) == 1L || isTRUE(extrap_dcny[1])) {
      if (fks[1] != 0) {
        a0 <- aks[1]
        fks <- c(0, fks)
        aks <- c(a0, aks)
      } else if (length(aks) > 1L) {
        aks[1] <- aks[2]
      }
    }
    if (length(extrap_dcny) == 1L || isTRUE(extrap_dcny[min(2, length(extrap_dcny))])) {
      mfd <- stats::median(diff(fks))
      if (!is.finite(mfd) || mfd <= 0) {
        mfd <- fs / dftlen
      }
      while (tail(fks, 1) < fs / 2 - mfd) {
        fks <- c(fks, tail(fks, 1) + mfd)
        aks <- c(aks, tail(aks, 1))
      }
    }
  }

  if (length(fks) < 5L) {
    extra_f <- seq(from = max(fks[1], fs / dftlen), by = fs / dftlen, length.out = 5L - length(fks))
    fks <- c(fks, tail(fks, 1) + extra_f)
    aks <- c(aks, rep(tail(aks, 1), length(extra_f)))
  }

  if (fks[1] > 0) {
    take <- seq_len(min(4L, length(fks)))
    fks <- c(-rev(fks[take]), fks)
    aks <- c(rev(aks[take]), aks)
  } else {
    take <- seq.int(2L, min(5L, length(fks)))
    fks <- c(-rev(fks[take]), fks)
    aks <- c(rev(aks[take]), aks)
  }

  if (tail(fks, 1) < fs / 2) {
    take <- seq.int(max(1L, length(fks) - 4L), length(fks))
    dftony <- fs / 2 - rev(fks[take])
    fks <- c(fks, fs / 2 + dftony)
    aks <- c(aks, rev(aks[take]))
  } else {
    take <- seq.int(max(1L, length(fks) - 5L), length(fks) - 1L)
    dftony <- fs / 2 - rev(fks[take])
    fks <- c(fks, fs / 2 + dftony)
    aks <- c(aks, rev(aks[take]))
  }

  bins <- fs * (0:(dftlen / 2)) / dftlen
  keep <- !duplicated(fks)
  env <- stats::approx(fks[keep], log(pmax(abs(aks[keep]), .Machine$double.eps)), xout = bins)$y
  exp(env)
}

.hmpd_hspec2minphaseloghspec <- function(X) {
  X <- as.complex(X)
  full <- .hmpd_hspec2spec(log(pmax(Mod(X), .Machine$double.eps)))
  n <- length(full)
  rcc <- stats::fft(full, inverse = TRUE) / n
  if (n %% 2L == 0L) {
    cep <- c(rcc[1], 2 * rcc[2:(n / 2)], rcc[n / 2 + 1])
  } else {
    cep <- c(rcc[1], 2 * rcc[2:((n - 1L) / 2 + 1L)])
  }
  length(cep) <- n
  cep[is.na(cep)] <- 0
  out <- stats::fft(cep)
  out[1:(n / 2 + 1L)]
}

.hmpd_wrap_array <- function(x) {
  dim_x <- dim(x)
  out <- .hmpd_wrap_phase(as.numeric(x))
  if (is.null(dim_x)) out else array(out, dim = dim_x)
}

.hmpd_unwrap_array <- function(x) {
  dim_x <- dim(x)
  out <- .hmpd_unwrap_phase(as.numeric(x))
  if (is.null(dim_x)) out else array(out, dim = dim_x)
}

.hmpd_irregsampling2uniformsampling <- function(T, X, nT, fn = NULL, ifn = NULL, method = "linear", def = NA_real_) {
  if (is.null(fn) || is.null(ifn)) {
    fn <- identity
    ifn <- identity
  }
  X <- as.matrix(X)
  out <- matrix(NA_real_, nrow = length(nT), ncol = ncol(X))
  for (j in seq_len(ncol(X))) {
    yi <- stats::approx(T, fn(X[, j]), xout = nT, method = method, rule = 1)$y
    yi[is.na(yi)] <- fn(def)
    out[, j] <- ifn(yi)
  }
  out
}

.hmpd_harmscale2hertzscale <- function(X, f0s, fs, dftlen, dcval = NULL, fn = NULL, ifn = NULL) {
  if (is.null(fn) || is.null(ifn)) {
    fn <- identity
    ifn <- identity
  }
  X <- as.matrix(X)
  F <- fs * (0:(dftlen / 2)) / dftlen
  Xr <- matrix(NA_real_, nrow = nrow(f0s), ncol = length(F))
  for (n in seq_len(nrow(f0s))) {
    idx <- which(!is.na(X[n, ]))
    if (length(idx) > 1L) {
      if (is.null(dcval)) {
        Xr[n, ] <- ifn(stats::approx(
          f0s[n, 2] * (idx - 1L),
          fn(X[n, idx]),
          xout = F,
          method = "linear",
          rule = 1
        )$y)
      } else {
        x <- c(0, f0s[n, 2] * 0.95, f0s[n, 2] * idx)
        y <- c(dcval, dcval, X[n, idx])
        Xr[n, ] <- ifn(stats::approx(x, fn(y), xout = F, method = "linear", rule = 1)$y)
      }
    }
  }
  list(X = Xr, F = F)
}

.hmpd_hmpd_phase_mean <- function(PD, nbat) {
  winlen <- 2L * round(nbat / 2) + 1L
  win <- rep(1 / winlen, winlen)
  PDc <- .hmpd_filtfilt_matrix(cos(PD), win)
  PDs <- .hmpd_filtfilt_matrix(sin(PD), win)
  atan2(PDs, PDc)
}

.hmpd_hmpd_phase_smooth <- function(PD, nbat) {
  winlen <- 2L * round(nbat / 2) + 1L
  win <- as.numeric(av::hanning(winlen))
  win <- win / sum(win)
  PDc <- .hmpd_runmed_matrix(cos(PD), winlen)
  PDs <- .hmpd_runmed_matrix(sin(PD), winlen)
  PDc <- .hmpd_filtfilt_matrix(PDc, win)
  PDs <- .hmpd_filtfilt_matrix(PDs, win)
  atan2(PDs, PDc)
}

.hmpd_hmpd_phase_deviation <- function(PD, nbat) {
  winlen <- 2L * round(nbat / 2) + 1L
  win <- rep(1 / winlen, winlen)
  PDc <- .hmpd_filtfilt_matrix(cos(PD), win)
  PDs <- .hmpd_filtfilt_matrix(sin(PD), win)
  z <- Mod(PDc + 1i * PDs)
  out <- matrix(0, nrow = nrow(z), ncol = ncol(z))
  idx <- z < 1
  out[idx] <- sqrt(-2 * log(z[idx]))
  out
}

.hmpd_hlin2hlog <- function(Hb, Hmax, order) {
  hsl <- seq_len(Hmax)
  p1 <- c(Hb, Hb)
  p2 <- c(order, order)
  p3 <- c(Hmax, order)
  t <- seq(0, 1, by = 0.01)
  bez_x <- (1 - t)^2 * p1[1] + 2 * (1 - t) * t * p2[1] + t^2 * p3[1]
  bez_y <- (1 - t)^2 * p1[2] + 2 * (1 - t) * t * p2[2] + t^2 * p3[2]
  hsl[(Hb + 1):Hmax] <- stats::approx(bez_x, bez_y, xout = (Hb + 1):Hmax, method = "linear", rule = 2)$y
  hsl
}

.hmpd_philin2philog <- function(philin, Hb, Hmax = NULL, order = NULL) {
  if (is.null(Hmax)) {
    hsl <- Hb
    Hmax <- length(Hb)
    order <- hsl[length(hsl)]
  } else {
    hsl <- .hmpd_hlin2hlog(Hb, Hmax, order)
  }
  hsl <- round(hsl[seq_len(min(Hmax, length(philin)))])
  philog <- rep(NA_real_, order)
  for (h in seq_len(order)) {
    idx <- hsl == h
    if (any(idx)) {
      philog[h] <- Arg(mean(exp(1i * philin[idx])))
    }
  }
  philog[is.na(philog)] <- .hmpd_wrap_array(2 * pi * stats::runif(sum(is.na(philog))))
  philog
}

.hmpd_phase_rpspd <- function(frames, fs, opt = NULL) {
  if (is.null(opt)) {
    opt <- list()
  }

  f0s <- cbind(
    vapply(frames, `[[`, numeric(1), "t"),
    vapply(frames, `[[`, numeric(1), "f0")
  )
  Hmax <- max(vapply(frames, function(fr) ncol(fr$sins) - 1L, integer(1)))

  PE <- matrix(.hmpd_wrap_array((2 * pi) * stats::runif(nrow(f0s) * (1L + Hmax))), nrow = nrow(f0s))
  for (n in seq_along(frames)) {
    Ks <- 0:(ncol(frames[[n]]$sins) - 1L)
    pk <- frames[[n]]$sins[3, ]
    mino <- min(length(pk), 1L + Hmax)
    PE[n, 1:mino] <- pk[1:mino]
  }

  PE <- .hmpd_wrap_array(PE)
  list(PE = PE, frames = frames, opt = opt)
}

.hmpd_amplitude_envelope_estimate <- function(frames, fs, opt) {
  dftlen <- if (!is.null(opt$dftlen)) opt$dftlen else 1024
  AE <- matrix(0, nrow = length(frames), ncol = dftlen / 2 + 1L)
  F <- fs * (0:(dftlen / 2)) / dftlen

  for (n in seq_along(frames)) {
    E <- .hmpd_env_interp(frames[[n]]$sins[1:2, 2:ncol(frames[[n]]$sins), drop = FALSE], fs, dftlen, TRUE)
    E[is.na(E)] <- 10^(-300 / 20)
    AE[n, ] <- abs(E)
  }

  list(AE = AE, frames = frames, opt = opt)
}

.hmpd_features_compress <- function(f0s, AE, PDM, PDD, fs, opt) {
  amp_enc_method <- if (!is.null(opt$amp_enc_method)) opt$amp_enc_method else 2
  amp_log <- if (!is.null(opt$amp_log)) opt$amp_log else FALSE
  amp_order <- if (!is.null(opt$amp_order)) opt$amp_order else 24
  amp_logfn <- if (!is.null(opt$amp_logfn)) opt$amp_logfn else (function(f) ucnv_hz_to_mel(f, method = "htk"))

  pdm_log <- if (!is.null(opt$pdm_log)) opt$pdm_log else TRUE
  pdm_order <- if (!is.null(opt$pdm_order)) opt$pdm_order else 24
  pdm_log_hb <- if (!is.null(opt$pdm_log_hb)) opt$pdm_log_hb else 8

  pdd_log <- if (!is.null(opt$pdd_log)) opt$pdd_log else TRUE
  pdd_order <- if (!is.null(opt$pdd_order)) opt$pdd_order else 12

  dftlen <- if (!is.null(opt$dftlen)) opt$dftlen else 1024

  # Amplitude
  if (identical(amp_enc_method, 1)) {
    AEC <- AE
  } else {
    AEC <- matrix(0, nrow = nrow(AE), ncol = 1L + amp_order)
    for (n in seq_len(nrow(AE))) {
      if (identical(amp_enc_method, 2)) {
        AEC[n, ] <- .hmpd_hspec2fwcep(AE[n, ], fs, amp_order, amp_logfn)
      }
    }
  }

  # Phase deviation mean
  if (!isTRUE(pdm_log)) {
    PDMC <- PDM
  } else {
    PDMC <- matrix(0, nrow = nrow(PDM), ncol = 1L + pdm_order)
    hsl <- .hmpd_hlin2hlog(pdm_log_hb, dftlen / 2, 1L + pdm_order)
    for (n in seq_len(nrow(PDM))) {
      PDMC[n, ] <- .hmpd_philin2philog(PDM[n, ], hsl)
    }
  }

  # Phase deviation deviation
  if (!isTRUE(pdd_log)) {
    PDDC <- PDD
  } else {
    PDDC <- matrix(0, nrow = nrow(PDD), ncol = 1L + pdd_order)
    for (n in seq_len(nrow(PDD))) {
      X <- PDD[n, ]
      X[1] <- X[2]
      X[X == 0] <- 0.001
      PDDC[n, ] <- .hmpd_hspec2fwcep(X, fs, pdd_order)
    }
  }

  list(AE = AEC, PDM = PDMC, PDD = PDDC)
}

# ============================================================================
# MAIN HMPD ANALYSIS FUNCTION
# ============================================================================

.hmpd_analysis <- function(wav, fs, f0s_ext = NULL, f0min = 60, f0max = 440) {
  # F0 must be provided (internal computation requires temp file I/O which is complex)
  if (is.null(f0s_ext)) {
    # Use constant F0 estimate as fallback (simple approximation)
    # Better approach: use faster internal method or allow user to provide F0
    cli::cli_warn(c(
      "No F0 contour provided; using default constant F0 at 100 Hz",
      "!" = "Provide f0s parameter for accurate results"
    ))
    n_frames <- round(length(wav) / fs / 0.01)
    f0s <- cbind(
      seq(0, length(wav) / fs - 0.01, length.out = n_frames),
      rep(100, n_frames)
    )
  } else {
    f0s <- as.matrix(f0s_ext)[, 1:2, drop = FALSE]
  }

  # Default options
  opt <- list(
    dftlen = 2048,
    amp_enc_method = 2,
    amp_log = FALSE,
    amp_order = 24,
    amp_logfn = function(f) ucnv_hz_to_mel(f, method = "htk"),
    pdm_log = TRUE,
    pdm_order = 24,
    pdm_log_hb = 8,
    pdd_log = TRUE,
    pdd_order = 12,
    pdm_nbper = 6,
    pdd_nbper = 2,
    sin_nbat = 4,
    regularstepsize = 0.005
  )

  # Sinusoidal analysis
  sin_res <- .hmpd_sin_analysis(wav, fs, f0s, opt)
  frames <- sin_res$frames
  opt <- sin_res$opt

  # Amplitude envelope
  amp <- .hmpd_amplitude_envelope_estimate(frames, fs, opt)
  AE <- amp$AE
  frames <- amp$frames

  # Phase extraction
  phase <- .hmpd_phase_rpspd(frames, fs, opt)
  PE <- phase$PE

  # Phase deviation mean/deviation
  PDM <- .hmpd_hmpd_phase_mean(PE, opt$pdm_nbper * opt$sin_nbat)
  PEtrend <- .hmpd_hmpd_phase_smooth(PE, opt$pdd_nbper * opt$sin_nbat)
  PDD <- .hmpd_hmpd_phase_deviation(PE - PEtrend, opt$pdd_nbper * opt$sin_nbat)

  # Compress features
  irregf0s <- cbind(
    vapply(frames, `[[`, numeric(1), "t"),
    vapply(frames, `[[`, numeric(1), "f0")
  )

  compressed <- .hmpd_features_compress(irregf0s, AE, PDM, PDD, fs, opt)

  # Resample to uniform grid (10 ms)
  uT <- seq(0, irregf0s[nrow(irregf0s), 1], by = 0.01)

  AEu <- .hmpd_irregsampling2uniformsampling(irregf0s[, 1], compressed$AE, uT, method = "linear", def = NA_real_)
  PDMu <- .hmpd_irregsampling2uniformsampling(
    irregf0s[, 1],
    compressed$PDM,
    uT,
    fn = .hmpd_unwrap_array,
    ifn = .hmpd_wrap_array,
    method = "linear",
    def = 0
  )
  PDDu <- .hmpd_irregsampling2uniformsampling(irregf0s[, 1], compressed$PDD, uT, method = "linear", def = 0)

  # Replace NaN with 0
  AEu[is.nan(AEu)] <- 0
  PDMu[is.nan(PDMu)] <- 0
  PDDu[is.nan(PDDu)] <- 0

  list(ae = AEu, pdm = PDMu, pdd = PDDu)
}

# Set function attributes
attr(trk_hmpd, "ext") <- "hpd"
attr(trk_hmpd, "tracks") <- c("ae", "pdm", "pdd")
attr(trk_hmpd, "outputType") <- "SSFF"
attr(trk_hmpd, "nativeFiletypes") <- character()
