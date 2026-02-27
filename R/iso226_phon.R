#' ISO 226:2023 Phon (Loudness Level) Conversions
#'
#' @description
#' Functions for converting between sound pressure level (dB) and loudness level (phon)
#' according to ISO 226:2023 "Acoustics - Normal equal-loudness-level contours".
#'
#' @details
#' The phon scale represents loudness level - the perceived loudness of a sound
#' compared to a 1000 Hz reference tone. Unlike other psychoacoustic scales
#' (Bark, ERB, Mel), phon is a **bivariate measure** requiring both frequency
#' and sound pressure level.
#'
#' **Key concept**: A sound at 40 phon sounds as loud as a 1000 Hz tone at 40 dB SPL,
#' but different frequencies require different SPLs to achieve the same loudness.
#'
#' **Valid ranges** (per ISO 226:2023):
#' - Frequency: 20 Hz to 12,500 Hz
#' - Loudness level: 20 phon to 90 phon (20-4000 Hz), 20-80 phon (5000-12500 Hz)
#' - Below 20 phon: Informative only (near hearing threshold)
#' - Above 90/80 phon: Limited experimental data
#'
#' **Implementation notes**:
#' - Based on ISO 226:2023 Formulas (1) and (2)
#' - Parameters from ISO 226:2023 Table 1
#' - Interpolation used for frequencies between standard 1/3-octave values
#' - Free-field listening conditions (frontal incidence, binaural)
#'
#' @references
#' ISO 226:2023, Acoustics - Normal equal-loudness-level contours.
#' International Organization for Standardization, Geneva, Switzerland.
#'
#' @name iso226_phon
NULL

#' ISO 226:2023 Parameter Table
#'
#' Parameters for calculating equal-loudness-level contours according to
#' ISO 226:2023 Table 1. These are used in the formulas for converting between
#' sound pressure level and loudness level.
#'
#' @format A data frame with 29 rows (one per 1/3-octave frequency) and 4 columns:
#' \describe{
#'   \item{freq_hz}{Frequency in Hz (20 to 12500 Hz)}
#'   \item{alpha_f}{Exponent for loudness perception at frequency f}
#'   \item{L_U}{Magnitude of linear transfer function normalized at 1000 Hz (dB)}
#'   \item{T_f}{Threshold of hearing at frequency f (dB)}
#' }
#'
#' @source ISO 226:2023, Table 1
#' @keywords internal
.iso226_params <- data.frame(
  freq_hz = c(20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315,
              400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000,
              5000, 6300, 8000, 10000, 12500),
  alpha_f = c(0.635, 0.602, 0.569, 0.537, 0.509, 0.482, 0.456, 0.433, 0.412,
              0.391, 0.373, 0.357, 0.343, 0.330, 0.320, 0.311, 0.303, 0.300,
              0.295, 0.292, 0.290, 0.290, 0.289, 0.289, 0.289, 0.293, 0.303,
              0.323, 0.354),
  L_U = c(-31.5, -27.2, -23.1, -19.3, -16.1, -13.1, -10.4, -8.2, -6.3, -4.6,
          -3.2, -2.1, -1.2, -0.5, 0.0, 0.4, 0.5, 0.0, -2.7, -4.2, -1.2, 1.4,
          2.3, 1.0, -2.3, -7.2, -11.2, -10.9, -3.5),
  T_f = c(78.1, 68.7, 59.5, 51.1, 44.0, 37.5, 31.5, 26.5, 22.1, 17.9, 14.4,
          11.4, 8.6, 6.2, 4.4, 3.0, 2.2, 2.4, 3.5, 1.7, -1.3, -4.2, -6.0,
          -5.4, -1.5, 6.0, 12.6, 13.9, 12.3),
  stringsAsFactors = FALSE
)

#' Get ISO 226 Parameters for a Frequency
#'
#' Retrieves or interpolates the ISO 226:2023 parameters (alpha_f, L_U, T_f)
#' for a given frequency. Uses linear interpolation (on log-frequency scale)
#' for frequencies between the standard 1/3-octave values.
#'
#' @param freq_hz Numeric; frequency in Hz (20 to 12500 Hz)
#' @return Named list with elements: alpha_f, L_U, T_f, T_r (reference threshold)
#' @keywords internal
.get_iso226_params <- function(freq_hz) {
  # Validate frequency range
  if (any(freq_hz < 20 | freq_hz > 12500)) {
    stop("Frequency must be between 20 Hz and 12500 Hz (ISO 226:2023 valid range)",
         call. = FALSE)
  }

  # Reference values at 1000 Hz (always the same)
  alpha_r <- 0.300
  T_r <- 2.4

  # Vectorized lookup/interpolation
  n <- length(freq_hz)
  result <- list(
    alpha_f = numeric(n),
    L_U = numeric(n),
    T_f = numeric(n),
    alpha_r = alpha_r,
    T_r = T_r
  )

  for (i in seq_along(freq_hz)) {
    f <- freq_hz[i]

    # Check if exact match in table
    exact_idx <- which(.iso226_params$freq_hz == f)

    if (length(exact_idx) > 0) {
      # Exact match - use table values directly
      result$alpha_f[i] <- .iso226_params$alpha_f[exact_idx]
      result$L_U[i] <- .iso226_params$L_U[exact_idx]
      result$T_f[i] <- .iso226_params$T_f[exact_idx]
    } else {
      # Interpolate on log-frequency scale
      log_f <- log10(f)
      log_freq_table <- log10(.iso226_params$freq_hz)

      # Find surrounding frequencies
      idx_upper <- which(log_freq_table > log_f)[1]
      idx_lower <- idx_upper - 1

      if (is.na(idx_upper)) {
        # Beyond upper limit - use last value (with warning)
        warning("Frequency ", f, " Hz is above highest standard frequency. ",
                "Using parameters from 12500 Hz.", call. = FALSE)
        result$alpha_f[i] <- .iso226_params$alpha_f[nrow(.iso226_params)]
        result$L_U[i] <- .iso226_params$L_U[nrow(.iso226_params)]
        result$T_f[i] <- .iso226_params$T_f[nrow(.iso226_params)]
      } else {
        # Linear interpolation on log-frequency scale
        log_f1 <- log_freq_table[idx_lower]
        log_f2 <- log_freq_table[idx_upper]
        weight <- (log_f - log_f1) / (log_f2 - log_f1)

        result$alpha_f[i] <- .iso226_params$alpha_f[idx_lower] +
          weight * (.iso226_params$alpha_f[idx_upper] - .iso226_params$alpha_f[idx_lower])

        result$L_U[i] <- .iso226_params$L_U[idx_lower] +
          weight * (.iso226_params$L_U[idx_upper] - .iso226_params$L_U[idx_lower])

        result$T_f[i] <- .iso226_params$T_f[idx_lower] +
          weight * (.iso226_params$T_f[idx_upper] - .iso226_params$T_f[idx_lower])
      }
    }
  }

  return(result)
}

#' Convert Sound Pressure Level and Frequency to Loudness Level (Phon)
#'
#' @description
#' Converts sound pressure level (dB) at a given frequency to loudness level (phon)
#' according to ISO 226:2023 equal-loudness-level contours.
#'
#' @param spl_db Numeric vector; sound pressure level in dB (re 20 μPa)
#' @param freq_hz Numeric vector; frequency in Hz (20 to 12500 Hz)
#'
#' @return Numeric vector; loudness level in phon
#'
#' @details
#' This function implements Formula (2) from ISO 226:2023 Section 4.2.
#'
#' **Loudness level (phon)** represents the perceived loudness of a sound.
#' By definition, a sound at N phon has the same loudness as a 1000 Hz tone
#' at N dB SPL.
#'
#' **Valid ranges**:
#' - Frequency: 20-12500 Hz (1/3-octave standard frequencies)
#' - SPL: Threshold of hearing to ~130 dB (pain threshold)
#' - Output: 20-90 phon (reliable data), <20 phon informative only
#'
#' **Interpolation**: Frequencies between standard 1/3-octave values are
#' interpolated on a log-frequency scale.
#'
#' **Vectorization**: Both spl_db and freq_hz can be vectors. If both are vectors,
#' they must be the same length, or one must be length 1 (recycled).
#'
#' @examples
#' \dontrun{
#' # 40 dB at 1000 Hz = 40 phon (by definition)
#' ucnv_db_and_hz_to_phon(40, 1000)
#'
#' # Same loudness at different frequencies
#' ucnv_db_and_hz_to_phon(60, 100)   # Low frequency needs more SPL
#' ucnv_db_and_hz_to_phon(40, 1000)  # Reference
#' ucnv_db_and_hz_to_phon(35, 4000)  # High frequency needs less SPL
#'
#' # Vectorized: multiple measurements
#' spl <- c(40, 50, 60, 70)
#' freq <- c(100, 500, 1000, 4000)
#' ucnv_db_and_hz_to_phon(spl, freq)
#'
#' # Single frequency, multiple SPLs
#' ucnv_db_and_hz_to_phon(c(30, 40, 50, 60), 1000)
#' }
#'
#' @references
#' ISO 226:2023, Acoustics - Normal equal-loudness-level contours.
#' Formula (2), Section 4.2, page 3.
#'
#' @seealso \code{\link{ucnv_phon_and_hz_to_db}} for the inverse conversion
#'
#' @export
ucnv_db_and_hz_to_phon <- function(spl_db, freq_hz) {
  # Input validation
  if (!is.numeric(spl_db) || !is.numeric(freq_hz)) {
    stop("Both spl_db and freq_hz must be numeric", call. = FALSE)
  }

  # Check vector lengths
  n_spl <- length(spl_db)
  n_freq <- length(freq_hz)

  if (n_spl != n_freq && n_spl != 1 && n_freq != 1) {
    stop("spl_db and freq_hz must have the same length, or one must be length 1",
         call. = FALSE)
  }

  # Recycle to common length
  n <- max(n_spl, n_freq)
  if (n_spl == 1) spl_db <- rep(spl_db, n)
  if (n_freq == 1) freq_hz <- rep(freq_hz, n)

  # Get parameters for each frequency
  params <- .get_iso226_params(freq_hz)

  # Extract parameters
  alpha_f <- params$alpha_f
  L_U <- params$L_U
  T_f <- params$T_f
  alpha_r <- params$alpha_r
  T_r <- params$T_r

  # Implement ISO 226:2023 Formula (2)
  # L_N = (100/3) * lg{ [10^(α_f*(L_f+L_U)/10dB) - 10^(α_f*T_f/10dB)] /
  #                     [(4*10^-10)^(0.3-α_f)] + 10^0.072 } phon

  # Calculate numerator terms
  term1 <- 10^(alpha_f * (spl_db + L_U) / 10)
  term2 <- 10^(alpha_f * T_f / 10)
  numerator <- term1 - term2

  # Calculate denominator
  base_term <- (4 * 10^(-10))^(0.3 - alpha_f)

  # Combined expression
  combined <- (numerator / base_term) + 10^0.072

  # Calculate loudness level in phon
  L_N <- (100 / 3) * log10(combined)

  # Warn if outside reliable range
  if (any(L_N < 20, na.rm = TRUE)) {
    n_low <- sum(L_N < 20, na.rm = TRUE)
    warning(n_low, " value(s) below 20 phon. Results are informative only ",
            "(near hearing threshold).", call. = FALSE)
  }

  # Check upper limits per ISO 226
  for (i in seq_along(L_N)) {
    f <- freq_hz[i]
    upper_limit <- if (f >= 5000) 80 else 90

    if (!is.na(L_N[i]) && L_N[i] > upper_limit) {
      warning("Frequency ", f, " Hz at ", round(L_N[i], 1),
              " phon exceeds reliable range (>", upper_limit,
              " phon). Limited experimental data.", call. = FALSE)
    }
  }

  return(L_N)
}

#' Convert Loudness Level (Phon) and Frequency to Sound Pressure Level
#'
#' @description
#' Converts loudness level (phon) at a given frequency to sound pressure level (dB)
#' according to ISO 226:2023 equal-loudness-level contours.
#'
#' @param phon Numeric vector; loudness level in phon (20 to 90/80 phon)
#' @param freq_hz Numeric vector; frequency in Hz (20 to 12500 Hz)
#'
#' @return Numeric vector; sound pressure level in dB (re 20 μPa)
#'
#' @details
#' This function implements Formula (1) from ISO 226:2023 Section 4.1.
#'
#' **Loudness level (phon)** represents the perceived loudness of a sound.
#' This function calculates the sound pressure level required at a given
#' frequency to achieve a specified loudness level.
#'
#' **Valid ranges**:
#' - Frequency: 20-12500 Hz (1/3-octave standard frequencies)
#' - Phon: 20-90 phon (20-4000 Hz), 20-80 phon (5000-12500 Hz)
#' - Below 20 phon: Informative only (near hearing threshold)
#' - Above limits: Limited experimental data
#'
#' **Interpolation**: Frequencies between standard 1/3-octave values are
#' interpolated on a log-frequency scale.
#'
#' **Vectorization**: Both phon and freq_hz can be vectors. If both are vectors,
#' they must be the same length, or one must be length 1 (recycled).
#'
#' @examples
#' \dontrun{
#' # 40 phon at 1000 Hz = 40 dB (by definition)
#' ucnv_phon_and_hz_to_db(40, 1000)
#'
#' # Same loudness (40 phon) at different frequencies
#' ucnv_phon_and_hz_to_db(40, 100)   # Low frequency needs higher SPL
#' ucnv_phon_and_hz_to_db(40, 1000)  # Reference: 40 dB
#' ucnv_phon_and_hz_to_db(40, 4000)  # High frequency needs lower SPL
#'
#' # Round-trip conversion check
#' spl_original <- 60
#' freq <- 1000
#' phon_calculated <- ucnv_db_and_hz_to_phon(spl_original, freq)
#' spl_recovered <- ucnv_phon_and_hz_to_db(phon_calculated, freq)
#' # spl_recovered ≈ spl_original
#'
#' # Calculate equal-loudness contour (40 phon)
#' frequencies <- c(100, 200, 500, 1000, 2000, 4000, 8000)
#' spls <- ucnv_phon_and_hz_to_db(40, frequencies)
#' plot(frequencies, spls, log = "x", type = "b",
#'      xlab = "Frequency (Hz)", ylab = "SPL (dB)",
#'      main = "40 Phon Equal-Loudness Contour")
#' }
#'
#' @references
#' ISO 226:2023, Acoustics - Normal equal-loudness-level contours.
#' Formula (1), Section 4.1, page 2.
#'
#' @seealso \code{\link{ucnv_db_and_hz_to_phon}} for the inverse conversion
#'
#' @export
ucnv_phon_and_hz_to_db <- function(phon, freq_hz) {
  # Input validation
  if (!is.numeric(phon) || !is.numeric(freq_hz)) {
    stop("Both phon and freq_hz must be numeric", call. = FALSE)
  }

  # Check vector lengths
  n_phon <- length(phon)
  n_freq <- length(freq_hz)

  if (n_phon != n_freq && n_phon != 1 && n_freq != 1) {
    stop("phon and freq_hz must have the same length, or one must be length 1",
         call. = FALSE)
  }

  # Recycle to common length
  n <- max(n_phon, n_freq)
  if (n_phon == 1) phon <- rep(phon, n)
  if (n_freq == 1) freq_hz <- rep(freq_hz, n)

  # Warn if outside reliable range
  if (any(phon < 20, na.rm = TRUE)) {
    n_low <- sum(phon < 20, na.rm = TRUE)
    warning(n_low, " value(s) below 20 phon. Results are informative only ",
            "(near hearing threshold).", call. = FALSE)
  }

  # Check upper limits per ISO 226
  for (i in seq_along(phon)) {
    f <- freq_hz[i]
    upper_limit <- if (f >= 5000) 80 else 90

    if (!is.na(phon[i]) && phon[i] > upper_limit) {
      warning("Frequency ", f, " Hz at ", round(phon[i], 1),
              " phon exceeds reliable range (>", upper_limit,
              " phon). Limited experimental data.", call. = FALSE)
    }
  }

  # Get parameters for each frequency
  params <- .get_iso226_params(freq_hz)

  # Extract parameters
  alpha_f <- params$alpha_f
  L_U <- params$L_U
  T_f <- params$T_f
  alpha_r <- params$alpha_r
  T_r <- params$T_r

  # Implement ISO 226:2023 Formula (1)
  # L_f = (10/alpha_f) * lg{ [(p0/pa)^2]^(alpha_r-alpha_f) *
  #       [10^(alpha_f*L_N/10phon) - 10^(alpha_f*T_r/10dB)] +
  #       10^(alpha_f*T_f+L_U/10dB) } dB - L_U

  # Simplified form from the standard:
  # L_f = (10/alpha_f) * lg{ (4*10^-10)^(0.3-alpha_f) *
  #       [10^(0.03*L_N/phon) - 10^0.072] + 10^(alpha_f*T_f+L_U/10dB) } dB - L_U

  # Calculate each term
  term1 <- (4 * 10^(-10))^(0.3 - alpha_f)
  term2 <- 10^(0.03 * phon)
  term3 <- 10^0.072
  term4 <- 10^((alpha_f * T_f + L_U) / 10)

  bracket <- term1 * (term2 - term3) + term4

  # Calculate SPL in dB
  L_f <- (10 / alpha_f) * log10(bracket) - L_U

  return(L_f)
}
