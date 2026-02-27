#' ISO 532 Sone (Loudness) Conversions
#'
#' @description
#' Functions for converting between phon (loudness level) and sone (loudness)
#' according to ISO 532 standards. Sone is a linear unit of perceived loudness,
#' where 1 sone = 40 phons by definition.
#'
#' @details
#' ## Standards
#'
#' Two ISO 532 methods are supported:
#'
#' - **ISO 532-1 (Zwicker method)**: Uses piecewise mathematical formulas
#' - **ISO 532-2 (Moore-Glasberg method)**: Uses lookup table with interpolation
#'
#' ## Zwicker Method (ISO 532-1)
#'
#' The Zwicker method uses different formulas depending on loudness level:
#'
#' **For sone < 1 (phon < 40):**
#' \deqn{L_N = 40 \times S^{0.35}}
#'
#' **For sone ≥ 1 (phon ≥ 40):**
#' \deqn{L_N = 40 + 10 \times \log_2(S)}
#'
#' Inverse formulas:
#'
#' **For phon < 40:**
#' \deqn{S = (L_N / 40)^{1/0.35} = (L_N / 40)^{2.857}}
#'
#' **For phon ≥ 40:**
#' \deqn{S = 2^{(L_N - 40) / 10}}
#'
#' ## Moore-Glasberg Method (ISO 532-2)
#'
#' Uses a lookup table with 23 reference points from 0.001 sone (0 phon)
#' to 337.6 sone (120 phon), with log-linear interpolation.
#'
#' ## Key Properties
#'
#' - **Reference:** 1 sone = 40 phons = 40 dB SPL at 1 kHz
#' - **Doubling:** Each 10 phon increase ≈ doubles loudness in sones
#' - **Linear perception:** Sones represent linear loudness (2 sones = 2× louder)
#' - **Stevens' power law:** Loudness ∝ intensity^0.3
#'
#' @name iso532-sone
#' @references
#' ISO 532-1:2017 - Acoustics - Methods for calculating loudness - Part 1: Zwicker method
#'
#' ISO 532-2:2017 - Acoustics - Methods for calculating loudness - Part 2: Moore-Glasberg method
#'
#' Stevens, S. S. (1936). "A scale for the measurement of a psychological magnitude: loudness".
#' Psychological Review. 43 (5): 405–416.
NULL

#' Moore-Glasberg Lookup Table (ISO 532-2)
#'
#' Reference table for sone-phon conversion using the Moore-Glasberg method.
#' Contains 23 calibration points from 0 to 120 phon.
#'
#' @format A data frame with 23 rows and 2 columns:
#' \describe{
#'   \item{phon}{Loudness level in phons (0 to 120)}
#'   \item{sone}{Loudness in sones (0.001 to 337.6)}
#' }
#'
#' @references
#' MATLAB Audio Toolbox documentation for sone2phon (ISO 532-2 reference table)
#'
#' @keywords internal
.moore_glasberg_table <- data.frame(
  phon = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
           65, 70, 75, 80, 85, 90, 95, 100, 110, 120),
  sone = c(0.001, 0.024, 0.061, 0.11, 0.18, 0.28, 0.44, 0.67, 1.00, 1.46,
           2.14, 3.14, 4.59, 6.73, 9.87, 14.47, 21.22, 31.12, 45.65,
           66.94, 98.18, 215.27, 337.60),
  stringsAsFactors = FALSE
)

#' Convert Phon to Sone
#'
#' Convert loudness level (phon) to loudness (sone) using ISO 532 methods.
#'
#' @param phon Numeric; loudness level in phons. Must be non-negative.
#' @param method Character; conversion method. One of:
#'   - `"zwicker"` (default): ISO 532-1 Zwicker method with piecewise formula
#'   - `"moore-glasberg"`: ISO 532-2 Moore-Glasberg method with lookup table
#'
#' @return Numeric vector of loudness values in sones
#'
#' @details
#' The Zwicker method (ISO 532-1) uses:
#' - For phon < 40: `sone = (phon / 40)^(1/0.35)`
#' - For phon ≥ 40: `sone = 2^((phon - 40) / 10)`
#'
#' The Moore-Glasberg method (ISO 532-2) uses a lookup table with 23 reference
#' points and log-linear interpolation.
#'
#' By definition, 1 sone = 40 phons (a 1 kHz tone at 40 dB SPL).
#'
#' ## Vectorization
#'
#' This function is fully vectorized and can process vectors of phon values
#' efficiently.
#'
#' @examples
#' # Reference value
#' phon_to_sone(40)  # Returns 1.0 (by definition)
#'
#' # Doubling property: +10 phon ≈ 2× loudness
#' phon_to_sone(50)  # Returns ~2.0
#' phon_to_sone(60)  # Returns ~4.0
#'
#' # Low levels (below 40 phon)
#' phon_to_sone(20)  # Returns ~0.15
#'
#' # Vectorized
#' phon_to_sone(c(20, 40, 60, 80))
#'
#' # Moore-Glasberg method
#' phon_to_sone(40, method = "moore-glasberg")
#'
#' @seealso [sone_to_phon()], [db_and_hz_to_phon()], [phon_and_hz_to_db()]
phon_to_sone <- function(phon, method = c("zwicker", "moore-glasberg")) {
  # Input validation
  if (!is.numeric(phon)) {
    stop("phon must be numeric", call. = FALSE)
  }

  if (any(is.na(phon))) {
    stop("phon contains NA values", call. = FALSE)
  }

  if (any(phon < 0)) {
    stop("phon must be non-negative (≥ 0)", call. = FALSE)
  }

  # Match method
  method <- match.arg(method)

  if (method == "zwicker") {
    # Zwicker method (ISO 532-1): piecewise formula
    sone <- numeric(length(phon))

    # Below 40 phon: S = (L_N / 40)^(1/0.35)
    below_40 <- phon < 40
    sone[below_40] <- (phon[below_40] / 40)^(1 / 0.35)

    # At or above 40 phon: S = 2^((L_N - 40) / 10)
    above_40 <- phon >= 40
    sone[above_40] <- 2^((phon[above_40] - 40) / 10)

  } else {
    # Moore-Glasberg method (ISO 532-2): lookup table with interpolation
    # Use log-linear interpolation for better accuracy across orders of magnitude

    # Get reference table
    ref_table <- .moore_glasberg_table

    # Log-transform sone values for interpolation
    log_sone_ref <- log10(ref_table$sone)

    # Interpolate in log space
    log_sone <- stats::approx(
      x = ref_table$phon,
      y = log_sone_ref,
      xout = phon,
      method = "linear",
      rule = 2  # Extrapolate for values outside range
    )$y

    # Back-transform to linear sone values
    sone <- 10^log_sone
  }

  # Add warning for very high values
  if (any(phon > 120)) {
    warning("Some phon values exceed 120 (maximum typical range). ",
            "Results are extrapolations.",
            call. = FALSE)
  }

  return(sone)
}

#' Convert Sone to Phon
#'
#' Convert loudness (sone) to loudness level (phon) using ISO 532 methods.
#'
#' @param sone Numeric; loudness in sones. Must be positive.
#' @param method Character; conversion method. One of:
#'   - `"zwicker"` (default): ISO 532-1 Zwicker method with piecewise formula
#'   - `"moore-glasberg"`: ISO 532-2 Moore-Glasberg method with lookup table
#'
#' @return Numeric vector of loudness level values in phons
#'
#' @details
#' The Zwicker method (ISO 532-1) uses:
#' - For sone < 1: `phon = 40 × sone^0.35`
#' - For sone ≥ 1: `phon = 40 + 10 × log₂(sone)`
#'
#' The Moore-Glasberg method (ISO 532-2) uses a lookup table with 23 reference
#' points and log-linear interpolation.
#'
#' By definition, 1 sone = 40 phons (a 1 kHz tone at 40 dB SPL).
#'
#' ## Vectorization
#'
#' This function is fully vectorized and can process vectors of sone values
#' efficiently.
#'
#' @examples
#' # Reference value
#' sone_to_phon(1)  # Returns 40 (by definition)
#'
#' # Doubling property: 2× loudness ≈ +10 phon
#' sone_to_phon(2)  # Returns ~50
#' sone_to_phon(4)  # Returns ~60
#'
#' # Low levels (below 1 sone)
#' sone_to_phon(0.15)  # Returns ~20
#'
#' # Vectorized
#' sone_to_phon(c(0.15, 1, 4, 16))
#'
#' # Moore-Glasberg method
#' sone_to_phon(1, method = "moore-glasberg")
#'
#' @seealso [phon_to_sone()], [db_and_hz_to_phon()], [phon_and_hz_to_db()]
sone_to_phon <- function(sone, method = c("zwicker", "moore-glasberg")) {
  # Input validation
  if (!is.numeric(sone)) {
    stop("sone must be numeric", call. = FALSE)
  }

  if (any(is.na(sone))) {
    stop("sone contains NA values", call. = FALSE)
  }

  if (any(sone <= 0)) {
    stop("sone must be positive (> 0)", call. = FALSE)
  }

  # Match method
  method <- match.arg(method)

  if (method == "zwicker") {
    # Zwicker method (ISO 532-1): piecewise formula
    phon <- numeric(length(sone))

    # Below 1 sone: L_N = 40 × S^0.35
    below_1 <- sone < 1
    phon[below_1] <- 40 * (sone[below_1]^0.35)

    # At or above 1 sone: L_N = 40 + 10 × log₂(S)
    above_1 <- sone >= 1
    phon[above_1] <- 40 + 10 * log2(sone[above_1])

  } else {
    # Moore-Glasberg method (ISO 532-2): lookup table with interpolation
    # Use log-linear interpolation for better accuracy

    # Get reference table
    ref_table <- .moore_glasberg_table

    # Log-transform sone values for interpolation
    log_sone_ref <- log10(ref_table$sone)
    log_sone_input <- log10(sone)

    # Interpolate in log space
    phon <- stats::approx(
      x = log_sone_ref,
      y = ref_table$phon,
      xout = log_sone_input,
      method = "linear",
      rule = 2  # Extrapolate for values outside range
    )$y
  }

  # Add warning for very high values
  if (any(sone > 340)) {
    warning("Some sone values exceed 340 (~120 phon maximum typical range). ",
            "Results are extrapolations.",
            call. = FALSE)
  }

  return(phon)
}

#' Convert dB and Hz Directly to Sone
#'
#' Convenience function to convert sound pressure level (dB) and frequency (Hz)
#' directly to loudness (sone), combining ISO 226:2023 and ISO 532 conversions.
#'
#' @param spl_db Numeric; sound pressure level in dB SPL (re 20 μPa)
#' @param freq_hz Numeric; frequency in Hz (20-12500 Hz)
#' @param method Character; ISO 532 method for phon-to-sone conversion.
#'   One of `"zwicker"` (default) or `"moore-glasberg"`.
#'
#' @return Numeric vector of loudness values in sones
#'
#' @details
#' This function combines two conversions:
#' 1. (SPL, frequency) → phon using ISO 226:2023 equal-loudness contours
#' 2. phon → sone using ISO 532 Zwicker or Moore-Glasberg method
#'
#' Equivalent to: `phon_to_sone(db_and_hz_to_phon(spl_db, freq_hz), method)`
#'
#' @examples
#' # At 1 kHz reference
#' db_and_hz_to_sone(40, 1000)  # Returns 1.0 (by definition)
#' db_and_hz_to_sone(50, 1000)  # Returns ~2.0 (+10 dB ≈ 2× loudness)
#'
#' # Low frequency requires more SPL for same loudness
#' db_and_hz_to_sone(60, 100)   # Lower sone value
#' db_and_hz_to_sone(60, 1000)  # Higher sone value
#'
#' @seealso [db_and_hz_to_phon()], [phon_to_sone()]
db_and_hz_to_sone <- function(spl_db, freq_hz, method = c("zwicker", "moore-glasberg")) {
  # Step 1: Convert (dB, Hz) to phon using ISO 226:2023
  phon <- db_and_hz_to_phon(spl_db, freq_hz)

  # Step 2: Convert phon to sone using ISO 532
  sone <- phon_to_sone(phon, method = method)

  return(sone)
}

#' Convert Sone and Hz to dB
#'
#' Convenience function to convert loudness (sone) and frequency (Hz) to
#' sound pressure level (dB), combining ISO 532 and ISO 226:2023 conversions.
#'
#' @param sone Numeric; loudness in sones. Must be positive.
#' @param freq_hz Numeric; frequency in Hz (20-12500 Hz)
#' @param method Character; ISO 532 method for sone-to-phon conversion.
#'   One of `"zwicker"` (default) or `"moore-glasberg"`.
#'
#' @return Numeric vector of sound pressure level in dB SPL (re 20 μPa)
#'
#' @details
#' This function combines two conversions:
#' 1. sone → phon using ISO 532 Zwicker or Moore-Glasberg method
#' 2. (phon, frequency) → SPL using ISO 226:2023 equal-loudness contours
#'
#' Equivalent to: `phon_and_hz_to_db(sone_to_phon(sone, method), freq_hz)`
#'
#' @examples
#' # At 1 kHz reference
#' sone_and_hz_to_db(1, 1000)  # Returns 40 dB (by definition)
#' sone_and_hz_to_db(2, 1000)  # Returns ~50 dB (2× loudness ≈ +10 dB)
#'
#' # Same loudness at different frequencies requires different SPL
#' sone_and_hz_to_db(2, 100)   # Higher dB (low freq needs more SPL)
#' sone_and_hz_to_db(2, 1000)  # Reference
#' sone_and_hz_to_db(2, 4000)  # Lower dB (high freq more sensitive)
#'
#' @seealso [sone_to_phon()], [phon_and_hz_to_db()]
sone_and_hz_to_db <- function(sone, freq_hz, method = c("zwicker", "moore-glasberg")) {
  # Step 1: Convert sone to phon using ISO 532
  phon <- sone_to_phon(sone, method = method)

  # Step 2: Convert (phon, Hz) to dB using ISO 226:2023
  spl_db <- phon_and_hz_to_db(phon, freq_hz)

  return(spl_db)
}
