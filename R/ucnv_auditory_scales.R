#' Convert Frequency to Bark Scale
#'
#' @description
#' Converts frequency values from Hz to the Bark scale, a psychoacoustic scale
#' that approximates the human ear's critical bands. The Bark scale is linear
#' up to 1 Bark per critical band (24 critical bands total).
#'
#' @param freq Numeric vector or units object with frequency values in Hz.
#'   If a units object is provided, it will be converted to Hz.
#' @param method Character string specifying the conversion formula to use.
#'   Options are:
#'   \itemize{
#'     \item \code{"traunmuller"} (default): Traunmüller (1990) approximation
#'     \item \code{"zwicker"}: Zwicker's formula
#'     \item \code{"wang"}: Wang, Sekey & Gersho (1992) formula
#'   }
#' @param as_units Logical. If TRUE, returns a units object with "Bark" units.
#'   If FALSE, returns a plain numeric vector. Default is TRUE if the units
#'   package is available and input has units.
#'
#' @details
#' Three different approximation formulas are available:
#'
#' **Traunmüller (1990)** (default):
#' \deqn{Bark = \frac{26.81 \cdot f}{1960 + f} - 0.53}
#'
#' Valid range: approximately 20-15500 Hz. Values below 2 Bark are adjusted
#' by \eqn{Bark_{adjusted} = Bark + 0.15 \cdot (2 - Bark)}.
#'
#' **Zwicker**:
#' \deqn{Bark = 13 \cdot \arctan(0.00076 \cdot f) + 3.5 \cdot \arctan\left(\frac{f}{7500}\right)^2}
#'
#' **Wang, Sekey & Gersho (1992)**:
#' \deqn{Bark = 6 \cdot \sinh^{-1}\left(\frac{f}{600}\right)}
#'
#' @return
#' If \code{as_units = TRUE} and the units package is available: a units object
#' with "Bark" units. Otherwise, a numeric vector of Bark values.
#'
#' @references
#' Traunmüller, H. (1990). Analytical expressions for the tonotopic sensory
#' scale. Journal of the Acoustical Society of America, 88, 97-100.
#'
#' Zwicker, E. (1961). Subdivision of the audible frequency range into critical
#' bands. Journal of the Acoustical Society of America, 33, 248.
#'
#' Wang, S., Sekey, A., & Gersho, A. (1992). An objective measure for
#' predicting subjective quality of speech coders. IEEE Journal on Selected
#' Areas in Communications, 10(5), 819-829.
#'
#' @examples
#' \dontrun{
#' # Basic usage with numeric input
#' ucnv_hz_to_bark(1000)
#' ucnv_hz_to_bark(c(100, 500, 1000, 2000, 4000))
#'
#' # With units package
#' library(units)
#' freq <- set_units(c(100, 500, 1000, 2000), Hz)
#' ucnv_hz_to_bark(freq)
#'
#' # Different methods
#' ucnv_hz_to_bark(1000, method = "zwicker")
#' ucnv_hz_to_bark(1000, method = "wang")
#'
#' # Without units in output
#' ucnv_hz_to_bark(1000, as_units = FALSE)
#' }
#'
#' @export
ucnv_hz_to_bark <- function(freq, method = c("traunmuller", "zwicker", "wang"),
                       as_units = NULL) {
  method <- match.arg(method)

  # Handle units if present
  has_units <- FALSE
  if (inherits(freq, "units")) {
    has_units <- TRUE
    # Convert to Hz (numeric)
    freq_hz <- as.numeric(units::set_units(freq, Hz))
  } else {
    freq_hz <- as.numeric(freq)
  }

  # Validate input
  if (any(freq_hz < 0, na.rm = TRUE)) {
    stop("Frequency values must be non-negative")
  }

  # Apply conversion formula
  bark <- switch(method,
    traunmuller = {
      # Traunmüller (1990)
      b <- ((26.81 * freq_hz) / (1960 + freq_hz)) - 0.53
      # Correction for low frequencies (< 2 Bark)
      low_freq <- b < 2
      b[low_freq] <- b[low_freq] + 0.15 * (2 - b[low_freq])
      b
    },
    zwicker = {
      # Zwicker's formula
      13 * atan(0.00076 * freq_hz) + 3.5 * atan((freq_hz / 7500)^2)
    },
    wang = {
      # Wang, Sekey & Gersho (1992)
      6 * asinh(freq_hz / 600)
    }
  )

  # Determine if we should return units
  if (is.null(as_units)) {
    as_units <- has_units && requireNamespace("units", quietly = TRUE)
  }

  # Add units if requested
  if (as_units) {
    if (!requireNamespace("units", quietly = TRUE)) {
      warning("units package not available, returning plain numeric vector")
      return(bark)
    }

    # Ensure Bark unit is installed
    .ensure_bark_unit()

    # Return with units by multiplying with unit
    bark_unit <- units::as_units("Bark")
    return(bark * bark_unit)
  }

  bark
}


#' Convert Bark Scale to Frequency
#'
#' @description
#' Converts values from the Bark scale to frequency in Hz. This is the inverse
#' of \code{\link{ucnv_hz_to_bark}}.
#'
#' @param bark Numeric vector or units object with Bark values.
#'   If a units object is provided with "Bark" units, it will be converted
#'   to numeric.
#' @param method Character string specifying which inverse formula to use.
#'   Must match the forward conversion method. Options are:
#'   \itemize{
#'     \item \code{"traunmuller"} (default): Inverse of Traunmüller (1990)
#'     \item \code{"zwicker"}: Inverse of Zwicker's formula (numerical)
#'     \item \code{"wang"}: Inverse of Wang, Sekey & Gersho (1992)
#'   }
#' @param as_units Logical. If TRUE, returns a units object with "Hz" units.
#'   If FALSE, returns a plain numeric vector. Default is TRUE if the units
#'   package is available and input has units.
#'
#' @details
#' The inverse formulas are derived analytically where possible:
#'
#' **Traunmüller (1990)** (default):
#' \deqn{f = \frac{1960 \cdot (Bark + 0.53)}{26.81 - (Bark + 0.53)}}
#'
#' Note: For Bark < 2, the input is adjusted to account for the low-frequency
#' correction in the forward transform.
#'
#' **Wang, Sekey & Gersho (1992)**:
#' \deqn{f = 600 \cdot \sinh\left(\frac{Bark}{6}\right)}
#'
#' **Zwicker**: Uses numerical root-finding (uniroot) since the inverse
#' is not analytically solvable.
#'
#' @return
#' If \code{as_units = TRUE} and the units package is available: a units object
#' with "Hz" units. Otherwise, a numeric vector of frequencies in Hz.
#'
#' @seealso \code{\link{ucnv_hz_to_bark}}
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' ucnv_bark_to_hz(10)
#' ucnv_bark_to_hz(c(5, 10, 15, 20))
#'
#' # Round trip conversion
#' freq <- 1000
#' bark <- ucnv_hz_to_bark(freq)
#' ucnv_bark_to_hz(bark)  # Should return ~1000
#'
#' # With units package
#' library(units)
#' b <- set_units(10, Bark)
#' ucnv_bark_to_hz(b)
#'
#' # Different methods
#' ucnv_bark_to_hz(10, method = "wang")
#' }
#'
#' @export
ucnv_bark_to_hz <- function(bark, method = c("traunmuller", "zwicker", "wang"),
                       as_units = NULL) {
  method <- match.arg(method)

  # Handle units if present
  has_units <- FALSE
  if (inherits(bark, "units")) {
    has_units <- TRUE
    # Just extract the numeric value - assume it's already in Bark
    bark_val <- as.numeric(bark)
  } else {
    bark_val <- as.numeric(bark)
  }

  # Validate input
  if (any(bark_val < 0 | bark_val > 24, na.rm = TRUE)) {
    warning("Bark values outside typical range [0, 24] detected")
  }

  # Apply inverse conversion formula
  freq_hz <- switch(method,
    traunmuller = {
      # Inverse of Traunmüller (1990)
      # First, undo low-frequency correction
      b_corrected <- bark_val
      low_freq <- bark_val < 2
      b_corrected[low_freq] <- (bark_val[low_freq] - 0.3) / 1.15

      # Now invert: Bark = (26.81*f)/(1960+f) - 0.53
      # => Bark + 0.53 = (26.81*f)/(1960+f)
      # => (Bark + 0.53)*(1960+f) = 26.81*f
      # => 1960*(Bark + 0.53) + f*(Bark + 0.53) = 26.81*f
      # => 1960*(Bark + 0.53) = f*(26.81 - (Bark + 0.53))
      # => f = 1960*(Bark + 0.53) / (26.81 - (Bark + 0.53))

      numerator <- 1960 * (b_corrected + 0.53)
      denominator <- 26.81 - (b_corrected + 0.53)
      numerator / denominator
    },
    zwicker = {
      # Zwicker's inverse requires numerical solution
      # Use vectorized approach
      sapply(bark_val, function(b) {
        if (is.na(b)) return(NA_real_)

        # Define the forward function minus target
        f <- function(freq) {
          13 * atan(0.00076 * freq) + 3.5 * atan((freq / 7500)^2) - b
        }

        # Use uniroot to find the frequency
        # Valid frequency range is roughly 20-15500 Hz
        tryCatch({
          uniroot(f, interval = c(20, 15500), tol = 0.01)$root
        }, error = function(e) {
          warning(sprintf("Failed to invert Bark=%g: %s", b, e$message))
          NA_real_
        })
      })
    },
    wang = {
      # Inverse of Wang, Sekey & Gersho (1992)
      # Bark = 6 * asinh(f/600)
      # => Bark/6 = asinh(f/600)
      # => sinh(Bark/6) = f/600
      # => f = 600 * sinh(Bark/6)
      600 * sinh(bark_val / 6)
    }
  )

  # Determine if we should return units
  if (is.null(as_units)) {
    as_units <- has_units && requireNamespace("units", quietly = TRUE)
  }

  # Add units if requested
  if (as_units) {
    if (!requireNamespace("units", quietly = TRUE)) {
      warning("units package not available, returning plain numeric vector")
      return(freq_hz)
    }

    # Return with units
    return(units::set_units(freq_hz, Hz))
  }

  freq_hz
}


#' Ensure Bark Unit is Installed
#'
#' @description
#' Internal helper function that ensures the "Bark" unit is registered
#' with the units package. This is called automatically by ucnv_hz_to_bark
#' and ucnv_bark_to_hz when needed.
#'
#' @details
#' The Bark scale is defined as a dimensionless psychoacoustic unit.
#' It represents critical bands in human auditory perception.
#'
#' @return NULL (called for side effects)
#' @keywords internal
#' @noRd
.ensure_bark_unit <- function() {
  # Check if already installed by trying to create the unit
  tryCatch({
    test <- units::as_units("Bark")
    # If we get here without error, it's already installed
    return(invisible(NULL))
  }, error = function(e) {
    # Not installed, install it
    tryCatch({
      units::install_unit("Bark", "", "bark")
    }, error = function(e2) {
      invisible(NULL)
    })
  })

  invisible(NULL)
}


#' Initialize Psychoacoustic Units
#'
#' @description
#' Registers psychoacoustic units (Bark, ERB, mel, semitone) with the units package.
#' This function is called automatically when superassp is loaded if the
#' units package is available.
#'
#' @details
#' Registers:
#' \itemize{
#'   \item Bark - Critical band rate scale
#'   \item ERB - Equivalent Rectangular Bandwidth rate scale
#'   \item mel - Mel scale (perceptual pitch)
#'   \item semitone - Musical semitone (ST)
#' }
#'
#' @return NULL (called for side effects)
#' @keywords internal
.onLoad_psychoacoustic_units <- function() {
  if (requireNamespace("units", quietly = TRUE)) {
    .ensure_bark_unit()
    .ensure_erb_unit()
    .ensure_mel_unit()
    .ensure_semitone_unit()
  }
}

# ============================================================================
# ERB Scale (Equivalent Rectangular Bandwidth)
# ============================================================================

#' Convert Frequency to ERB-rate Scale
#'
#' @description
#' Converts frequency values from Hz to the ERB-rate scale (also called ERBS or Cams).
#' The ERB-rate scale is based on the Equivalent Rectangular Bandwidth of auditory
#' filters and provides a perceptually uniform frequency scale.
#'
#' @param freq Numeric vector or units object with frequency values in Hz.
#'   If a units object is provided, it will be converted to Hz.
#' @param method Character string specifying the conversion formula to use.
#'   Options are:
#'   \itemize{
#'     \item \code{"glasberg1990"} (default): Glasberg & Moore (1990) linear approximation
#'     \item \code{"moore1983"}: Moore & Glasberg (1983) polynomial approximation
#'   }
#' @param as_units Logical. If TRUE, returns a units object with "ERB" units.
#'   If FALSE, returns a plain numeric vector. Default is TRUE if the units
#'   package is available and input has units.
#'
#' @details
#' Two different approximation formulas are available:
#'
#' **Glasberg & Moore (1990)** (default):
#' \deqn{ERBS(f) = 21.4 \cdot \log_{10}(1 + 0.00437 \cdot f)}
#'
#' Valid range: 100-10,000 Hz. This is the most commonly used formula.
#'
#' **Moore & Glasberg (1983)**:
#' \deqn{ERBS(f) = 11.17 \cdot \ln\left(\frac{f + 312}{f + 14675}\right) + 43.0}
#'
#' where f is in Hz. Valid range: 100-6500 Hz.
#'
#' @return
#' If \code{as_units = TRUE} and the units package is available: a units object
#' with "ERB" units. Otherwise, a numeric vector of ERB-rate values.
#'
#' @references
#' Glasberg, B. R., & Moore, B. C. (1990). Derivation of auditory filter shapes
#' from notched-noise data. Hearing Research, 47(1-2), 103-138.
#'
#' Moore, B. C., & Glasberg, B. R. (1983). Suggested formulae for calculating
#' auditory-filter bandwidths and excitation patterns. The Journal of the
#' Acoustical Society of America, 74(3), 750-753.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' ucnv_hz_to_erb(1000)
#' ucnv_hz_to_erb(c(100, 500, 1000, 2000, 4000))
#'
#' # With units package
#' library(units)
#' freq <- set_units(1000, Hz)
#' ucnv_hz_to_erb(freq)
#'
#' # Different methods
#' ucnv_hz_to_erb(1000, method = "moore1983")
#' }
#'
#' @export
ucnv_hz_to_erb <- function(freq, method = c("glasberg1990", "moore1983"),
                      as_units = NULL) {
  method <- match.arg(method)

  # Handle units if present
  has_units <- FALSE
  if (inherits(freq, "units")) {
    has_units <- TRUE
    freq_hz <- as.numeric(units::set_units(freq, Hz))
  } else {
    freq_hz <- as.numeric(freq)
  }

  # Validate input
  if (any(freq_hz < 0, na.rm = TRUE)) {
    stop("Frequency values must be non-negative")
  }

  # Apply conversion formula
  erb <- switch(method,
    glasberg1990 = {
      # Glasberg & Moore (1990) - most common
      21.4 * log10(1 + 0.00437 * freq_hz)
    },
    moore1983 = {
      # Moore & Glasberg (1983) - polynomial form
      # Convert to kHz for this formula
      f_khz <- freq_hz / 1000
      11.17 * log((f_khz + 0.312) / (f_khz + 14.675)) + 43.0
    }
  )

  # Determine if we should return units
  if (is.null(as_units)) {
    as_units <- has_units && requireNamespace("units", quietly = TRUE)
  }

  # Add units if requested
  if (as_units) {
    if (!requireNamespace("units", quietly = TRUE)) {
      warning("units package not available, returning plain numeric vector")
      return(erb)
    }

    # Ensure ERB unit is installed
    .ensure_erb_unit()

    # Return with units
    erb_unit <- units::as_units("ERB")
    return(erb * erb_unit)
  }

  erb
}


#' Convert ERB-rate Scale to Frequency
#'
#' @description
#' Converts values from the ERB-rate scale to frequency in Hz.
#' This is the inverse of \code{\link{ucnv_hz_to_erb}}.
#'
#' @param erb Numeric vector or units object with ERB-rate values.
#' @param method Character string specifying which inverse formula to use.
#'   Must match the forward conversion method.
#' @param as_units Logical. If TRUE, returns a units object with "Hz" units.
#'
#' @return
#' If \code{as_units = TRUE}: a units object with "Hz" units.
#' Otherwise, a numeric vector of frequencies in Hz.
#'
#' @seealso \code{\link{ucnv_hz_to_erb}}
#'
#' @examples
#' \dontrun{
#' # Round trip
#' erb <- ucnv_hz_to_erb(1000)
#' ucnv_erb_to_hz(erb)  # Should return ~1000
#' }
#'
#' @export
ucnv_erb_to_hz <- function(erb, method = c("glasberg1990", "moore1983"),
                      as_units = NULL) {
  method <- match.arg(method)

  # Handle units if present
  has_units <- FALSE
  if (inherits(erb, "units")) {
    has_units <- TRUE
    erb_val <- as.numeric(erb)
  } else {
    erb_val <- as.numeric(erb)
  }

  # Apply inverse conversion formula
  freq_hz <- switch(method,
    glasberg1990 = {
      # Inverse of Glasberg & Moore (1990)
      # ERBS = 21.4 * log10(1 + 0.00437 * f)
      # => 10^(ERBS/21.4) = 1 + 0.00437 * f
      # => f = (10^(ERBS/21.4) - 1) / 0.00437
      (10^(erb_val / 21.4) - 1) / 0.00437
    },
    moore1983 = {
      # Inverse of Moore & Glasberg (1983)
      # ERBS = 11.17 * ln((f_khz + 0.312)/(f_khz + 14.675)) + 43.0
      # This requires numerical solution
      sapply(erb_val, function(e) {
        if (is.na(e)) return(NA_real_)

        f <- function(freq_khz) {
          11.17 * log((freq_khz + 0.312) / (freq_khz + 14.675)) + 43.0 - e
        }

        tryCatch({
          uniroot(f, interval = c(0.02, 20), tol = 0.001)$root * 1000
        }, error = function(err) {
          warning(sprintf("Failed to invert ERB=%g: %s", e, err$message))
          NA_real_
        })
      })
    }
  )

  # Determine if we should return units
  if (is.null(as_units)) {
    as_units <- has_units && requireNamespace("units", quietly = TRUE)
  }

  # Add units if requested
  if (as_units) {
    if (!requireNamespace("units", quietly = TRUE)) {
      warning("units package not available, returning plain numeric vector")
      return(freq_hz)
    }
    return(units::set_units(freq_hz, Hz))
  }

  freq_hz
}


# ============================================================================
# Mel Scale
# ============================================================================

#' Convert Frequency to Mel Scale
#'
#' @description
#' Converts frequency values from Hz to the mel scale, a perceptual pitch scale
#' where equal distances sound equally different to human listeners.
#'
#' @param freq Numeric vector or units object with frequency values in Hz.
#' @param method Character string specifying the conversion formula to use.
#'   Options are:
#'   \itemize{
#'     \item \code{"htk"} (default): HTK formula using base-10 logarithm
#'     \item \code{"slaney"}: Slaney's formula (linear below 1000 Hz)
#'   }
#' @param as_units Logical. If TRUE, returns a units object with "mel" units.
#'
#' @details
#' Two formulas are available:
#'
#' **HTK formula** (default):
#' \deqn{mel = 2595 \cdot \log_{10}(1 + f/700)}
#'
#' This gives exactly 1000 mels at 1000 Hz.
#'
#' **Slaney formula**:
#' Linear below 1000 Hz, logarithmic above. Used in some audio processing libraries.
#'
#' @return
#' If \code{as_units = TRUE}: a units object with "mel" units.
#' Otherwise, a numeric vector of mel values.
#'
#' @references
#' O'Shaughnessy, D. (1987). Speech Communication: Human and Machine.
#' Addison-Wesley.
#'
#' @examples
#' \dontrun{
#' # 1000 Hz = 1000 mels (approximately)
#' ucnv_hz_to_mel(1000)
#'
#' # Vector conversion
#' ucnv_hz_to_mel(c(100, 500, 1000, 2000, 4000))
#'
#' # With units
#' library(units)
#' freq <- set_units(1000, Hz)
#' ucnv_hz_to_mel(freq)
#' }
#'
#' @export
ucnv_hz_to_mel <- function(freq, method = c("htk", "slaney"), as_units = NULL) {
  method <- match.arg(method)

  # Handle units if present
  has_units <- FALSE
  if (inherits(freq, "units")) {
    has_units <- TRUE
    freq_hz <- as.numeric(units::set_units(freq, Hz))
  } else {
    freq_hz <- as.numeric(freq)
  }

  # Validate input
  if (any(freq_hz < 0, na.rm = TRUE)) {
    stop("Frequency values must be non-negative")
  }

  # Apply conversion formula
  mel <- switch(method,
    htk = {
      # HTK formula (most common)
      2595 * log10(1 + freq_hz / 700)
    },
    slaney = {
      # Slaney formula - linear below 1000 Hz, log above
      f_min <- 0
      f_sp <- 200 / 3  # ~66.67
      min_log_hz <- 1000
      min_log_mel <- (min_log_hz - f_min) / f_sp
      step_log <- 0.068751777  # log(6.4) / 27
      logstep <- log(6.4) / 27

      ifelse(freq_hz < min_log_hz,
             # Linear part
             (freq_hz - f_min) / f_sp,
             # Logarithmic part
             min_log_mel + log(freq_hz / min_log_hz) / logstep
      )
    }
  )

  # Determine if we should return units
  if (is.null(as_units)) {
    as_units <- has_units && requireNamespace("units", quietly = TRUE)
  }

  # Add units if requested
  if (as_units) {
    if (!requireNamespace("units", quietly = TRUE)) {
      warning("units package not available, returning plain numeric vector")
      return(mel)
    }

    .ensure_mel_unit()
    mel_unit <- units::as_units("mel")
    return(mel * mel_unit)
  }

  mel
}


#' Convert Mel Scale to Frequency
#'
#' @description
#' Converts values from the mel scale to frequency in Hz.
#' This is the inverse of \code{\link{ucnv_hz_to_mel}}.
#'
#' @param mel Numeric vector or units object with mel values.
#' @param method Character string specifying which inverse formula to use.
#' @param as_units Logical. If TRUE, returns a units object with "Hz" units.
#'
#' @return
#' If \code{as_units = TRUE}: a units object with "Hz" units.
#' Otherwise, a numeric vector of frequencies in Hz.
#'
#' @seealso \code{\link{ucnv_hz_to_mel}}
#'
#' @examples
#' \dontrun{
#' # Round trip
#' mel <- ucnv_hz_to_mel(1000)
#' ucnv_mel_to_hz(mel)  # Should return ~1000
#' }
#'
#' @export
ucnv_mel_to_hz <- function(mel, method = c("htk", "slaney"), as_units = NULL) {
  method <- match.arg(method)

  # Handle units if present
  has_units <- FALSE
  if (inherits(mel, "units")) {
    has_units <- TRUE
    mel_val <- as.numeric(mel)
  } else {
    mel_val <- as.numeric(mel)
  }

  # Apply inverse conversion formula
  freq_hz <- switch(method,
    htk = {
      # Inverse of HTK formula
      # mel = 2595 * log10(1 + f/700)
      # => 10^(mel/2595) = 1 + f/700
      # => f = 700 * (10^(mel/2595) - 1)
      700 * (10^(mel_val / 2595) - 1)
    },
    slaney = {
      # Inverse of Slaney formula
      f_min <- 0
      f_sp <- 200 / 3
      min_log_hz <- 1000
      min_log_mel <- (min_log_hz - f_min) / f_sp
      logstep <- log(6.4) / 27

      ifelse(mel_val < min_log_mel,
             # Linear part
             f_min + f_sp * mel_val,
             # Logarithmic part
             min_log_hz * exp(logstep * (mel_val - min_log_mel))
      )
    }
  )

  # Determine if we should return units
  if (is.null(as_units)) {
    as_units <- has_units && requireNamespace("units", quietly = TRUE)
  }

  # Add units if requested
  if (as_units) {
    if (!requireNamespace("units", quietly = TRUE)) {
      warning("units package not available, returning plain numeric vector")
      return(freq_hz)
    }
    return(units::set_units(freq_hz, Hz))
  }

  freq_hz
}


# ============================================================================
# Semitone Scale (Musical)
# ============================================================================

#' Convert Frequency to Semitones
#'
#' @description
#' Converts frequency values to semitones relative to a reference frequency.
#' In equal temperament, each semitone represents a frequency ratio of 2^(1/12).
#'
#' @param freq Numeric vector or units object with frequency values in Hz.
#' @param ref_freq Reference frequency in Hz. Default is NULL, which uses the
#'   reference determined by \code{ref_source}. Can also be a units object.
#'   If specified, overrides \code{ref_source}.
#' @param ref_source Character string specifying the reference standard.
#'   Options are:
#'   \itemize{
#'     \item \code{"UEP83"}: UEP 1983 standard (110 Hz = A2/A₁ in Helmholtz notation)
#'     \item \code{"Praat"}: Praat convention (100 Hz arbitrary reference)
#'     \item \code{"A4"} (default): Concert pitch A4 = 440 Hz
#'   }
#'   Ignored if \code{ref_freq} is explicitly specified.
#' @param as_units Logical. If TRUE, returns a units object with "semitone" units.
#'
#' @details
#' The conversion formula is:
#' \deqn{ST = 12 \cdot \log_2(f / f_{ref})}
#'
#' where f is the frequency and f_ref is the reference frequency.
#'
#' **Reference Standards:**
#'
#' - **UEP 1983** (Schutte & Seidner): Uses 110 Hz (A2, or A₁ in Helmholtz notation)
#'   as the reference for voice range profiles (phonetograms). This standard is
#'   commonly used in clinical phoniatrics and voice assessment.
#'
#' - **Praat**: Uses 100 Hz as an arbitrary reference frequency for semitone-based
#'   pitch analysis. This provides a convenient round number for relative pitch
#'   measurements in speech analysis.
#'
#' - **A4** (default): Uses 440 Hz (A4, concert pitch) as reference, following
#'   standard musical tuning conventions.
#'
#' This gives the number of semitones above (positive) or below (negative) the
#' reference frequency.
#'
#' @return
#' If \code{as_units = TRUE}: a units object with "semitone" units.
#' Otherwise, a numeric vector of semitone values.
#'
#' @references
#' Schutte, H.K., & Seidner, W. (1983). Recommendation by the Union of European
#' Phoniatricians (UEP): Standardizing Voice Area Measurement/Phonetography.
#' Folia Phoniatrica, 35, 286-288.
#'
#' @examples
#' \dontrun{
#' # Using different reference standards
#' ucnv_hz_to_semitone(220, ref_source = "UEP83")   # 12 ST above 110 Hz
#' ucnv_hz_to_semitone(200, ref_source = "Praat")   # 12 ST above 100 Hz
#' ucnv_hz_to_semitone(880, ref_source = "A4")      # 12 ST above 440 Hz
#'
#' # Explicit reference frequency (overrides ref_source)
#' ucnv_hz_to_semitone(880, ref_freq = 440)
#'
#' # Musical notes relative to A4
#' ucnv_hz_to_semitone(c(440, 494, 523.25, 587.33, 659.25, 698.46, 783.99, 880))
#'
#' # With units
#' library(units)
#' freq <- set_units(880, Hz)
#' ucnv_hz_to_semitone(freq, ref_source = "A4")
#' }
#'
#' @export
ucnv_hz_to_semitone <- function(freq, ref_freq = NULL, ref_source = c("A4", "UEP83", "Praat"),
                            as_units = NULL) {
  ref_source <- match.arg(ref_source)

  # Determine reference frequency
  if (is.null(ref_freq)) {
    ref_freq <- switch(ref_source,
      "UEP83" = 110,    # A2 (A₁ in Helmholtz notation) - UEP 1983 standard
      "Praat" = 100,    # Praat's arbitrary reference for semitones
      "A4" = 440        # Concert pitch A4
    )
  }
  # Handle units if present
  has_units <- FALSE
  freq_names <- names(freq)  # Preserve names
  if (inherits(freq, "units")) {
    has_units <- TRUE
    freq_hz <- as.numeric(units::set_units(freq, Hz))
  } else {
    freq_hz <- as.numeric(freq)
  }

  # Handle reference frequency units
  if (inherits(ref_freq, "units")) {
    ref_freq_hz <- as.numeric(units::set_units(ref_freq, Hz))
  } else {
    ref_freq_hz <- as.numeric(ref_freq)
  }

  # Validate input
  if (any(freq_hz <= 0, na.rm = TRUE)) {
    stop("Frequency values must be positive")
  }
  if (ref_freq_hz <= 0) {
    stop("Reference frequency must be positive")
  }

  # Apply conversion formula
  # ST = 12 * log2(f / f_ref)
  semitone <- 12 * log2(freq_hz / ref_freq_hz)

  # Restore names
  if (!is.null(freq_names)) {
    names(semitone) <- freq_names
  }

  # Determine if we should return units
  if (is.null(as_units)) {
    as_units <- has_units && requireNamespace("units", quietly = TRUE)
  }

  # Add units if requested
  if (as_units) {
    if (!requireNamespace("units", quietly = TRUE)) {
      warning("units package not available, returning plain numeric vector")
      return(semitone)
    }

    .ensure_semitone_unit()
    st_unit <- units::as_units("semitone")
    return(semitone * st_unit)
  }

  semitone
}


#' Convert Semitones to Frequency
#'
#' @description
#' Converts semitone values to frequency in Hz relative to a reference frequency.
#' This is the inverse of \code{\link{ucnv_hz_to_semitone}}.
#'
#' @param semitone Numeric vector or units object with semitone values.
#' @param ref_freq Reference frequency in Hz. Default is NULL, which uses the
#'   reference determined by \code{ref_source}. Can also be a units object.
#'   If specified, overrides \code{ref_source}.
#' @param ref_source Character string specifying the reference standard.
#'   Options are:
#'   \itemize{
#'     \item \code{"UEP83"}: UEP 1983 standard (110 Hz = A2)
#'     \item \code{"Praat"}: Praat convention (100 Hz)
#'     \item \code{"A4"} (default): Concert pitch A4 = 440 Hz
#'   }
#'   Ignored if \code{ref_freq} is explicitly specified.
#' @param as_units Logical. If TRUE, returns a units object with "Hz" units.
#'
#' @details
#' The conversion formula is:
#' \deqn{f = f_{ref} \cdot 2^{ST/12}}
#'
#' See \code{\link{ucnv_hz_to_semitone}} for details on reference standards.
#'
#' @return
#' If \code{as_units = TRUE}: a units object with "Hz" units.
#' Otherwise, a numeric vector of frequencies in Hz.
#'
#' @seealso \code{\link{ucnv_hz_to_semitone}}
#'
#' @examples
#' \dontrun{
#' # Using different reference standards
#' ucnv_semitone_to_hz(12, ref_source = "UEP83")   # 220 Hz (octave above 110)
#' ucnv_semitone_to_hz(12, ref_source = "Praat")   # 200 Hz (octave above 100)
#' ucnv_semitone_to_hz(12, ref_source = "A4")      # 880 Hz (octave above 440)
#'
#' # Musical scale from A4
#' ucnv_semitone_to_hz(0:12, ref_freq = 440)
#' }
#'
#' @export
ucnv_semitone_to_hz <- function(semitone, ref_freq = NULL, ref_source = c("A4", "UEP83", "Praat"),
                            as_units = NULL) {
  ref_source <- match.arg(ref_source)

  # Determine reference frequency
  if (is.null(ref_freq)) {
    ref_freq <- switch(ref_source,
      "UEP83" = 110,    # A2 (A₁ in Helmholtz notation) - UEP 1983 standard
      "Praat" = 100,    # Praat's arbitrary reference for semitones
      "A4" = 440        # Concert pitch A4
    )
  }
  # Handle units if present
  has_units <- FALSE
  st_names <- names(semitone)  # Preserve names
  if (inherits(semitone, "units")) {
    has_units <- TRUE
    st_val <- as.numeric(semitone)
  } else {
    st_val <- as.numeric(semitone)
  }

  # Handle reference frequency units
  if (inherits(ref_freq, "units")) {
    ref_freq_hz <- as.numeric(units::set_units(ref_freq, Hz))
  } else {
    ref_freq_hz <- as.numeric(ref_freq)
  }

  # Validate input
  if (ref_freq_hz <= 0) {
    stop("Reference frequency must be positive")
  }

  # Apply inverse conversion formula
  # f = f_ref * 2^(ST/12)
  freq_hz <- ref_freq_hz * 2^(st_val / 12)

  # Restore names
  if (!is.null(st_names)) {
    names(freq_hz) <- st_names
  }

  # Determine if we should return units
  if (is.null(as_units)) {
    as_units <- has_units && requireNamespace("units", quietly = TRUE)
  }

  # Add units if requested
  if (as_units) {
    if (!requireNamespace("units", quietly = TRUE)) {
      warning("units package not available, returning plain numeric vector")
      return(freq_hz)
    }
    return(units::set_units(freq_hz, Hz))
  }

  freq_hz
}


# ============================================================================
# Helper functions for unit registration
# ============================================================================

#' @keywords internal
#' @noRd
.ensure_erb_unit <- function() {
  tryCatch({
    test <- units::as_units("ERB")
    return(invisible(NULL))
  }, error = function(e) {
    tryCatch({
      units::install_unit("ERB", "", "erb")
    }, error = function(e2) {
      invisible(NULL)
    })
  })
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.ensure_mel_unit <- function() {
  tryCatch({
    test <- units::as_units("mel")
    return(invisible(NULL))
  }, error = function(e) {
    # Mel unit doesn't exist, install it
    tryCatch({
      units::install_unit("mel", "")
    }, error = function(e2) {
      # Ignore if already installed by another process
      invisible(NULL)
    })
  })
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.ensure_semitone_unit <- function() {
  tryCatch({
    test <- units::as_units("semitone")
    return(invisible(NULL))
  }, error = function(e) {
    tryCatch({
      units::install_unit("semitone", "", "ST")
    }, error = function(e2) {
      invisible(NULL)
    })
  })
  invisible(NULL)
}
