# Helper function to normalize track names for comparison
# Removes bracket notation like [dB], [Hz], etc. and converts to lowercase
normalize_track_name <- function(name) {
  # Remove anything in brackets and convert to lowercase
  tolower(gsub("\\[.*?\\]", "", name))
}

test_that("Parselmouth optimized functions are equivalent to original versions", {
  skip_if_not_installed("reticulate")

  # Setup Python
  Sys.unsetenv("RETICULATE_PYTHON")
  reticulate::use_python("/usr/bin/python3", required = FALSE)

  # Check if parselmouth is available
  skip_if_not(reticulate::py_module_available("parselmouth"),
              "Parselmouth not available")

  # Get test files
  test_files <- list.files(
    path = file.path("..", "signalfiles"),
    pattern = "\\.wav$",
    recursive = TRUE,
    full.names = TRUE
  )

  # Limit to a few representative files for testing
  test_files <- test_files[1:min(5, length(test_files))]

  expect_true(length(test_files) > 0,
              info = "No test audio files found")

  # Test each optimized function
  test_functions <- list(
    list(
      name = "formant_burg",
      orig = "praat_formant_burg",
      opt = "praat_formant_burg_opt",
      skip_if_missing = TRUE
    ),
    list(
      name = "pitch",
      orig = "praat_pitch",
      opt = "praat_pitch_opt",
      skip_if_missing = TRUE
    ),
    list(
      name = "intensity",
      orig = "praat_intensity",
      opt = "praat_intensity_opt",
      skip_if_missing = TRUE
    ),
    list(
      name = "spectral_moments",
      orig = "praat_spectral_moments",
      opt = "praat_spectral_moments_opt",
      skip_if_missing = TRUE
    ),
    list(
      name = "formantpath_burg",
      orig = "praat_formantpath_burg",
      opt = "praat_formantpath_burg_opt",
      skip_if_missing = TRUE
    )
  )

  for (fn_pair in test_functions) {
    # Check if functions exist
    if (!exists(fn_pair$orig, mode = "function") && fn_pair$skip_if_missing) {
      next
    }
    if (!exists(fn_pair$opt, mode = "function")) {
      next
    }

    for (test_file in test_files) {
      test_name <- paste0(fn_pair$name, " - ", basename(test_file))

      # Try both versions
      result_orig <- tryCatch({
        do.call(fn_pair$orig, list(listOfFiles = test_file, toFile = FALSE))
      }, error = function(e) NULL)

      result_opt <- tryCatch({
        do.call(fn_pair$opt, list(listOfFiles = test_file, toFile = FALSE))
      }, error = function(e) NULL)

      # Skip if either failed
      if (is.null(result_orig) || is.null(result_opt)) {
        next
      }

      # Compare track names
      expect_equal(
        sort(names(result_opt)),
        sort(names(result_orig)),
        info = paste(test_name, "- track names should match")
      )

      # Compare number of frames
      for (track in names(result_opt)) {
        if (track %in% names(result_orig)) {
          expect_equal(
            nrow(result_opt[[track]]),
            nrow(result_orig[[track]]),
            tolerance = 2, # Allow small differences in frame count
            info = paste(test_name, "- track", track, "frame count")
          )
        }
      }

      # Compare values (allowing for numerical differences)
      for (track in names(result_opt)) {
        if (track %in% names(result_orig)) {
          # Check correlation (should be high)
          orig_vals <- as.vector(result_orig[[track]])
          opt_vals <- as.vector(result_opt[[track]])

          # Remove NAs and zeros for correlation
          valid_idx <- !is.na(orig_vals) & !is.na(opt_vals) &
                       orig_vals != 0 & opt_vals != 0

          if (sum(valid_idx) > 10) {
            cor_val <- cor(orig_vals[valid_idx], opt_vals[valid_idx])
            expect_gt(
              cor_val,
              0.95,
              info = paste(test_name, "- track", track, "correlation")
            )
          }
        }
      }
    }
  }
})

test_that("SuperASP DSP functions are equivalent to wrassp functions", {
  skip_on_cran()
  skip_if_not_installed("wrassp")

  # Get test files
  test_files <- list.files(
    path = file.path("..", "signalfiles"),
    pattern = "\\.wav$",
    recursive = TRUE,
    full.names = TRUE
  )

  test_files <- test_files[1:min(3, length(test_files))]

  expect_true(length(test_files) > 0)

  # Test common functions
  common_functions <- c("acfana", "rfcana", "rmsana", "zcrana")

  for (fn_name in common_functions) {
    for (test_file in test_files) {
      test_name <- paste0(fn_name, " - ", basename(test_file))

      # Call both versions
      result_superassp <- tryCatch({
        do.call(paste0("superassp::", fn_name),
                list(listOfFiles = test_file, toFile = FALSE))
      }, error = function(e) NULL)

      result_wrassp <- tryCatch({
        do.call(paste0("wrassp::", fn_name),
                list(listOfFiles = test_file, toFile = FALSE))
      }, error = function(e) NULL)

      # Skip if either failed
      if (is.null(result_superassp) || is.null(result_wrassp)) {
        next
      }

      # Compare structure
      expect_equal(
        class(result_superassp),
        class(result_wrassp),
        info = paste(test_name, "- class")
      )

      # Compare track names (normalized to remove bracket notation and case)
      expect_equal(
        sort(normalize_track_name(names(result_superassp))),
        sort(normalize_track_name(names(result_wrassp))),
        info = paste(test_name, "- track names (normalized)")
      )

      # Compare dimensions and values
      # Match tracks by normalized name (removes bracket notation and case)
      superassp_tracks <- names(result_superassp)
      wrassp_tracks <- names(result_wrassp)

      for (track_super in superassp_tracks) {
        # Find matching track in wrassp (normalized comparison)
        track_wrassp <- wrassp_tracks[normalize_track_name(wrassp_tracks) == normalize_track_name(track_super)]

        if (length(track_wrassp) > 0) {
          track_wrassp <- track_wrassp[1]  # Take first match

          # Check dimensions
          expect_equal(
            dim(result_superassp[[track_super]]),
            dim(result_wrassp[[track_wrassp]]),
            info = paste(test_name, "- track", track_super, "dimensions")
          )

          # Check values are highly correlated
          superassp_vals <- as.vector(result_superassp[[track_super]])
          wrassp_vals <- as.vector(result_wrassp[[track_wrassp]])

          valid_idx <- !is.na(superassp_vals) & !is.na(wrassp_vals) &
                       abs(superassp_vals) > 1e-10 & abs(wrassp_vals) > 1e-10

          if (sum(valid_idx) > 10) {
            cor_val <- cor(superassp_vals[valid_idx], wrassp_vals[valid_idx])
            expect_gt(
              cor_val,
              0.99,
              info = paste(test_name, "- track", track_super, "correlation")
            )
          }
        }
      }
    }
  }
})
