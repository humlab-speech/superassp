#' Vowel space analysis (F1×F2 area ratio)
#'
#' Computes vowel space area from formant data using k-means clustering and
#' reference vowel comparison. The vowel space ratio is computed as the ratio
#' of the actual vowel space convex hull area to a reference vowel space area
#' for the same gender.
#'
#' @param formant_data Data frame or matrix with columns F1, F2 (or first two columns).
#'   F1 and F2 should be in Hz.
#' @param gender Speaker gender for reference vowel selection:
#'   0 = female, 1 = male (default), 2 = child
#' @param mode Vowel space computation mode:
#'   "triangle" = use 3 corner vowels (default, more robust)
#'   "polygon" = use 4 vowels (may be less stable)
#' @param scaling Apply Bark frequency scaling (default: FALSE)
#' @param plot_formants Show scatter plot of formant points (default: FALSE)
#'
#' @return
#' List with elements:
#' - `vowel_space_ratio`: Numeric, ratio of actual to reference vowel space area (0–2+ range)
#' - `centroids`: Matrix of estimated vowel formant centroids from k-means clustering
#' - `n_frames`: Number of frames used in computation
#'
#' @details
#' **Algorithm**:
#' 1. Filter formant data by gender-specific frequency ranges to identify candidate vowels
#' 2. Run k-means clustering using reference vowel centroids as initialization
#' 3. Match estimated centroids to reference vowels using Mahalanobis distance
#' 4. Compute convex hull of matched vowel points (actual vowel space)
#' 5. Compute convex hull of reference vowels
#' 6. Ratio = actual area / reference area
#'
#' **Interpretation**:
#' - Ratio = 1.0: perfect agreement with reference vowel space
#' - Ratio < 1.0: speaker has reduced vowel space (dysarthria, age, etc.)
#' - Ratio > 1.0: speaker has expanded vowel space (hyperarticulation, singing)
#' - NaN/error: insufficient vowel tokens (<1000 frames) in valid range
#'
#' **Reference vowel sets**:
#' - Triangle: /i/, /a/, /u/ (corner vowels)
#' - Polygon: /i/, /ε/, /a/, /u/ (4-vowel system)
#'
#' **Minimum frame requirement**: 1000 frames must be within gender-specific frequency bounds,
#' or function returns 0 ratio (insufficient data).
#'
#' @examples
#' \dontrun{
#' # Compute vowel space from file batch
#' files <- c("speaker1.wav", "speaker2.wav")
#' formants <- trk_formant(files, toFile = FALSE)
#'
#' # Extract F1-F2 and compute vowel space
#' formant_df <- data.frame(F1 = formants$fm1, F2 = formants$fm2)
#' vs <- lst_vowel_space(formant_df, gender = 1)
#' cat("Vowel space ratio:", vs$vowel_space_ratio, "\n")
#'
#' # For female speaker
#' vs_female <- lst_vowel_space(formant_df, gender = 0)
#'
#' # Using Bark scaling
#' vs_bark <- lst_vowel_space(formant_df, gender = 1, scaling = TRUE)
#' }
#'
#' @export
lst_vowel_space <- function(formant_data,
                            gender = 1,
                            mode = "triangle",
                            scaling = FALSE,
                            plot_formants = FALSE) {

  # Input validation
  if (!gender %in% c(0, 1, 2)) {
    cli::cli_abort("gender must be 0 (female), 1 (male), or 2 (child).")
  }
  if (!mode %in% c("triangle", "polygon")) {
    cli::cli_abort("mode must be 'triangle' or 'polygon'.")
  }

  # Convert to matrix and extract F1-F2
  formant_data <- as.matrix(formant_data)
  if (ncol(formant_data) < 2L) {
    cli::cli_abort("formant_data must have at least two columns (F1 and F2).")
  }

  formant_data <- formant_data[, 1:2, drop = FALSE]

  # Reference vowel formants (Peterson & Barney 1952 acoustic vowel space)
  # Vowel order: /i/, /ɪ/, /ɛ/, /æ/, /a/, /ɑ/, /ɔ/, /o/, /ʊ/, /u/, /ʌ/, /ə/
  men_f1 <- c(342, 427, 476, 580, 588, 768, 652, 497, 469, 378, 623, 474)
  men_f2 <- c(2322, 2034, 2089, 1799, 1952, 1333, 997, 910, 1122, 997, 1200, 1379)

  fem_f1 <- c(437, 483, 536, 731, 669, 936, 781, 555, 519, 459, 753, 523)
  fem_f2 <- c(2761, 2365, 2530, 2058, 2349, 1551, 1136, 1035, 1225, 1105, 1426, 1588)

  child_f1 <- c(452, 511, 564, 749, 717, 1002, 803, 597, 568, 494, 794, 586)
  child_f2 <- c(3081, 2552, 2656, 2267, 2501, 1688, 1210, 1137, 1490, 1345, 1546, 1719)

  # Bark frequency scaling (perceptual frequency warping)
  bark_scaling <- function(x) {
    13 * atan(0.00076 * x) + 3.5 * atan((x / 7500)^2)
  }

  # Select target vowels based on mode
  # Triangle: /i/ (1), /a/ (5 = /ɑ/), /u/ (10 = /u/)
  # Polygon: /i/ (1), /ɛ/ (3), /a/ (5), /u/ (10)
  target_vowels <- if (mode == "polygon") c(1, 3, 5, 10) else c(1, 5, 10)

  # Select gender-specific parameters
  if (gender == 1) {
    # Male
    init_formants <- cbind(men_f1, men_f2)
    # Frequency ranges for male speakers
    keep <- (formant_data[, 1] < 900 & formant_data[, 2] > 800 &
             formant_data[, 1] > 250 & formant_data[, 2] < 2450)
  } else if (gender == 0) {
    # Female
    init_formants <- cbind(fem_f1, fem_f2)
    # Frequency ranges for female speakers (higher F1/F2)
    keep <- (formant_data[, 1] < 1000 & formant_data[, 2] > 1000 &
             formant_data[, 1] > 350 & formant_data[, 2] < 2800)
  } else {
    # Child
    init_formants <- cbind(child_f1, child_f2)
    # Frequency ranges for children (even higher)
    keep <- (formant_data[, 1] < 1000 & formant_data[, 2] > 1000 &
             formant_data[, 1] > 350 & formant_data[, 2] < 2800)
  }

  formant_data_filtered <- formant_data[keep, 1:2, drop = FALSE]

  # Check minimum sample size (1000 frames)
  if (nrow(formant_data_filtered) < 1000L) {
    return(list(
      vowel_space_ratio = 0,
      centroids = matrix(numeric(), ncol = 2L),
      n_frames = nrow(formant_data_filtered)
    ))
  }

  # Apply Bark scaling if requested
  if (scaling) {
    formant_data_filtered <- bark_scaling(formant_data_filtered)
    init_formants <- bark_scaling(init_formants)
  }

  # K-means clustering (with robustness to empty clusters)
  km <- tryCatch(
    stats::kmeans(formant_data_filtered,
                  centers = init_formants,
                  iter.max = 100,
                  nstart = 1,
                  algorithm = "Hartigan-Wong"),
    error = function(e) {
      # If standard kmeans fails, try Lloyd's algorithm
      stats::kmeans(formant_data_filtered,
                    centers = init_formants,
                    iter.max = 100,
                    nstart = 1,
                    algorithm = "Lloyd")
    }
  )
  centroids <- km$centers

  # Mahalanobis distance matching between estimated and reference vowels
  sds <- apply(rbind(init_formants, centroids), 2, stats::sd)
  sds[sds == 0] <- 1

  pairwise <- outer(
    seq_len(nrow(init_formants)),
    seq_len(nrow(centroids)),
    Vectorize(function(i, j) {
      sqrt(sum(((init_formants[i, ] - centroids[j, ]) / sds)^2))
    })
  )

  # Match target vowel reference to estimated clusters
  idx <- apply(pairwise[target_vowels, , drop = FALSE], 1, which.min)
  conv_hull_points <- centroids[idx, , drop = FALSE]
  vowel_space_points <- init_formants[target_vowels, , drop = FALSE]

  # Compute convex hulls
  cur_space_idx <- grDevices::chull(conv_hull_points[, 1], conv_hull_points[, 2])
  ref_space_idx <- grDevices::chull(vowel_space_points[, 1], vowel_space_points[, 2])

  # Polygon area computation (shoelace formula)
  polygon_area <- function(pts, ord) {
    pts_ordered <- pts[c(ord, ord[1]), , drop = FALSE]
    abs(sum(pts_ordered[-1, 1] * pts_ordered[-nrow(pts_ordered), 2] -
            pts_ordered[-nrow(pts_ordered), 1] * pts_ordered[-1, 2])) / 2
  }

  # Compute ratio
  cur_area <- polygon_area(conv_hull_points, cur_space_idx)
  ref_area <- polygon_area(vowel_space_points, ref_space_idx)
  ratio <- cur_area / ref_area

  if (plot_formants) {
    plot(formant_data_filtered[, 1], formant_data_filtered[, 2],
         pch = 16, cex = 0.3,
         xlab = "F1 (Hz)", ylab = "F2 (Hz)",
         main = "Vowel formants with cluster centroids")
    points(centroids[, 1], centroids[, 2], col = "red", pch = 1, cex = 1.2)
  }

  list(
    vowel_space_ratio = ratio,
    centroids = centroids,
    n_frames = nrow(formant_data_filtered)
  )
}
