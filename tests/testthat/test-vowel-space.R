test_that("lst_vowel_space handles valid formant data", {
  skip_if_not_installed("superassp")

  # Create synthetic formant data from realistic vowel space
  set.seed(42)

  # Male vowel ranges: F1 < 900, F2 > 800, F1 > 250, F2 < 2450
  # Simulate points around /i/, /a/, /u/ with tight distributions to pass filters
  f1_vals <- c(
    rnorm(350, mean = 350, sd = 20),   # /i/ - low F1, tight
    rnorm(350, mean = 700, sd = 30),   # /a/ - mid-high F1
    rnorm(350, mean = 380, sd = 20)    # /u/ - low F1, tight
  )

  f2_vals <- c(
    rnorm(350, mean = 2000, sd = 80),  # /i/ - high F2
    rnorm(350, mean = 1200, sd = 70),  # /a/ - low-mid F2
    rnorm(350, mean = 1000, sd = 70)   # /u/ - low F2, tight
  )

  formant_data <- data.frame(F1 = f1_vals, F2 = f2_vals)

  result <- lst_vowel_space(formant_data, gender = 1, mode = "triangle")

  expect_is(result, "list")
  expect_true("vowel_space_ratio" %in% names(result))
  expect_true("centroids" %in% names(result))
  expect_true("n_frames" %in% names(result))
  expect_true(is.numeric(result$vowel_space_ratio))
  expect_true(result$vowel_space_ratio >= 0)
  # Should have at least some frames passing the frequency filters
  expect_true(result$n_frames >= 500)
})

test_that("lst_vowel_space respects gender parameter", {
  skip_if_not_installed("superassp")

  set.seed(42)

  # Create realistic male vowel space data
  f1_male <- c(
    rnorm(300, mean = 380, sd = 30),
    rnorm(300, mean = 650, sd = 40),
    rnorm(300, mean = 400, sd = 30),
    rnorm(100, mean = 500, sd = 100)
  )

  f2_male <- c(
    rnorm(300, mean = 2100, sd = 100),
    rnorm(300, mean = 1100, sd = 80),
    rnorm(300, mean = 1050, sd = 80),
    rnorm(100, mean = 1500, sd = 200)
  )

  formant_male <- data.frame(F1 = f1_male, F2 = f2_male)

  result_male <- lst_vowel_space(formant_male, gender = 1)

  # Male result should be valid
  expect_is(result_male, "list")
  expect_true(result_male$vowel_space_ratio >= 0)
})

test_that("lst_vowel_space returns zero ratio for insufficient data", {
  skip_if_not_installed("superassp")

  # Only 100 data points, but some might pass frequency filters
  # Just check that if we get zero ratio, the function handles it gracefully
  formant_data <- data.frame(
    F1 = rnorm(100, mean = 2000, sd = 50),  # Too high F1 for male (> 900)
    F2 = rnorm(100, mean = 500, sd = 100)   # Too low F2 for male (< 800)
  )

  result <- lst_vowel_space(formant_data, gender = 1)

  # Should return zero ratio due to insufficient data
  expect_equal(result$vowel_space_ratio, 0)
  # n_frames should be very small or zero
  expect_true(result$n_frames < 100)
})

test_that("lst_vowel_space respects mode parameter", {
  skip_if_not_installed("superassp")

  set.seed(42)

  # Create realistic male vowel space data
  f1_vals <- c(
    rnorm(300, mean = 380, sd = 30),
    rnorm(300, mean = 650, sd = 40),
    rnorm(300, mean = 400, sd = 30),
    rnorm(100, mean = 500, sd = 100)
  )

  f2_vals <- c(
    rnorm(300, mean = 2100, sd = 100),
    rnorm(300, mean = 1100, sd = 80),
    rnorm(300, mean = 1050, sd = 80),
    rnorm(100, mean = 1500, sd = 200)
  )

  formant_data <- data.frame(F1 = f1_vals, F2 = f2_vals)

  result_tri <- lst_vowel_space(formant_data, gender = 1, mode = "triangle")
  result_poly <- lst_vowel_space(formant_data, gender = 1, mode = "polygon")

  # Both should be valid
  expect_true(result_tri$vowel_space_ratio >= 0)
  expect_true(result_poly$vowel_space_ratio >= 0)
})

test_that("lst_vowel_space returns centroids matrix", {
  skip_if_not_installed("superassp")

  set.seed(42)

  # Create realistic male vowel space data with more samples per vowel
  f1_vals <- c(
    rnorm(400, mean = 380, sd = 20),
    rnorm(400, mean = 650, sd = 30),
    rnorm(400, mean = 400, sd = 20),
    rnorm(400, mean = 500, sd = 50)
  )

  f2_vals <- c(
    rnorm(400, mean = 2100, sd = 80),
    rnorm(400, mean = 1100, sd = 60),
    rnorm(400, mean = 1050, sd = 60),
    rnorm(400, mean = 1500, sd = 150)
  )

  formant_data <- data.frame(F1 = f1_vals, F2 = f2_vals)

  result <- lst_vowel_space(formant_data, gender = 1)

  expect_is(result$centroids, "matrix")
  expect_equal(ncol(result$centroids), 2)
  # Kmeans will produce centroids for all reference vowels (12 for male)
  expect_equal(nrow(result$centroids), 12)
})

test_that("lst_vowel_space requires valid gender parameter", {
  skip_if_not_installed("superassp")

  set.seed(42)
  f1_vals <- c(
    rnorm(300, mean = 380, sd = 30),
    rnorm(300, mean = 650, sd = 40),
    rnorm(300, mean = 400, sd = 30),
    rnorm(100, mean = 500, sd = 100)
  )
  f2_vals <- c(
    rnorm(300, mean = 2100, sd = 100),
    rnorm(300, mean = 1100, sd = 80),
    rnorm(300, mean = 1050, sd = 80),
    rnorm(100, mean = 1500, sd = 200)
  )
  formant_data <- data.frame(F1 = f1_vals, F2 = f2_vals)

  expect_error(lst_vowel_space(formant_data, gender = 3))
  expect_error(lst_vowel_space(formant_data, gender = -1))
})

test_that("lst_vowel_space requires valid mode parameter", {
  skip_if_not_installed("superassp")

  set.seed(42)
  f1_vals <- c(
    rnorm(300, mean = 380, sd = 30),
    rnorm(300, mean = 650, sd = 40),
    rnorm(300, mean = 400, sd = 30),
    rnorm(100, mean = 500, sd = 100)
  )
  f2_vals <- c(
    rnorm(300, mean = 2100, sd = 100),
    rnorm(300, mean = 1100, sd = 80),
    rnorm(300, mean = 1050, sd = 80),
    rnorm(100, mean = 1500, sd = 200)
  )
  formant_data <- data.frame(F1 = f1_vals, F2 = f2_vals)

  expect_error(lst_vowel_space(formant_data, mode = "invalid"))
})
