test_that("lf_area_newton returns finite alpha/epsi for typical parameters", {
  res <- .vat_lf_area_newton(Tc = 0.008, fs = 16000, Tp = 0.005,
                            Te = 0.006, Ta = 0.0005, EE = 1.0)
  expect_true(is.finite(res$alpha))
  expect_true(is.finite(res$epsi))
  expect_gt(res$epsi, 0)
})

test_that("rd2r returns plausible R-parameters", {
  r <- .vat_rd2r(Rd = 1.0, EE = 1.0, F0 = 120)
  expect_gt(r$Ra, 0); expect_lt(r$Ra, 0.2)
  expect_gt(r$Rk, 0); expect_lt(r$Rk, 1)
  expect_gt(r$Rg, 0)
})

test_that("lf_cont generates a non-degenerate pulse of expected length", {
  r <- .vat_rd2r(1.0, 1.0, 120)
  pulse <- .vat_lf_cont(120, 16000, r$Ra, r$Rk, r$Rg, EE = 1.0)
  expect_gt(length(pulse), 50)
  expect_lt(length(pulse), 16000 / 120 + 5)
  # Should have a clear negative excitation peak
  expect_lt(min(pulse), -0.5)
  expect_true(all(is.finite(pulse)))
})

test_that("lf_cont varies meaningfully with Rd", {
  pulse_tense   <- {
    r <- .vat_rd2r(0.5, 1.0, 120); .vat_lf_cont(120, 16000, r$Ra, r$Rk, r$Rg, 1.0)
  }
  pulse_breathy <- {
    r <- .vat_rd2r(2.5, 1.0, 120); .vat_lf_cont(120, 16000, r$Ra, r$Rk, r$Rg, 1.0)
  }
  # The two pulses should be measurably different
  expect_true(length(pulse_tense) != length(pulse_breathy) ||
              abs(min(pulse_tense) - min(pulse_breathy)) > 1e-3)
})
