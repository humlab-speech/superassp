test_that("lst_vq_vat returns plausible NAQ/QOQ/H1H2/HRF arrays", {
  skip_if_not_installed("voiceanalysis")
  wav <- system.file("samples/sustained/a1.wav", package = "superassp")
  skip_if(wav == "", "a1.wav not bundled")
  out <- lst_vq_vat(wav, toFile = FALSE)
  expect_named(out, c("gci_time", "NAQ", "QOQ", "H1H2", "HRF", "meta"))
  expect_gt(length(out$gci_time), 10)
  expect_equal(length(out$gci_time), length(out$NAQ))
  expect_equal(length(out$gci_time), length(out$QOQ))

  naq_nz <- out$NAQ[out$NAQ != 0]
  qoq_nz <- out$QOQ[out$QOQ != 0]
  if (length(naq_nz) > 5) expect_true(all(naq_nz > 0 & naq_nz < 1))
  if (length(qoq_nz) > 5) expect_true(all(qoq_nz >= 0 & qoq_nz <= 1))
})

test_that("lst_lf_vat_synthesis generates a pulse of finite samples", {
  skip_if_not_installed("voiceanalysis")
  r <- lst_lf_vat_synthesis(Rd = 1.0, F0 = 120, fs = 16000)
  expect_named(r, c("pulse", "Rd", "Ra", "Rk", "Rg", "F0", "fs", "EE"))
  expect_gt(length(r$pulse), 50)
  expect_true(all(is.finite(r$pulse)))
  # LF derivative pulse has a negative excitation peak
  expect_lt(min(r$pulse), 0)
})
