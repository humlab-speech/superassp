# The libassp wrappers run performAsspMemory(), which never writes a file
# itself; processMediaFiles_LoadAndProcess() must do it when toFile = TRUE.

assp_wrappers <- c(
  trk_rms = "rms", trk_zcr = "zcr", trk_acf = "acf", trk_afdiff = "dif",
  trk_cepstrum = "cep", trk_css_spectrum = "css", trk_dft_spectrum = "dft",
  trk_lps_spectrum = "lps", trk_formant_forest = "fms",
  trk_pitch_ksv = "fo", trk_pitch_mhs = "pit"
)

test_that("libassp wrappers write the in-memory result when toFile = TRUE", {
  wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(wav == "", "Test file not found")
  out <- withr::local_tempdir()

  for (fn in names(assp_wrappers)) {
    f <- get(fn, envir = asNamespace("superassp"))
    in_mem <- suppressWarnings(f(wav, toFile = FALSE, verbose = FALSE))
    n <- suppressWarnings(f(wav, toFile = TRUE, outputDirectory = out,
                            verbose = FALSE))
    path <- file.path(out, paste0("a1.", assp_wrappers[[fn]]))

    expect_equal(n, 1, info = fn)
    expect_true(file.exists(path), info = fn)
    # Files keep libassp's native track names; only the in-memory object is
    # renamed (e.g. "rms" -> "RMS[dB]"), so compare tracks by position.
    on_disk <- read_ssff(path)
    expect_equal(length(on_disk), length(in_mem), info = fn)
    for (i in seq_along(in_mem)) {
      expect_equal(unclass(on_disk[[i]]), unclass(in_mem[[i]]),
                   ignore_attr = TRUE, tolerance = 1e-6,
                   info = paste(fn, names(in_mem)[i]))
    }
  }
})

test_that("toFile = TRUE counts every file written in a batch", {
  wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(wav == "", "Test file not found")
  src <- withr::local_tempdir()
  files <- file.path(src, c("x.wav", "y.wav"))
  file.copy(wav, files)

  n <- trk_rms(files, toFile = TRUE, verbose = FALSE)
  expect_equal(n, 2)
  expect_true(all(file.exists(file.path(src, c("x.rms", "y.rms")))))
})
