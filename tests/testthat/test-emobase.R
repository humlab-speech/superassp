emobase_raw_arff <- function(wav) {
  smile <- system.file("opensmile", "bin", "SMILExtract", package = "superassp")
  conf <- system.file("opensmile", "config", "emobase", "emobase.conf",
                      package = "superassp")
  out <- tempfile(fileext = ".arff")
  on.exit(unlink(out))
  status <- system(sprintf("%s -C %s -I %s -O %s -instname off -l 0",
                           shQuote(smile), shQuote(conf), shQuote(wav),
                           shQuote(out)),
                   ignore.stdout = TRUE, ignore.stderr = TRUE)
  stopifnot(status == 0)
  lines <- readLines(out)
  nm <- sub("^@attribute\\s+(\\S+)\\s+.*$", "\\1",
            grep("^@attribute", lines, value = TRUE))
  data <- lines[(which(lines == "@data") + 1):length(lines)]
  fields <- strsplit(data[nzchar(trimws(data))][1], ",")[[1]]
  keep <- !nm %in% c("name", "class")
  stats::setNames(as.numeric(fields[keep]), nm[keep])
}

test_that("lst_emobase labels every feature with its own ARFF attribute", {
  wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(wav == "", "Test file not found")
  skip_if(system.file("opensmile", "bin", "SMILExtract", package = "superassp") == "",
          "SMILExtract not built")

  res <- lst_emobase(wav, verbose = FALSE)
  raw <- emobase_raw_arff(wav)
  expect_length(res, 988L)
  expect_identical(names(res), names(raw))
  expect_equal(unlist(res), raw)
  expect_true("pcm_intensity_sma_max" %in% names(res))
  expect_false("class" %in% names(res))
})

test_that("lst_emobase reads non-WAV and >16-bit input like the 16-bit WAV", {
  wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(wav == "", "Test file not found")
  ref <- lst_emobase(wav, verbose = FALSE)
  # the same recording stored losslessly in other containers / sample formats
  for (f in c("a2.flac", "a3.au", "a32b.wav")) {
    other <- system.file("samples", "sustained", f, package = "superassp")
    expect_equal(lst_emobase(other, verbose = FALSE), ref, info = f)
  }
  win <- lst_emobase(wav, beginTime = 1, endTime = 2, verbose = FALSE)
  expect_identical(names(win), names(ref))
  expect_false(isTRUE(all.equal(win, ref)))
})
