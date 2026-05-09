# Enforce the export policy from the consistency refactor:
# user-exported functions must be one of trk_*, lst_*, ucnv_*, read_*, write_*,
# or one of the named class-method generics (dur, numRecs, rate, startTime, tracks).

test_that("only allowed function families are exported", {
  exported <- getNamespaceExports("superassp")
  allowed_re <- "^(trk_|lst_|ucnv_|read_|write_|dur$|numRecs$|rate$|startTime$|tracks$)"
  bad <- exported[!grepl(allowed_re, exported)]
  expect_equal(
    bad, character(0),
    info = paste("Disallowed exports:", paste(bad, collapse = ", "))
  )
})

test_that("read_/write_ exports are paired with their counterpart", {
  exported <- getNamespaceExports("superassp")
  reads  <- sub("^read_",  "", grep("^read_",  exported, value = TRUE))
  writes <- sub("^write_", "", grep("^write_", exported, value = TRUE))

  # `read_audio` has no `write_audio` counterpart by design (audio output is
  # via the source file or via av) — exempt it.
  reads_paired  <- setdiff(reads, c("audio"))
  writes_paired <- writes

  expect_setequal(reads_paired, writes_paired)
})
