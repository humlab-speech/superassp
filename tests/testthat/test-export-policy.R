# Enforce the export policy from the consistency refactor:
# user-exported functions must be one of trk_*, lst_*, ucnv_*, read_*, write_*,
# the new snake_case generics (sample_rate, n_records, etc.), or the deprecated
# legacy generics (dur, numRecs, rate, startTime, tracks) kept for 2.8.x compat.

test_that("only allowed function families are exported", {
  exported <- getNamespaceExports("superassp")
  allowed_re <- paste0(
    "^(trk_|lst_|ucnv_|read_|write_",
    "|sample_rate$|n_records$|signal_duration$|start_time$|track_names$",
    "|file_path$|track_formats$",
    "|dur$|numRecs$|rate$|startTime$|tracks$",
    ")"
  )
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
