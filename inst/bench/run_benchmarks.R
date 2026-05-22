#!/usr/bin/env Rscript
# Performance benchmark harness for superassp v2.7.1.
#
# Times the user-facing trk_*/lst_* wrappers against a fixed test corpus.
# Designed to be run before and after a perf change so the delta is visible.
#
# Usage:
#   Rscript inst/bench/run_benchmarks.R [--repeats=N] [--out=path.rds]
#
# The output .rds contains a data.frame with columns:
#   fn, n_files, repeats, elapsed_s, user_s, system_s
#
# To compare baseline vs post:
#   pre  <- readRDS("tests/testthat/_snaps/perf/baseline.rds")
#   post <- readRDS("tests/testthat/_snaps/perf/post.rds")
#   merge(pre, post, by = "fn", suffixes = c("_pre", "_post"))

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}

repeats <- as.integer(parse_arg("repeats", 3L))
out     <- parse_arg("out", "bench-output.rds")

suppressMessages({
  library(superassp)
})

wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
if (!nzchar(wav)) stop("Sample WAV not available — install superassp first.")

corpus_small  <- rep(wav, 1L)
corpus_medium <- rep(wav, 4L)
corpus_large  <- rep(wav, 16L)

bench_one <- function(label, fn, listOfFiles) {
  gc()
  t <- system.time(for (i in seq_len(repeats)) fn(listOfFiles))
  data.frame(
    fn        = label,
    n_files   = length(listOfFiles),
    repeats   = repeats,
    elapsed_s = t[["elapsed"]],
    user_s    = t[["user.self"]],
    system_s  = t[["sys.self"]],
    stringsAsFactors = FALSE
  )
}

results <- list()
add <- function(...) results[[length(results) + 1L]] <<- bench_one(...)

cat("Benchmarking with repeats =", repeats, "...\n")

# Add specific wrappers worth tracking. Skip any that error (missing deps).
try_add <- function(label, fn, listOfFiles) {
  tryCatch(add(label, fn, listOfFiles), error = function(e) {
    cat(sprintf("  SKIP %s — %s\n", label, conditionMessage(e)))
  })
}

try_add("trk_pitch_rapt (1)",   function(x) trk_pitch_rapt(x, toFile = FALSE, verbose = FALSE), corpus_small)
try_add("trk_pitch_rapt (4)",   function(x) trk_pitch_rapt(x, toFile = FALSE, verbose = FALSE), corpus_medium)
try_add("trk_pitch_rapt (16)",  function(x) trk_pitch_rapt(x, toFile = FALSE, verbose = FALSE), corpus_large)
try_add("trk_rms (16)",         function(x) trk_rms(x, toFile = FALSE, verbose = FALSE),        corpus_large)
try_add("read_audio (16)",      function(x) lapply(x, read_audio),                              corpus_large)

df <- do.call(rbind, results)
print(df, row.names = FALSE)

saveRDS(df, file = out)
cat("\nSaved benchmark results to:", out, "\n")
