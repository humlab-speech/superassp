library(testthat)

# CI diagnostic: windows-latest test runs have been crashing/hanging (R CMD
# check reports a bare ERROR after several minutes with no testthat output
# at all -- the OS-buffered stdout is lost when the R process dies, so the
# check log never shows which test was running). This reporter prints and
# force-flushes "Start test: ..." before each test, so a crash mid-run still
# leaves a trail in 00check.log. Remove once the windows crash is
# root-caused and fixed; see 00check.log from a failing windows-latest run.
FlushingLocationReporter <- R6::R6Class("FlushingLocationReporter",
  inherit = testthat::LocationReporter,
  public = list(
    start_test = function(context, test) {
      super$start_test(context, test)
      flush(stdout())
    }
  )
)
test_check("superassp",
  reporter = testthat::MultiReporter$new(list(
    testthat::CheckReporter$new(),
    FlushingLocationReporter$new()
  ))
)

logger::log_threshold(logger::WARN)