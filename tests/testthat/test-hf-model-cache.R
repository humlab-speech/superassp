test_that("revision must be a pinned tag/commit, not left to default", {
  expect_error(
    superassp:::.hf_get_cached_model("someuser/some-model", "model.onnx", "unit-test"),
    class = "rlang_error"
  )
  expect_error(
    superassp:::.hf_get_cached_model("someuser/some-model", "model.onnx", "unit-test", revision = ""),
    class = "rlang_error"
  )
})

test_that("a cached file is returned without touching the network", {
  tmp <- withr::local_tempdir()
  testthat::local_mocked_bindings(
    R_user_dir = function(...) tmp,
    .package = "tools"
  )

  cache_dir <- file.path(tmp, "onnx", "unit-test")
  dir.create(cache_dir, recursive = TRUE)
  writeLines("fake model bytes", file.path(cache_dir, "model.onnx"))

  path <- superassp:::.hf_get_cached_model(
    "someuser/some-model", "model.onnx", "unit-test",
    revision = "v1"
  )
  expect_equal(path, file.path(cache_dir, "model.onnx"))
})

test_that("a cache miss without huggingfaceR installed aborts clearly", {
  skip_if(requireNamespace("huggingfaceR", quietly = TRUE),
          "huggingfaceR is installed; this test targets the not-installed path")

  tmp <- withr::local_tempdir()
  testthat::local_mocked_bindings(
    R_user_dir = function(...) tmp,
    .package = "tools"
  )

  expect_error(
    superassp:::.hf_get_cached_model(
      "someuser/some-model", "model.onnx", "unit-test",
      revision = "v1"
    ),
    "huggingfaceR"
  )
})

test_that("a cache miss with huggingfaceR installed downloads into the cache dir", {
  skip_if_not_installed("huggingfaceR")

  tmp <- withr::local_tempdir()
  testthat::local_mocked_bindings(
    R_user_dir = function(...) tmp,
    .package = "tools"
  )

  cache_dir <- file.path(tmp, "onnx", "unit-test")
  called_with <- NULL
  testthat::local_mocked_bindings(
    hf_hub_download = function(repo_id, filename, repo_type, revision, dest, ...) {
      called_with <<- list(repo_id = repo_id, filename = filename,
                            repo_type = repo_type, revision = revision, dest = dest)
      dir.create(dest, recursive = TRUE, showWarnings = FALSE)
      out <- file.path(dest, basename(filename))
      writeLines("fake model bytes", out)
      out
    },
    .package = "huggingfaceR"
  )

  path <- superassp:::.hf_get_cached_model(
    "someuser/some-model", "model.onnx", "unit-test",
    revision = "v1"
  )

  expect_equal(path, file.path(cache_dir, "model.onnx"))
  expect_true(file.exists(path))
  expect_equal(called_with$repo_id, "someuser/some-model")
  expect_equal(called_with$revision, "v1")
})
