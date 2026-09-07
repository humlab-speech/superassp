# Get a cached path to a Hugging Face Hub model file, downloading if needed

Downloads `filename` from the Hugging Face Hub repo `repo_id` on first
use and caches it under the R user cache directory
(`tools::R_user_dir("superassp", "cache")`); subsequent calls return the
cached path without any network access.

## Usage

``` r
.hf_get_cached_model(repo_id, filename, subdir, revision, repo_type = "model")
```

## Arguments

- repo_id:

  Hugging Face Hub repo id, e.g. `"username/model-name"`.

- filename:

  Path to the file inside the repo, e.g. `"model.onnx"`.

- subdir:

  Cache subdirectory name for this model (keeps different models from
  colliding when their filenames share a basename), e.g. `"swift-f0"`.

- revision:

  Git revision (tag or commit SHA) to pin. Required, and deliberately
  has no default – pinning to `"main"` would let the model change under
  the package without a version bump, undermining DSP output
  faithfulness.

- repo_type:

  One of `"model"`, `"dataset"`, `"space"`. Default `"model"`.

## Value

Character path to the cached local file.
