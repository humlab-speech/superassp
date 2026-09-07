# Download and install ONNX Runtime native library

Download and install ONNX Runtime native library

## Usage

``` r
.ort_install(version = "1.24.3", path = NULL, gpu = FALSE)
```

## Arguments

- version:

  ONNX Runtime version to install. Default: "1.24.3".

- path:

  Installation directory. Default: platform-appropriate user cache.

- gpu:

  Logical. If TRUE, install GPU-enabled variant (Linux/Windows only).

## Value

Invisible path to installed library directory.
