# Ensure ONNX Runtime is available, installing automatically if needed

Called internally by any function that requires ONNX Runtime inference
(e.g. `trk_pitch_crepe`). On first call the function checks for a cached
installation; if none is found it downloads and installs the runtime
automatically with an informative message.

## Usage

``` r
ensure_onnx(version = "1.24.3", gpu = FALSE)
```

## Arguments

- version:

  ONNX Runtime version to install if not already present. Default:
  `"1.24.3"`.

- gpu:

  Logical. If TRUE, install GPU variant (Linux/Windows only).

## Value

Invisible TRUE.
