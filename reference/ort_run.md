# Run ONNX Runtime inference

Run ONNX Runtime inference

## Usage

``` r
ort_run(session, inputs, shapes, output_names = NULL)
```

## Arguments

- session:

  An ORT session from
  [`ort_session()`](https://humlab-speech.github.io/superassp/reference/ort_session.md).

- inputs:

  Named list of input data. Each element is a numeric vector.

- shapes:

  List of integer vectors specifying the shape of each input. Must be in
  the same order as `inputs`.

- output_names:

  Character vector of output tensor names to fetch. NULL = fetch all
  outputs.

## Value

Named list of output tensors (numeric vectors with "shape" attribute).
