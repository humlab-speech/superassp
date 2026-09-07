# Create an ONNX Runtime inference session

Create an ONNX Runtime inference session

## Usage

``` r
ort_session(model_path, num_threads = 0L)
```

## Arguments

- model_path:

  Path to an ONNX model file (.onnx).

- num_threads:

  Number of intra-op threads. 0 = auto (ORT default).

## Value

An external pointer to the ORT session (class "ort_session").

## Details

The session holds the loaded model and is reusable across multiple
inference calls. It is automatically released when garbage collected.

## Examples

``` r
if (FALSE) { # \dontrun{
sess <- ort_session("model.onnx")
result <- ort_run(sess, list(input = rnorm(1024)), list(c(1L, 1024L)))
} # }
```
