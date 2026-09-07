# Emit a consistent "Applying ()" progress message

Prints either: "Applying `fun()` to N recording(s)" or, when a time
window is active: "Applying `fun()` to a X.X second long portion of N
recording(s)"

## Usage

``` r
format_apply_msg(fun_name, n_files, beginTime = NULL, endTime = NULL)
```

## Arguments

- fun_name:

  Character. Public function name, e.g. "trk_pitch_rapt".

- n_files:

  Integer. Number of files being processed.

- beginTime:

  Numeric vector of begin times (seconds).

- endTime:

  Numeric vector of end times (seconds; 0 = end of file).
