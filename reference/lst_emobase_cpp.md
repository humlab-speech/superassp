# Compute emobase Features via SMILExtract (C++ Implementation)

Extracts emobase features using the SMILExtract command-line tool. This
approach is used because emobase's frameMode=full is incompatible with
the external audio source architecture used by other feature sets.

## Usage

``` r
lst_emobase_cpp(file, beginTime = 0, endTime = 0, verbose = FALSE)
```

## Arguments

- file:

  Path to audio file

- beginTime:

  Start time in seconds (default: 0)

- endTime:

  End time in seconds (default: 0 = end of file)

- verbose:

  Print processing information (default: FALSE)

## Value

Named list with 988 emobase features
