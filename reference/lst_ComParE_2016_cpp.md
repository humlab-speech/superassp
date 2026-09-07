# Compute ComParE 2016 Features (C++ Implementation)

Extracts ComParE 2016 features using OpenSMILE C++ library.

## Usage

``` r
lst_ComParE_2016_cpp(file, beginTime = 0, endTime = 0, verbose = FALSE)
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

Named list with 6373 ComParE features
