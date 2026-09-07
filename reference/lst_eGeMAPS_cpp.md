# Compute eGeMAPS Features (C++ Implementation)

Extracts Extended Geneva Minimalistic Acoustic Parameter Set (eGeMAPS)
features using OpenSMILE C++ library.

## Usage

``` r
lst_eGeMAPS_cpp(file, beginTime = 0, endTime = 0, verbose = FALSE)
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

Named list with 88 eGeMAPS features
