# Generic OpenSMILE Feature Extraction using C++

Internal function that extracts OpenSMILE features using the C++
implementation

## Usage

``` r
opensmile_extract_generic(
  file,
  config_name,
  config_dir,
  feature_set_name = "features",
  beginTime = 0,
  endTime = 0,
  verbose = FALSE
)
```

## Arguments

- file:

  Path to audio file

- config_name:

  Name of the config file (without extension)

- config_dir:

  Directory path relative to inst/opensmile/config

- feature_set_name:

  Name for verbose output

- beginTime:

  Start time in seconds

- endTime:

  End time in seconds

- verbose:

  Print processing information

## Value

Named list with acoustic features
