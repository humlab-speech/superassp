# Convert input media files (Rcpp-optimized)

Convert input media files (Rcpp-optimized)

Convert input media files (Rcpp-optimized)

## Usage

``` r
convertInputMediaFiles(
  listOfFiles,
  beginTime,
  endTime,
  windowShift = 5,
  nativeFiletypes,
  preferedFiletype,
  knownLossless,
  funName,
  keepConverted,
  verbose
)

convertInputMediaFiles(
  listOfFiles,
  beginTime,
  endTime,
  windowShift = 5,
  nativeFiletypes,
  preferedFiletype,
  knownLossless,
  funName,
  keepConverted,
  verbose
)
```

## Arguments

- listOfFiles:

  Character vector of input file paths

- beginTime:

  Numeric vector of begin times

- endTime:

  Numeric vector of end times

- windowShift:

  Numeric window shift in milliseconds

- nativeFiletypes:

  Character vector of natively supported formats

- preferedFiletype:

  Character string of preferred conversion format

- knownLossless:

  Character vector of known lossless formats

- funName:

  Character name of calling function

- keepConverted:

  Logical whether to keep converted files

- verbose:

  Logical whether to show progress messages

## Value

List with processed file information and files to clean

List with processed file information and files to clean
