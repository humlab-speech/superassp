# Process media files with performAssp using av package (load-and-process pattern)

This internal function replaces the convert-then-process pattern with a
load-and-process pattern. Instead of converting media files to WAV on
disk, it loads them directly into memory using av, then processes them
with performAssp.

## Usage

``` r
processMediaFiles_LoadAndProcess(
  listOfFiles,
  beginTime,
  endTime,
  nativeFiletypes,
  fname,
  toFile = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- listOfFiles:

  Character vector of input file paths

- beginTime:

  Numeric vector of begin times (seconds)

- endTime:

  Numeric vector of end times (seconds, 0 = end of file)

- nativeFiletypes:

  Character vector of natively supported formats

- fname:

  Character name of performAssp function to call

- toFile:

  Logical whether to write output files

- verbose:

  Logical whether to show progress messages

- ...:

  Additional parameters to pass to performAssp

- parallel:

  Logical whether to use parallel processing (default TRUE for \>1
  files)

- n_cores:

  Integer number of cores to use (default: detectCores() - 1)

## Value

List with:

- externalRes: Results from performAssp

- listOfFilesDF: Data frame with file processing information

- processed_native: Logical vector indicating which files were native

## Details

Supports parallel processing for batch operations on multi-core systems.
