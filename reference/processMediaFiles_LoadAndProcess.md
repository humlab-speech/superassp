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

  Additional parameters to pass to performAssp. Parallel processing is
  controlled through this argument: `parallel` (logical; default TRUE
  for more than one file) and `n_cores` (integer; default
  `detectCores() - 1`). Both are extracted from `...` before the
  remaining parameters are forwarded, so they never reach the underlying
  DSP routine.

## Value

List with:

- externalRes: Results from performAssp

- listOfFilesDF: Data frame with file processing information

- processed_native: Logical vector indicating which files were native

## Details

Supports parallel processing for batch operations on multi-core systems.
