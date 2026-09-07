# Setup S7 Method Dispatch for DSP Functions

Internal function called during .onLoad() to set up S7 dispatch for all
lst\_\* and trk\_\* functions. This enables them to work with AVAudio
objects while maintaining backward compatibility with file paths.

## Usage

``` r
.setup_s7_methods()
```

## Value

NULL (called for side effects)
