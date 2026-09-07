# Get external pointer from pladdrr R6 object

Extract the underlying C pointer from a pladdrr R6 object. Used
internally for passing objects to direct API functions.

## Usage

``` r
get_pladdrr_ptr(obj)
```

## Arguments

- obj:

  pladdrr R6 object (Sound, Pitch, Formant, etc.)

## Value

External pointer to the Praat C object

## Details

pladdrr uses R6 function factory pattern where all objects have a .xptr
field containing the C pointer.

## Examples

``` r
if (FALSE) { # \dontrun{
sound <- pladdrr::Sound("audio.wav")
ptr <- get_pladdrr_ptr(sound)
} # }
```
