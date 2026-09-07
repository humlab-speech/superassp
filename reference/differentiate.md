# Derivation of SSFF track objects

This function takes an SSFF object or file and computes the `order`th
derivative of the tracks in it. The user may also specify a lag of
differentiation. If `lag=1`, then ordinary differences between
consecutive values are computed. If `lag=2`, then the difference between
the 1st and 3rd value will be returned, and so on. The user may specify
an order of differentiation too, and thereby cause the differentiation
to be conducted in multiple iterations.

## Usage

``` r
differentiate(
  inSSFF,
  order = 1,
  onlyTracks = NULL,
  padLeft = TRUE,
  toFile = TRUE,
  explicitExt = NULL,
  overwrite = FALSE
)
```

## Arguments

- inSSFF:

  The SSFF object, or a full path to a file that contains an SSFF object
  and may be read as such by
  [read.AsspDataObj](https://humlab-speech.github.io/superassp/reference/read.AsspDataObj.md).

- order:

  The number of iterations in which the vector will be differentiated.
  The first order differentiation gives the size of changes in
  consecutive values (with an indicated lag). The second order
  differentiation gives the rate of change, and so on.

- onlyTracks:

  Only differentiate certain tracks, and leave the others as is.
  Defaults to process all tracks.

- padLeft:

  Should initial zeros be inserted into the vector from the left?

- toFile:

  Should the resulting SSFF object be written to file, or returned?

- explicitExt:

  By default, a character "d" will be prepended to the file name suffix
  when writing the output to file. The user can also specify an explicit
  extension which will be used instead.

- overwrite:

  Should an existing file be overwritten when writing the output?

## Value

The function will return an SSFF object if `toFile` is `TRUE`.
Otherwise, nothing is returned.

## Details

Differentiation always results in loss of data, and the user may specify
how to align the differentiation output. Initial zero padding values
will be inserted so that the vector length of the input and the output
will always be the same. If `padLeft=TRUE` (the default) the initial
zero values will be inserted into the tracks so that the differentiation
result aligns in time with the occurrence of a value change. That is, if
`lag=2` a value in the output vector indicates that at that place in the
input vector a change has happened of the indicated size compared to the
value two positions back in the vector. If `padLeft=FALSE` and `lag=2`
then the value indicates the change that will have occurred when looking
two positions forward in the vector. This is likely an unusual use case,
and therefore not the default behavior.

Padding the signal with zeros is performed after all iterations of
differentiation have been performed completely, and the padding zeros
will therefore never be differentiated themselves.
