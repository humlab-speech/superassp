# read.AsspDataObj from a signal/parameter file

read.AsspDataObj creates an object of class dobj from a signal or
parameter file readable by the ASSP Library (WAVE, SSFF, AU, ...)

## Usage

``` r
read.AsspDataObj(fname, begin = 0, end = 0, samples = FALSE)
```

## Arguments

- fname:

  filename of the signal or parameter file

- begin:

  begin time (default is in seconds) of segment to retrieve

- end:

  end time (default is in seconds) of segment to retrieve

- samples:

  (BOOL) if set to false seconds values of begin/end are sample numbers

## Value

list object containing file data

## Author

Lasse Bombien
