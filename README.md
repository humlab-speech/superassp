# An extended wrassp package

The idea is to make a package that has all the functionality of wrassp, and extend it with analyses made avaiable in Praat or MATLAB. The added functions should behave in a wrassp-like manner, and thereby be callable in a similar way in the `emuR` framwork.

The `praat_formant_burg` provides an illustration of how a Praat script that extracts formant values may be wrapped inside of an R function and produce a SSFF formant track file. 

## Details
By loading this package, you also get all the functions exported by the `wrassp` package into your namespace. This is achieved by the `superassp` package being *Depending*  the `wrassp` package (rather than *Importing*, which is usually the preferred way of creating depmendencies between R packages).

## Indications of performance of Praat and wrassp functions


```r
library(microbenchmark)
mb <- microbenchmark(
  praat_formant_burg=praat_formant_burg(
    file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE),
  forest=forest(
    file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE)
  )

```

which results in 

```
Unit: milliseconds
               expr       min        lq     mean    median        uq       max neval
 praat_formant_burg 113.57067 121.33921 128.5937 123.89103 133.10573 165.52178   100
             forest  25.92412  27.22489  28.1092  28.06385  28.65024  31.58972   100
             
```
Getting an SSFF file from a wrassp function rather than a wrapped Praat call (which also involves the parsing of a csv file) will normally be significantly faster.