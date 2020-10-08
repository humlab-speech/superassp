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

# Steps to implement a new Praat function

1. Indentify what the output of the function will be
    a. A signal track (or tracks) that follows the original sound wave
    b. A value (or a limited list of values) that summarises the acoustic properties of a wav file, and can therefore not sensibly be shown alongside the sound wave.
2. Implement the core analysis in a Praat script file, and place it in `inst/praat`.
    a. In the case where track(s) that follow the sound wave file are returned, the Praat function should write the output to a CSV table file and return the name of that table. The Praat script should also take the desired output table file name (inkluding full path) as an argument. Please refere to `praat/formant_burg.praat` for some example code that computes formants and bandwidths for them for a (possibly windowed) sound file and writes them to a table.
3. Make a copy of the suitable template function, rename it (please keep the praat_ prefix for clarity) and make modifications to the code to suit the new track computed by Praat. You will need to think about what the tracks should be called in the SSFF file and document your choice.
    a. For a function that computes a sound wave following signal track (or tracks), use the code of `praat_formant_burg` as a template. Please refer to a suitable function in wrassp for inspiration on what to call sets of tracks. (The `praat_formant_burg` outputs and "fm" and "bw" set, for formant frequencies and formant bandwidths respectivelly)
    b. For single value (or list of values) output, there is currently no template function implemented, but please note that the `tjm.praat::wrap_praat_script()` has an option to return the "Info window" of Praat, which opens up lots of possibilities.
4. There are many moving parts to this whole package, so make sure to contruct a test file and a test suit for the new function to make sure that it works. 