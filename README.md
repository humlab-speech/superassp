# An extended wrassp package

The idea is to make a package that has all the functionality of wrassp, and extend it with analyses made avaiable in Praat or MATLAB. The added functions should behave in a wrassp-like manner, and thereby be callable in a similar way in the `emuR` framwork.

The `praat_formant_burg` provides an illustration of how a Praat script that extracts formant values may be wrapped inside of an R function and produce a SSFF formant track file. 

## Details
By loading this package, you also get all the functions exported by the `wrassp` package into your namespace. This is achieved by the `superassp` package being *Depending*  the `wrassp` package (rather than *Importing*, which is usually the preferred way of creating depmendencies between R packages).

## Installation

The package requires the Praat program to be installed in the user's PATH (or in '/Applications' on Mac OS).

Then simply install the package using
```r
install.packages("devtools") # If not installed already
devtools::install_github("humlab-speech/superassp",dependencies = "Imports")
```

## Indications of performance of Praat and wrassp functions


```r
library(microbenchmark)
microbenchmark(
  "wrassp::forest"=wrassp::forest(
    file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE),
  praat_formant_burg=praat_formant_burg(
    file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE),

  praat_sauce=praat_sauce(
    file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE,pitchTracking = FALSE,spectralMeasures = FALSE)
)
```

which results in 

```
Unit: milliseconds
               expr       min        lq      mean    median        uq       max neval
     wrassp::forest  24.73287  26.97004  27.64405  27.70719  28.21722  30.35287   100
 praat_formant_burg 274.22872 287.31310 300.76562 294.81349 305.57010 389.37702   100
        praat_sauce 607.88763 642.12963 680.40282 686.59430 703.72266 834.41275   100
             
```
Getting an SSFF file from a `wrassp` function rather than a wrapped Praat call (which also involves the parsing of a csv file) will normally be significantly faster. Also, even it is adviced that even though the functions `praat_sauce` does compute formant values as well, it is not really efficiently implemented and is really mostly there for correction of harmonic amplitudes. And, an additional factor to consider is that the formant tracks will be stored by the `praat_sauce` function in a field in the same file as all the other tracks computed by the function, which will likelly result in a performance issue when working with the tracks.

So, if you need only formant estimates, then you should really use one of the other functions instead.

# Steps to implement a new Praat function

1. Indentify what the output of the function will be
    * A signal track (or tracks) that follows the original sound wave
    * A value (or a limited list of values) that summarises the acoustic properties of a wav file, and can therefore not sensibly be shown alongside the sound wave.
2. Implement the core analysis in a Praat script file, and place it in `inst/praat`.
    * In the case where track(s) that follow the sound wave file are returned, the Praat function should write the output to a CSV table file and return the name of that table. The Praat script should also take the desired output table file name (inkluding full path) as an argument. Please refer to `praat/formant_burg.praat` for some example code that computes formants and bandwidths for them for a (possibly windowed) sound file and writes them to a table.
3. Make a copy of the suitable template function, rename it (please keep the praat_ prefix for clarity) and make modifications to the code to suit the new track computed by Praat. You will need to think about what the tracks should be called in the SSFF file and document your choice.
    * For a function that computes a sound wave following signal track (or tracks), use the code of `praat_formant_burg` as a template. Please refer to a suitable function in wrassp for inspiration on what to call sets of tracks. (The `praat_formant_burg` outputs and "fm" and "bw" set, for formant frequencies and formant bandwidths respectivelly)
    * For single value (or list of values) output, there is currently no template function implemented, but please note that the `tjm.praat::wrap_praat_script()` has an option to return the "Info window" of Praat, which opens up lots of possibilities.
4. There are many moving parts to this whole package, so make sure to contruct a test file and a test suit for the new function to make sure that it works. 