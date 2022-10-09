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
 praat_formantpath_burg=praat_formantpath_burg(
    file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE),
  praat_sauce=praat_sauce(
    file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE),
 times=100
)
```

which results in 

```
Unit: milliseconds
                   expr        min         lq       mean     median        uq       max neval
         wrassp::forest   26.42033   28.25807   28.93884   28.74673   29.4339   34.7792   100
     praat_formant_burg  520.85644  556.63421  596.37337  578.95567  628.1025  777.4130   100
 praat_formantpath_burg  669.06082  708.40986  751.45798  733.72906  776.1413 1170.8230   100
            praat_sauce 3247.77570 3400.72715 3668.55957 3577.48971 3895.4160 4753.4219   100
             
```
Getting an SSFF file from a `wrassp` function rather than the `praat_formant_burg` function, which is wrapped call of Praat call and which also involves the parsing of a csv file. Since the parsing of input and output in the `praat_formant_burg` Praat calls already slows computation down considerably, the function also computes formant amplitudes (L) before returning the output to increase the usefulness of the function. The `praat_formantpath_burg` function is of course an additional bit slower than method of computing formant frequencies as multiple formant tracks are computed and compared when this function is used. 

Also, even it is adviced that even though the functions `praat_sauce` does compute formant tracks (F and B properties) as well, it is not really efficiently implemented and is really mostly there for correction of harmonic amplitudes. And, an additional factor to consider is that the formant tracks will be stored by the `praat_sauce` function in a field in the same file as all the other tracks computed by the function, which will likely result in a performance issue when working with the tracks.

So, if you need only formant frequency and bandwidth estimations, then you should really use one of the other functions instead.

Similarly, f0 computation using functions that call Praat or python are considerably slower than their `wrassp` counterparts:

```r
library(microbenchmark)
microbenchmark(
  "wrassp::ksvF0"=wrassp::ksvF0(
    file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE,windowShift=5),
  "wrassp::mhsF0"=wrassp::mhsF0(
    file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE,windowShift=5),
  "praat_pitch ac & cc"=praat_pitch(
    file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE,corr.only=TRUE,windowShift=5),
   "praat_pitch all methods"=praat_pitch(
    file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE,corr.only=FALSE,windowShift=5),
    "rapt"=rapt(
    file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE,windowShift=5),
  "swipe"=swipe(
    file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE,windowShift=5),
  "reaper"=reaper(
    file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE,windowShift=5),
  "yin"=yin(
    file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE,windowShift=5),
  "pyin"=pyin(
    file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE,windowShift=5),
  "dio"=dio(
    file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE,windowShift=5),
  "crepe"=crepe(
    file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE,windowShift=5),
  "harvest"=harvest(
    file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE,windowShift=5),
  "yaapt"=yaapt(
    file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE,windowShift=5),
 times=10
) 
```

```
Unit: milliseconds
                    expr         min           lq         mean      median           uq          max neval
           wrassp::ksvF0     2.25950     2.270251     2.284413     2.28298     2.285875     2.349209    10
           wrassp::mhsF0    15.88817    16.153917    16.310271    16.29860    16.481251    16.852375    10
                     dio    62.42633    62.860292    96.337913    65.58179    68.612251   373.232251    10
                    rapt    84.97538    87.076334    93.813997    88.50365    90.148542   146.366625    10
                   swipe   102.42483   102.572876   109.453142   104.94073   108.750542   148.430667    10
                     yin   144.30888   145.829543   153.094672   146.54019   156.877501   188.481876    10
                  reaper   201.77229   203.018459   210.588767   204.25125   208.676584   262.837876    10
                 harvest   301.28638   301.704042   307.238076   303.94225   306.866293   337.149209    10
                    pyin   373.87946   374.840542   855.573847   377.29808   383.740626  5146.458709    10
                   yaapt   464.19521   483.067584   545.585967   507.16013   567.194834   797.772959    10
                   crepe   551.81225   559.429667   779.727342   566.27196   583.744167  2668.921042    10
     praat_pitch ac & cc   621.59875   628.186917   640.165342   634.75581   650.736709   688.405375    10
 praat_pitch all methods 13851.93675 13872.469501 13927.357080 13894.94540 13964.398751 14078.464918    10
```
I have rearranged the output so that the algorithms are roughly ordered by (median) time used to compute output tracks.

Please note that these relative timings are not necessarily indicative of the relative efficiency of the algorithms themselves.
The communication between R and Praat / python has a severe impact on performance, so the benchmarks above indicate only the relative performance in the current version of `superassp`. 

It should also be noted that as the computation is already slow due to the process of calling Praat the `superassp` functions instead takes the opportunity to return more information once processing a file. For instance, `praat_pitch` returns up to two or four tracks in which f_0 was estimated and may therefore be worth the wait. The `swipe` estimates an additional "pitch" track, and `reaper` computes and returns also normalized cross-correlation.

# Steps to implement a new Praat function

1. Indentify what the output of the function will be
    * A signal track (or tracks) that follows the original sound wave
    * A value (or a limited list of values) that summarises the acoustic properties of a wav file, and can therefore not sensibly be shown alongside the sound wave.
2. Implement the core analysis in a Praat script file, and place it in `inst/praat`.
    * In the case where track(s) that follow the sound wave file are returned, the Praat function should write the output to a CSV table file and return the name of that table. The Praat script should also take the desired output table file name (including full path) as an argument. Please refer to `praat/formant_burg.praat` for some example code that computes formants and bandwidths for them for a (possibly windowed) sound file and writes them to a table.
3. Make a copy of the suitable template function, rename it (please keep the praat_ prefix for clarity) and make modifications to the code to suit the new track computed by Praat. You will need to think about what the tracks should be called in the SSFF file and document your choice.
    * For a function that computes a sound wave following signal track (or tracks), use the code of `praat_formant_burg` as a template. Please refer to a suitable function in wrassp for inspiration on what to call sets of tracks. (The `praat_formant_burg` outputs and "fm" and "bw" set, for formant frequencies and formant bandwidths respectivelly)
    * For single value (or list of values) output, there is currently no template function implemented, but please note that the `tjm.praat::wrap_praat_script()`, which `cs_wrap_praat_script` is a revised version of, has an option to return the "Info window" of Praat, which opens up lots of possibilities.
4. There are many moving parts to this whole package, so make sure to contruct a test file and a test suit for the new function to make sure that it works. 
