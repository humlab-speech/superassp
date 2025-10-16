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
  "wrassp::forest"=wrassp::forest(file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE,verbose=FALSE),
  "superassp::forest"=superassp::forest(file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE,verbose=FALSE),
  "superassp::praat_formant_burg"=superassp::praat_formant_burg(file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE),
 "superassp::praat_formantpath_burg"=superassp::praat_formantpath_burg(file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE),
 "superassp::praat_sauce"=superassp::praat_sauce(file.path(getwd(),"tests/signalfiles/msajc003.wav"),toFile=FALSE),
 
 times=10
)
```

which results in 

```
Unit: milliseconds
                              expr       min        lq      mean    median        uq       max neval
                    wrassp::forest  21.78219  21.82512  22.56147  22.46408  22.73331  24.64198    10
                 superassp::forest  32.05052  32.26893  32.67404  32.63813  33.05957  33.28942    10
     superassp::praat_formant_burg 566.44296 592.93524 595.09192 596.84491 601.73589 609.05947    10
 superassp::praat_formantpath_burg 486.35340 487.72140 489.20697 489.29277 491.10624 492.59876    10
            superassp::praat_sauce 569.54039 571.38022 576.49244 574.12995 579.52352 590.93038    10
       
```
The `uperassp::forest` performs the same signal processing as `wrassp::forest` but can take any media file supported by the `av` R package, and the conversion of data takes some additional time.
The  `praat_formant_burg` and `praat_formantpath_burg` functions has same flexible media loading ability, but calls praat processing routines through python and the parselmouth module, whcih takes some extra time.

The `praat_sauce` computes many extra 

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
   "kaldi pitch tracker"=kaldi_pitch(
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
                    expr          min           lq         mean       median           uq          max neval
           wrassp::ksvF0     2.233792     2.281709     2.335946     2.299771     2.428917     2.459084    10
           wrassp::mhsF0    16.071417    16.318417    16.747580    16.629147    16.911459    18.199459    10
                     dio    96.165001    98.599834   113.813180   103.385646   109.233959   196.728834    10
                    rapt   117.702001   123.111418   130.628734   125.266938   132.193834   172.862917    10
                   swipe   132.854251   142.303001   154.309442   146.920792   149.362251   203.044126    10
                     yin   174.342459   176.832625   184.944901   179.869314   185.566626   227.427792    10
     kaldi pitch tracker   195.539792   198.606459   217.634276   203.737271   211.397167   336.594000    10
                  reaper   232.866792   240.993793   817.140772   242.397271   244.391001  5989.057042    10
                 harvest   327.358626   334.249459   341.861072   337.557001   353.730001   368.915709    10
                    pyin   403.986376   410.861001   439.136880   440.595814   455.601042   500.009751    10
     praat_pitch ac & cc   516.680042   536.404792   543.467217   546.130938   552.376292   557.215459    10
                   crepe   566.955126   605.125417   853.842680   612.232230   650.363376  3008.081876    10
                   yaapt   526.269001   575.632584   632.245809   626.762709   662.692500   849.112584    10
 praat_pitch all methods 13428.942209 13825.150959 13894.481517 13884.433500 13964.615626 14224.678042    10
```
I have rearranged the output so that the algorithms are roughly ordered by (median) time used to compute output tracks.

Please note that these relative timings are not necessarily indicative of the relative efficiency of the algorithms themselves.
The communication between R and Praat / python has a severe impact on performance, so the benchmarks above indicate only the relative performance in the current version of `superassp`. 

It should also be noted that as the computation is already slow due to the process of calling Praat the `superassp` functions instead takes the opportunity to return more information once processing a file. For instance, `praat_pitch` returns up to two or four tracks in which f_0 was estimated and may therefore be worth the wait. The `swipe` estimates an additional "pitch" track, and `reaper` and `kaldi_pitch` computes and returns also normalized cross-correlation.

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
