

form Ange parametrar till analysen
	sentence SoundFile /Users/frkkan96/Desktop/aaa.wav
	real Start_time_(s) 0.0
	real End_time_(s) 0.0
	comment Parameters for extracting PowerCepstrogram
	real left_Pitch_range_(Hz) 60.0
	real right_Pitch_range_(Hz) 333.0
	real Time_step_(s) 0.002
	boolean Get_CPPS_for_each_time_step 0
	real Maximum_frequency_(Hz) 5000.0
	real Preemphasis_from_(Hz) 50
	comment Windowing when extracting part of a file
	#word Window_type Hanning
	optionmenu Window_type 4
		button rectangular
		button triangular
		button parabolic
		button Hanning
		button Hamming
		button Gaussian1
		button Gaussian2
		button Gaussian3
		button Gaussian4
		button Gaussian5
		button Kaiser1
		button Kaiser2
	real Window_width 1.0
	comment CPPS extraction paramers
	boolean	Detrend 1
	real Time_averaging_window_(s) 0.02
	real Quefrency_averaging_window_(s) 0.0005
	real Tolerance_(0-1) 0.05
	optionmenu Interpolation 2
				button None
				button Parabolic
				button Cubic
				button Sinc70
	real left_Trend_line_quefrency_range_(s) 0.001
	real right_Trend_line_quefrency_range_(s) 0.05
	optionmenu Trend_type 2
				button Straight
				button Exponential decay
	optionmenu Fit_method 1
				button Robust
				button Least squares
				button Robust slow	
endform

#For benchmarking
#stopwatch

interp1$ = "None"
interp2$ = "Parabolic"
interp3$ = "Cubic"
interp4$ = "Sinc70"

trendt1$ = "Straight"
trendt2$ = "Exponential decay"

fitm1$ = "Robust"
fitm2$ = "Least squares"
fitm3$ = "Robust slow"

wt1$ = "rectangular"
wt2$ = "triangular"
wt3$ = "parabolic"
wt4$ = "Hanning"
wt5$ = "Hamming"
wt6$ = "Gaussian1"
wt7$ = "Gaussian2"
wt8$ = "Gaussian3"
wt9$ = "Gaussian4"
wt10$ = "Gaussian5"
wt11$ = "Kaiser1"
wt12$ = "Kaiser2"

i$ = interp'interpolation'$
t$ = trendt'trend_type'$
f$ = fitm'fit_method'$
w$ = wt'window_type'$

currSound = Read from file: soundFile$
soundStart = Get start time
soundEnd = Get end time

if start_time = 0.0
	start_time = soundStart 
endif

if end_time = 0.0
	end_time = soundEnd 
endif

outStr$ = ""

Extract part: start_time, end_time, w$, window_width, 0
currPowerCepstr= noprogress To PowerCepstrogram: left_Pitch_range, time_step, maximum_frequency, preemphasis_from

if get_CPPS_for_each_time_step == 0

	cpps = noprogress Get CPPS: detrend, time_averaging_window, quefrency_averaging_window, left_Pitch_range, right_Pitch_range, tolerance, i$, left_Trend_line_quefrency_range, right_Trend_line_quefrency_range, t$, f$
	outStr$ = "'start_time';'cpps'"
else
	select currPowerCepstr
	noFrames = Get number of frames
	for frame from 1 to noFrames
		if outStr$ != ""
			outStr$ = outStr$ + ";"
		endif
		select currPowerCepstr
		frameTime = Get time from frame number: frame
		cepstrumSlice = To PowerCepstrum (slice): frameTime
		cpps = Get peak prominence: left_Pitch_range, right_Pitch_range, i$, 0.001, 0.05, t$, f$
		outStr$ = outStr$ + "'frameTime';'cpps'"
	endfor
endif

writeInfoLine: outStr$
#For benchmarking in Praat 
#time = stopwatch
#appendInfoLine: "('i$','t$','f$'), Time='time's : CPPS='outStr$'"


