#!/usr/bin/perl

# intsint	perl script
$version = 2.11 ;

#	author: Daniel Hirst
#	email:  daniel.hirst@lpl.univ-aix.fr

# purpose: 
#			input list of F0 targets as first argument
#			output list of optimised INTSINT targets and estimated target
#
#	input format [cb]:
#				time (ms) target (Hz)
#
#	output format [int]:
#				time (secs) tone (INTSINT) target (Hz) estimated_target (Hz)

$verbose = 1;

#input and ouput filename extensions
$inext = "momel";
$outext ="intsint";

#parameters for data checking
$MIN_F0=60;	#Hz
$MAX_F0=600;	#Hz

#parameters for optimisation
$MIN_PAUSE = 0.5;	#seconds
$MIN_RANGE = 0.5;	#octaves
$MAX_RANGE = 2.5;	#octaves
$STEP_RANGE = 0.1;	#octaves
$MEAN_SHIFT = 50;	#Hertz
$STEP_SHIFT=1;		#Hertz
$BIG_NUMBER = 9999;

#parameters for target estimation (distance towards T or B respectively)
$HIGHER = 0.5;
$LOWER = 0.5;
$UP = 0.25;
$DOWN = 0.25;

#define path separator
$PS = "\/";
if ($^O eq "MacOS") {$PS = ":"};
($script_name = $0) =~ s/^.*$PS//;
$verbose and print "Perl script $script_name version [$version] running under $^O\n";

# input file 1st argument
$infile = $ARGV[0];

$now = localtime;
if ($infile !~ /\.$inext/) {
			die("Error on input file - it should have the extension \[\.$inext\]\n");
}
open (IN, $infile) or die("cannot open $infile");
($outfile = $infile) =~ s/$inext$/$outext/;
open (OUT, ">$outfile");
$verbose and print "reading from $infile\n";
$verbose and print "writing to $outfile\n";

# name of input and output files
($infile_name =$infile) =~ s/^.*$PS//;
($outfile_name =$outfile) =~ s/^.*$PS//;


$nval=0;
$sum_f0=0;

LINE: while (<IN>) {
			($t, $f0) = split;
#make sure F0 is withing limits
			if($f0 < $MIN_F0) {
				$f0 = $MIN_F0;
			} elsif ($f0 > $MAX_F0 ) {
			    $f0 = $MAX_F0;
			}# if
			$f0[$nval] = &octave($f0); # f0 converted to octave scale
			$t[$nval] = $t/1000; # t converted to seconds
			$sum_f0 += $f0[$nval];			
			$nval++;
} #LINE: while

$mean_f0 = $sum_f0/$nval;

$linear_mean_f0 = &round(&linear($mean_f0));
print OUT "; $outfile_name created on $now by $script_name $version\n";
print OUT "; from $infile_name\n";
print OUT ";    $nval values  mean = $linear_mean_f0\n";

$min_mean = $linear_mean_f0 - $MEAN_SHIFT;
$max_mean = $linear_mean_f0 + $MEAN_SHIFT;
$min_ss_error = $BIG_NUMBER;
@list_tones = ("M", "S", "T",  "B", "H", "L", "U", "D");


$last_estimate = 0;

$verbose and print "optimising range between $MIN_RANGE and $MAX_RANGE (in octaves) by step $STEP_RANGE\n";
$verbose and print "optimising mid between $min_mean and $max_mean (in Hz) by step $STEP_SHIFT\n";

for ($range = $MIN_RANGE; $range < $MAX_RANGE; $range += $STEP_RANGE) {

			for ($lm = $min_mean; $lm < $max_mean; $lm += $STEP_SHIFT) {
#$verbose and print "range $range ; key $lm \n";
						$mid = &octave($lm);
						&optimise($mid, $range);
			} #for $lm

} #for $range

$mid = $best_mid;
$key = &linear($mid);
$bottom = $mid - $best_range/2;
$top = $mid + $best_range/2;
print OUT "<parameter range=$best_range>\n";
print OUT "<parameter key=$key>\n";

for $i (0..$nval-1) {
			$f0 = int(linear($f0[$i])+0.5);
			$estimate = int(linear($best_estimate[$i])+0.5);
			print OUT "$t[$i] $best_intsint[$i] $f0 $estimate\n";
} #for $i

sub optimise {
				my ($mid, $range) = @_;
				$top = $mid + $range/2;
				$bottom = $mid - $range/2;
				$f0 = $f0[0];

				if ($top-$f0 < abs($f0-$mid)) {
						$intsint[0] = "T";
				} elsif ($f0-$bottom < abs($f0-$mid)) {
						$intsint[0] = "B";
				} else { $intsint[0] = "M" }

				$estimate = &estimate($intsint[0], $last_estimate);
				$estimate[0] = $estimate;
				$error = abs($estimate-$f0[0]);
				$ss_error = $error*$error;
				$last_estimate = $estimate;

				for $i (1..$nval-1) {
							$target = $f0[$i];

# after pause choose from [MTB]
							if ($t[$i]-$t[$i-1] > $MIN_PAUSE) { 
										if ($top-$target < abs($target-$mid)) {
														$intsint[$i] = "T";
										} elsif ($target-$bottom < abs($target-$mid)) {
														$intsint[$i] = "B";
										} else { $intsint[$i] = "M" }

#elsewhere any tone except M
							} else { 
										$min_difference = $BIG_NUMBER;
										$best_tone = "";
										foreach $tone (@list_tones) {
														unless ($tone eq "M") {
														    $estimate = &estimate($tone, $last_estimate);
																		$difference = abs($target-$estimate);										
														    if ($difference < $min_difference) {
																						$min_difference = $difference;
																						$best_tone = $tone;
														    } #if
														} #unless
										} #foreach $tone

										$intsint[$i] = $best_tone;
							} #else

							$estimate[$i] = &estimate($intsint[$i], $last_estimate);
							$error = abs($estimate[$i]-$f0[$i]);
							$ss_error += $error*$error;
							$last_estimate = $estimate;
				} #for $i

				if ($ss_error < $min_ss_error) {
						$min_ss_error = $ss_error;
						$best_range = $range;
						$best_mid = $mid;
						@best_intsint = @intsint;
						@best_estimate = @estimate;
   	} #if

} #sub optimise

sub estimate {
		my ($tone, $last_target) = @_;
		SWITCH: for ($tone) {
				/M/ && do {$estimate = $mid; last SWITCH;};
				/S/ && do {$estimate = $last_target; last SWITCH;};
				/T/ && do {$estimate = $top; last SWITCH;};
				/H/ && do {$estimate = $last_target+($top-$last_target)*$HIGHER; last SWITCH;};
				/U/ && do {$estimate = $last_target + ($top-$last_target)*$UP; last SWITCH;};
				/B/ && do {$estimate = $bottom; last SWITCH;};
				/L/ && do {$estimate = $last_target - ($last_target - $bottom)*$LOWER; last SWITCH;};
				/D/ && do {$estimate = $last_target - ($last_target - $bottom)*$DOWN; last SWITCH;};
		} #SWITCH
		return $estimate;
} #sub

sub octave {
		my ($value) = @_;
		return log($value)/log(2);
}
sub linear {
		my ($value) = @_;
		return 2**$value;
}

sub round {
		my ($value) = @_;
		return int($value+0.5);
}

#version history
#2.11	29nov2006 changed input and output extensions to ".momel" and ".intsint"
#2.10	16may2003 targets (.cb) read in ms and converted to seconds
#2.9    16may2003 added constraint on F0 targets forced within interval [min_F0, max_f0]
#2.8 	13may2003 corrected error in minpause (was still in ms) (merci Robert!)
#2.7 	25jun2002 L and H set to halfway towards B and T
#				  D and U quarter of way towards  B and T (=halfway towards L and H)
#				  changed order of output : time tone target estimate
#2.6 	15may2002	fixed bug which allowed H or U after T and L or D after B
#2.5	9may2002	minor changes - output file extension ".int"
#2.3	04sep2000	time values in seconds
#2.2 	19aug2000	equal space in range given to T, H, U, S, D, L, B on octave scale
#						optimise put in subroutine
#2.1 	02aug2000	changed all f0 calculations to octave scale
#2.0	19jul2000  	Perl droplet working
#1.2	06sep1995	C version intsint.c 


