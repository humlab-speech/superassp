#version: 2021-05-06

Add menu command... Objects Help "Momel-Intsint Manual" "" 0 plugin_manual.praat

Add menu command... Objects New "--Momel-Intsint--" "" 0
Add menu command... Objects New "Momel-Intsint" "Momel-Intsint" 0
Add menu command... Objects New "Help" "Momel-Intsint" 1 plugin_manual.praat
Add menu command... Objects New "Create recording directories..." "Momel-Intsint" 1 analysis/create_recording_directories.praat
Add menu command... Objects New "Detect f0..." "Momel-Intsint" 1 analysis/detect_f0.praat
Add menu command... Objects New "Calculate Momel..." "Momel-Intsint" 1 analysis/calculate_momel_targets_extracts.praat
Add menu command... Objects New "Calculate Intsint..." "Momel-Intsint" 1 analysis/calculate_intsint_labels.praat

Add menu command... Objects New "Correct Momel..." "Momel-Intsint" 1 analysis/correct_momel_targets.praat
Add menu command... Objects New "Manipulation..." "Momel-Intsint" 1 analysis/manipulation_momel-intsint.praat

Add menu command... Objects New "Test Momel-Intsint functions" "Momel-Intsint" 1 analysis/test_Momel-Intsint_plugin.praat
Add menu command... Objects New "Reset test data" "Momel-Intsint" 1 analysis/reset_test_data.praat

Add action command... Sound 1 "" 0 "" 0 "To Pitch (auto max/min)..." "To Pitch..." 1 analysis/automatic_min_max_f0.praat
Add action command... Pitch 1 "" 0 "" 0 "Detect Momel anchor pointss..." "" 0 analysis/momel_single_file.praat
Add action command... PitchTier 1 "" 0 "" 0 "Code anchor points with INTSINT..." "" 0 analysis/code_with_intsint.praat
Add action command... Sound 1 PitchTier 1 "" 0 "Draw sound with Momel..." "View & Edit" 0 analysis/draw_sound_f0_&_momel.praat