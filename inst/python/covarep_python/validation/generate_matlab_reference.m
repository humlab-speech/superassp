% MATLAB Reference Output Generator for Python Validation
% Run this in MATLAB COVAREP to generate reference data

% Setup COVAREP path
run('../startup.m');

% Load test audio
[wave, fs] = audioread('0011.arctic_bdl1.wav');

%% F0 Tracking with SRH
fprintf('Computing F0 with SRH...\n');
f0min = 50;
f0max = 500;
hopsize = 10; % ms

[F0s, VUVDecisions, SRHVal, time] = pitch_srh(wave, fs, f0min, f0max, hopsize);

% Save F0 results
f0_data = [time(:), F0s(:), double(VUVDecisions(:)), SRHVal(:)];
save('matlab_f0_reference.mat', 'F0s', 'VUVDecisions', 'SRHVal', 'time');
dlmwrite('matlab_f0_reference.txt', f0_data, 'delimiter', '\t', 'precision', 6);
fprintf('✓ Saved MATLAB F0 reference\n');

%% IAIF on sample frame
fprintf('\nComputing IAIF on sample frame...\n');
% Extract 30ms frame at 1 second
start_sample = round(1.0 * fs);
frame_len = round(0.03 * fs);
frame = wave(start_sample:start_sample+frame_len-1);

% Run IAIF
[g, dg, a, ag] = iaif(frame, fs);

% Save IAIF results
iaif_data = struct('g', g, 'dg', dg, 'a', a, 'ag', ag, 'fs', fs);
save('matlab_iaif_reference.mat', '-struct', 'iaif_data');

% Save as text too
dlmwrite('matlab_iaif_glottal_flow.txt', g, 'delimiter', '\t', 'precision', 8);
dlmwrite('matlab_iaif_flow_derivative.txt', dg, 'delimiter', '\t', 'precision', 8);

fprintf('✓ Saved MATLAB IAIF reference\n');
fprintf('\n');
fprintf('Reference files created:\n');
fprintf('  - matlab_f0_reference.mat/.txt\n');
fprintf('  - matlab_iaif_reference.mat\n');
fprintf('  - matlab_iaif_glottal_flow.txt\n');
fprintf('  - matlab_iaif_flow_derivative.txt\n');
