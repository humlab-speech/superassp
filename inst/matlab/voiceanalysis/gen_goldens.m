% gen_goldens.m — Generate MATLAB reference outputs for testthat goldens.
% Runs the full Voice Analysis Toolkit pipeline on a1.wav and dumps all
% intermediate signals + final parameters to tests/testthat/golden/a1.mat.
% Also flattens the creak ANN weights to plain matrices so they can be
% loaded from R without the NN toolbox.
%
% Run from repo root via:
%   /Volumes/Fredrik_Nylen_0705853304/Applications/MATLAB_R2025b.app/bin/matlab \
%     -nosplash -nodesktop -nojvm -noawt \
%     -r "run('voiceanalysis/inst/matlab/gen_goldens.m'); exit"

vat_root = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), '..');
% Add only canonical VAT dirs to avoid Octave_VAT/ shadowing main code.
for d = {'SE_VQ', 'creak_fcns', 'general_fcns', 'peakSlope', 'MDQ_recovered'}
    addpath(fullfile(vat_root, d{1}));
end

% Audio
audio_path = '/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/superassp/samples/sustained/a1.wav';
[x, fs0] = audioread(audio_path);
if size(x, 2) > 1, x = mean(x, 2); end
if fs0 ~= 16000
    x = resample(x, 16000, fs0);
end
fs = 16000;
x = x(:);

% Settings (match test_Voice_Analysis_Toolkit.m)
F0min = 20; F0max = 500;

% --- Pipeline stages ---
fprintf('Running SRH pitch tracker...\n');
[f0, VUV, SRHVal, t_f0] = SRH_PitchTracking(x, fs, F0min, F0max);

fprintf('Running creak detector...\n');
[creak_post, creak_dec, t_creak, H2H1, res_p] = CreakyDetection_CompleteDetection(x, fs);
% Also dump the 36-column feature matrix used by the ANN (for parity tests).
creak_FeatTot = get_ALL_creak_features(x, fs);
creak_pp = zeros(1, length(x));
for k = 1:length(creak_dec)
    s = max(1, round((t_creak(k) - 0.005) * fs));
    e = min(length(x), round((t_creak(k) + 0.005) * fs));
    creak_pp(s:e) = creak_dec(k);
end

fprintf('Running SE-VQ (fixed F0) GCI, no creak post-proc...\n');
[GCI, rep, res, MBS] = SE_VQ(x, fs, f0, VUV);
fprintf('  detected %d GCIs\n', length(GCI));

fprintf('Running SE-VQ (var F0) GCI...\n');
try
    [GCI_varF0, ~, ~, f0_samp, ~, ~] = SE_VQ_varF0(x, fs, f0, VUV);
    fprintf('  detected %d GCIs (var F0)\n', length(GCI_varF0));
catch ME
    fprintf('SE_VQ_varF0 failed: %s\n', ME.message);
    GCI_varF0 = GCI;
    f0_samp   = zeros(size(x));
end

fprintf('Running IAIF...\n');
[g_iaif, ar_lpc, e_lpc] = IAIF(x, fs, GCI);

fprintf('Running voice quality...\n');
try
    [NAQ, QOQ, H1H2_param, HRF] = get_NAQ_QOQ_H1H2(g_iaif, fs, GCI);
catch ME
    fprintf('VQ failed (likely findpeaks API mismatch): %s\n', ME.message);
    NAQ = zeros(1, length(GCI));
    QOQ = zeros(1, length(GCI));
    H1H2_param = zeros(1, length(GCI));
    HRF = zeros(1, length(GCI));
end

fprintf('Running peakSlope...\n');
ps = get_peakSlope(x, fs);

fprintf('Running MDQ...\n');
try
    mdq = get_MDQ(res, fs, GCI(:)');
catch ME
    fprintf('MDQ failed: %s\n', ME.message);
    fprintf('GCI length = %d, sample range [%d, %d]\n', length(GCI), min(GCI), max(GCI));
    mdq = zeros(1, length(GCI));
end

% --- Save goldens ---
golden_dir = fullfile(vat_root, 'voiceanalysis', 'tests', 'testthat', 'golden');
if ~exist(golden_dir, 'dir'), mkdir(golden_dir); end

save(fullfile(golden_dir, 'a1.mat'), '-v7', ...
     'x', 'fs', 'f0', 'VUV', 'SRHVal', 't_f0', ...
     'creak_post', 'creak_dec', 't_creak', 'creak_pp', 'H2H1', 'res_p', ...
     'creak_FeatTot', ...
     'GCI', 'GCI_varF0', 'f0_samp', 'rep', 'res', 'MBS', ...
     'g_iaif', 'ar_lpc', 'e_lpc', ...
     'NAQ', 'QOQ', 'H1H2_param', 'HRF', 'ps', 'mdq');
fprintf('Wrote %s/a1.mat\n', golden_dir);

% --- Flatten creak ANN weights to plain arrays ---
ann_dir = fullfile(vat_root, 'creak_fcns');
ANN_net   = load(fullfile(ann_dir, 'SystemNet_creak.mat'));
ANN_maxis = load(fullfile(ann_dir, 'Maxis_creak.mat'));
ANN_minis = load(fullfile(ann_dir, 'Minis_creak.mat'));
net = ANN_net.net;
IW  = cell2mat(net.IW(1,1));
LW  = cell2mat(net.LW(2,1));
b_h = cell2mat(net.b(1));
b_o = cell2mat(net.b(2));
Maxis = ANN_maxis.Maxis(:);
Minis = ANN_minis.Minis(:);

% Capture output mapminmax denormalization parameters (sim() applies this).
out_xmin = NaN; out_xmax = NaN; out_ymin = -1; out_ymax = 1;
try
    pf = net.outputs{2}.processFcns;
    ps = net.outputs{2}.processSettings;
    for k = 1:length(pf)
        if strcmp(pf{k}, 'mapminmax')
            out_xmin = ps{k}.xmin;
            out_xmax = ps{k}.xmax;
            out_ymin = ps{k}.ymin;
            out_ymax = ps{k}.ymax;
            break;
        end
    end
catch
    fprintf('  could not extract output mapminmax params\n');
end
fprintf('  output mapminmax: xmin=%g xmax=%g ymin=%g ymax=%g\n', ...
        out_xmin, out_xmax, out_ymin, out_ymax);

save(fullfile(golden_dir, 'creak_ann_flat.mat'), '-v7', ...
     'IW', 'LW', 'b_h', 'b_o', 'Maxis', 'Minis', ...
     'out_xmin', 'out_xmax', 'out_ymin', 'out_ymax');
fprintf('Wrote %s/creak_ann_flat.mat\n', golden_dir);
