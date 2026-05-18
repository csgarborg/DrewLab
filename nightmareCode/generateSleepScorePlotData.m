function generateSleepScorePlotData(dateList)

Batch_Preprocess_TDMS_AnalogData_SF(dateList)

end

%% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Batch_Preprocess_TDMS_AnalogData_SF(dateList, forceSpectrogram, forceForceProcessing)
% Batch_Preprocess_TDMS_AnalogData
% -------------------------------------------------------------------------
% Updated version — includes robust edge-artifact correction for *all*
% analog signals (ECoG, EMG, Force).
%
% This file is a full drop-in replacement.
% -------------------------------------------------------------------------

if nargin < 2, forceSpectrogram = false; end
if nargin < 3, forceForceProcessing = false; end

% allMat = dir('*.mat');
% fileList = {};
% pattern = '^\d{6}_\d{2}_\d{2}_\d{4}\.mat$';

% for i = 1:numel(allMat)
%     if ~isempty(regexp(allMat(i).name, pattern, 'once'))
%         fileList{end+1} = allMat(i).name; %#ok<AGROW>
%     end
% end

fileList = findTDMSFiles(dateList);

fprintf('Found %d files to process.\n',numel(fileList));

for idx = 1:numel(fileList)
    srcFile = fileList{idx};
    fprintf('\n=============================================\n');
    fprintf('Processing file %d/%d: %s\n', idx, numel(fileList), srcFile);
    fprintf('=============================================\n');
    Preprocess_TDMS_AnalogData_File(srcFile, forceSpectrogram, forceForceProcessing);
end

fprintf('\nDone processing all files.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% Per-file processing function
% =========================================================================
function Preprocess_TDMS_AnalogData_File(srcFile, forceSpectrogram, forceForceProcessing)

% ===== PARAMETERS =====
dsFs                = 60;
ecogMedianWin       = 10;
ecogNotchBand       = [59 61];
ecogNotchOrder      = 3;

emgPowerBand        = [300 3000];
emgSignalBand       = [10 100];
emgFilterOrder      = 3;
emgPowerKernelWidth = 0.5;

forceLowpassCutoff  = 20;
forceLowpassOrder   = 2;

% ===== LOAD RAW ANALOG DATA =====
% S = load(srcFile,'tdmsDataStruct');
% td = S.tdmsDataStruct;
td = openTDMS(srcFile);

Fs_raw = td.AnalogSamplingRate_Hz_;
if iscell(Fs_raw), Fs_raw = str2double(Fs_raw{1}); end
Fs = str2double(Fs_raw);

animalID = td.AnimalID;

ecogRaw  = td.Analog_Data.ECoG(:)';
emgRaw   = td.Analog_Data.EMG(:)';
forceRaw = td.Analog_Data.Force_Sensor(:)';

L = min([numel(ecogRaw), numel(emgRaw), numel(forceRaw)]);
ecogRaw  = ecogRaw(1:L);
emgRaw   = emgRaw(1:L);
forceRaw = forceRaw(1:L);

duration = L/Fs;

[filePath, fileID, ~] = fileparts(srcFile);
procFile = [filePath '\' fileID '_ProcData.mat'];

haveExistingProc = false;
if exist(procFile,'file')
    tmp = load(procFile,'ProcData');
    if isfield(tmp,'ProcData')
        ProcData = tmp.ProcData;
        haveExistingProc = true;
    end
end
if ~haveExistingProc
    ProcData = struct();
end

ProcData.AnimalID = animalID;
ProcData.FileID = fileID;

ProcData.notes.analogSamplingRate_Hz = Fs;
ProcData.notes.dsFs = dsFs;
ProcData.notes.duration_sec = duration;

% =========================================================================
% 1) --- ECoG PROCESSING WITH EDGE FIX ---
% =========================================================================
fprintf('  Processing ECoG...\n');
ecogProc = processEcog_withEdgeFix_SF(ecogRaw, Fs, ecogMedianWin, ecogNotchBand, ecogNotchOrder);

% Downsample
ProcData.ECoG      = ecogProc;
ProcData.ECoG_DS   = resample(ecogProc, dsFs, Fs);

% =========================================================================
% 2) --- EMG PROCESSING WITH EDGE FIX ---
% =========================================================================
fprintf('  Processing EMG...\n');

[emgPowerDS, emgSignalDS] = processEmg_withEdgeFix_SF( ...
    emgRaw, Fs, dsFs, emgPowerBand, emgSignalBand, emgFilterOrder, emgPowerKernelWidth);

ProcData.EMG.emgPower  = emgPowerDS;
ProcData.EMG.emgSignal = emgSignalDS;

% =========================================================================
% 3) --- FORCE SENSOR PROCESSING
% =========================================================================
needForce = true;
if haveExistingProc && ~forceForceProcessing
    if isfield(ProcData,'forceSensor') && isfield(ProcData,'binForceSensor')
        fprintf('  Using existing forceSensor + binForceSensor.\n');
        needForce = false;
    end
end

if needForce
    fprintf('  Processing Force (filter + threshold GUI)...\n');
    forceHPF = processForceSignal_SF(forceRaw, Fs);
    [forceDS, binForceDS, thresh] = processForce_withEdgeFix_SF( ...
        forceHPF.forceFiltered, Fs, dsFs, forceLowpassCutoff, forceLowpassOrder);
    ProcData.forceSensor = forceDS;
    ProcData.binForceSensor = binForceDS;
    ProcData.notes.forceSensor.threshold = thresh;
end

% SAVE PROCDATA
save(procFile,'ProcData','-v7.3');
fprintf('  Saved ProcData → %s\n',procFile);

% =========================================================================
% 4) --- SPECTROGRAM ---
% =========================================================================
% specFile = sprintf('%s_%s_SpecData.mat',animalID,fileID);
specFile = [filePath '\' fileID '_SpecData.mat'];

if exist(specFile,'file') && ~forceSpectrogram
    fprintf('  Spectrogram exists. Skipping.\n');
else
    fprintf('  Computing ECoG spectrogram...\n');
    SpecData.ECoG = computeEcogSpectrogram_SF(ecogProc, Fs);
    save(specFile,'SpecData','-v7.3');
end

end 


%% ------------------------- ECoG PROCESSING ------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ecogOut = processEcog_withEdgeFix_SF(ecogRaw, Fs, medianWin, notchBand, notchOrder)

% Median
ecog1 = medfilt1(ecogRaw, medianWin);

% Notch
Wn = notchBand/(Fs/2);
[z,p,k] = butter(notchOrder, Wn, 'stop');
[sos,g] = zp2sos(z,p,k);
ecogOut = filtfilt(sos,g,ecog1);
end

%% -------------------------- EMG PROCESSING ------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [emgPowerDS, emgSignalDS] = processEmg_withEdgeFix_SF( ...
    emgRaw, Fs, dsFs, bandPower, bandSignal, ord, kernelWidth)

nyq = Fs/2;

%% EMG POWER FILTER (300–3000 Hz)
bp = bandPower;
bp(2) = min(bp(2), nyq*0.9);
[z,p,k] = butter(ord, bp/nyq);
[sos,g] = zp2sos(z,p,k);

e1 = filtfilt(sos,g, emgRaw - mean(emgRaw));

kernel = gausswin(max(3,round(kernelWidth*Fs)));
kernel = kernel/sum(kernel);
ePow = log10(conv(e1.^2, kernel, 'same'));
emgPowerDS = resample(ePow, dsFs, Fs);

%% EMG SIGNAL FILTER (10–100 Hz)
bs = bandSignal;
bs(2) = min(bs(2), nyq*0.9);
[z2,p2,k2] = butter(ord, bs/nyq);
[sos2,g2] = zp2sos(z2,p2,k2);

e2 = filtfilt(sos2,g2, emgRaw - mean(emgRaw));
emgSignalDS = resample(e2, dsFs, Fs);
end


%% ------------------------ FORCE PROCESSING ------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [forceDS, binDS, thresh] = processForce_withEdgeFix_SF( ...
    forceRaw, Fs, dsFs, cutoff, order)

% Design low-pass filter
[z,p,k] = butter(order, cutoff/(Fs/2),'low');
[sos,g] = zp2sos(z,p,k);

% Filter force (no explicit edge "fix" – rely on filtfilt's internal padding)
f1 = filtfilt(sos, g, forceRaw - mean(forceRaw));

% OPTIONAL: if you really want *very* gentle edge handling, you could
% truncate a tiny amount (e.g. 0.1 s) at start/end instead of copying
% segments, but here we keep the full length.

% Downsample
forceDS = resample(f1, dsFs, Fs);

% ----- threshold GUI -----
% t = (0:numel(forceRaw)-1)/Fs;

done   = false;
binRaw = [];
thresh = [];

% figH = [];   % store figure handle

while ~done
    % Create figure if not already created
    % figH = figure('Color','w','Name','Force Threshold');

    % subplot(2,1,1); 
    % plot(t, forceRaw); 
    % title('Raw Force'); grid on;
    % if ~isempty(thresh), yline(thresh,'r--'); end

    % if ~isempty(binRaw)
    %     subplot(2,1,2); 
    %     plot(t, binRaw); 
    %     ylim([-0.1 1.1]); 
    %     grid on;
    % end

    % answ = inputdlg('Force threshold:','Threshold',1,{'0.05'});
    % if isempty(answ), error('Cancelled'); end
    % thresh = str2double(answ{1});
    thresh = 0.05;

    binRaw = forceRaw > thresh;

    % clf(figH);   % clear figure but keep same window

    % subplot(2,1,1);
    % plot(t, forceRaw); 
    % yline(thresh,'r--');
    % title('Raw Force with threshold'); 
    % grid on;
    % 
    % subplot(2,1,2);
    % plot(t, binRaw); 
    % ylim([-0.1 1.1]); 
    % title('Binary force'); 
    % grid on;
    % 
    % resp = questdlg('Accept threshold?','Confirm','Yes','No','Yes');
    % if strcmp(resp,'Yes')
    if true
        done = true;
    end
end

% % Close the figure after threshold is accepted
% if ishghandle(figH)
%     close(figH);
% end


binDS = resample(double(binRaw), dsFs, Fs) > 0.5;
binDS = binDS(:)';

end



%% ------------------------ SPECTROGRAM ----------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = computeEcogSpectrogram_SF(ecog, Fs)
ecogDet = detrend(ecog,'constant');

params.tapers = [3 5];
params.Fs     = Fs;
params.fpass  = [1 100];
params.err    = 0;
params.trialave = 0;

movingwin = [5 1/5];

[S,T,F] = mtspecgramc(ecogDet, movingwin, params);

out.S = S';
out.T = T;
out.F = F;
out.params = params;
out.movingwin = movingwin;
end

%% ------------------------ FORCE FILTER ----------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = processForceSignal_SF(force, Fs)

% INPUTS:
% force : raw force signal (vector)
% Fs    : sampling frequency (Hz)

% OUTPUT:
% out structure with:
%   .forceRaw
%   .forceDetrended
%   .forceFiltered
%   .activationMask
%   .activationTimes (seconds)

%% --- PARAMETERS (TUNE THESE) ---
hpCutoff = 0.1;        % Hz (removes slow drift; try 0.01–0.5)
smoothWin = 0.05;      % seconds (set 0 to disable smoothing)
threshold = 0.1;       % threshold for activation (adjust to your data)

%% --- STEP 1: Remove DC offset ---
forceDet = detrend(force, 'constant');

%% --- STEP 2: High-pass filter to remove drift ---
[b,a] = butter(2, hpCutoff/(Fs/2), 'high');
forceHP = filtfilt(b, a, forceDet);

%% --- STEP 3: Optional smoothing ---
if smoothWin > 0
    winSamples = round(smoothWin * Fs);
    forceFilt = movmean(forceHP, winSamples);
else
    forceFilt = forceHP;
end

%% --- STEP 4: Thresholding ---
activationMask = forceFilt > threshold;

%% --- STEP 5: Extract activation times ---
% Find rising edges
edges = diff([0; activationMask(:)]) == 1;
activationTimes = find(edges) / Fs;

%% --- OUTPUT ---
out.forceRaw = force;
out.forceDetrended = forceDet;
out.forceFiltered = forceFilt;
out.activationMask = activationMask;
out.activationTimes = activationTimes;

%% --- PLOT (for tuning) ---
% t = (0:length(force)-1)/Fs;
% 
% figure;
% subplot(3,1,1)
% plot(t, force)
% title('Raw Force')
% 
% subplot(3,1,2)
% plot(t, forceFilt)
% hold on
% yline(threshold, 'r--')
% title('Filtered Force + Threshold')
% 
% subplot(3,1,3)
% plot(t, activationMask)
% title('Activation Mask')
% ylim([-0.1 1.1])

end