function detectREMfromExcel(excelPath, params)

%% =======================
% Default parameters
%% =======================
if nargin < 2 || isempty(params)
    params = defaultREMparams();
end

%% Load Excel (folders)
T = readtable(excelPath, 'ReadVariableNames', false);
folders = T{:,1};

for iF = 1:length(folders)
    folder = folders{iF};
    fprintf('\nProcessing folder: %s\n', folder);

    tdmsFiles = dir(fullfile(folder, '*.tdms'));

    for iFile = 1:length(tdmsFiles)
    
        if contains(tdmsFiles(iFile).name,'converted')
            continue;
        end

        tdmsPath = fullfile(folder, tdmsFiles(iFile).name);
        [~, baseName, ~] = fileparts(tdmsFiles(iFile).name);

        outPath = fullfile(folder, [baseName '_remTimestamps.mat']);

        % Skip if already processed
        if exist(outPath, 'file')
            fprintf('Skipping (exists): %s\n', baseName);
            continue;
        end

        fprintf('Processing file: %s\n', baseName);

        %% =======================
        % Load TDMS
        tdmsDataStruct = openTDMS(tdmsPath);

        %% =======================
        % USER-DEFINED SIGNAL EXTRACTION
        % >>> EDIT THIS FUNCTION <<<
        [ecog, emg, motion, fs] = extractSignalsFromTDMS(tdmsDataStruct);

        %% =======================
        % Feature extraction
        features = computeREMfeatures(ecog, emg, motion, fs, params);

        %% =======================
        % REM detection
        isREM = detectREM(features, params);

        %% =======================
        % Convert to timestamps
        REM_events = getREMtimestamps(isREM, features.time, params);

        %% Save
        save(outPath, 'REM_events');

    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function params = defaultREMparams()

params.win_sec  = 2;     % window size (sec)
params.step_sec = 1;     % step size (sec)

% Threshold percentiles (adaptive)
params.thetaDelta_prctile = 70;
params.emg_prctile        = 30;
params.motion_prctile     = 30;

% Absolute override (optional, set [] to ignore)
params.thetaDelta_thresh = [];
params.emg_thresh        = [];
params.motion_thresh     = [];

% REM rules
params.minREM_sec = 10;

% Frequency bands
params.thetaBand = [6 10];
params.deltaBand = [0.5 4];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ecog, emg, motion, fs] = extractSignalsFromTDMS(tdmsDataStruct)

ecog   = tdmsDataStruct.Analog_Data.ECoG;
emg    = tdmsDataStruct.Analog_Data.EMG;
motion = tdmsDataStruct.Analog_Data.Force_Sensor;

fs = str2double(tdmsDataStruct.AnalogSamplingRate_Hz_);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function features = computeREMfeatures(ecog, emg, motion, fs, params)

win  = round(params.win_sec * fs);
step = round(params.step_sec * fs);

n = floor((length(ecog) - win) / step);

theta_delta = zeros(n,1);
emg_rms     = zeros(n,1);
motion_rms  = zeros(n,1);
timeVec     = zeros(n,1);

for i = 1:n
    idx = (i-1)*step + (1:win);

    ecog_seg = ecog(idx);
    emg_seg  = emg(idx);
    mot_seg  = motion(idx);

    theta = bandpower(ecog_seg, fs, params.thetaBand);
    delta = bandpower(ecog_seg, fs, params.deltaBand);

    theta_delta(i) = theta / (delta + eps);
    emg_rms(i)     = rms(emg_seg);
    motion_rms(i)  = rms(mot_seg);

    timeVec(i) = idx(1) / fs;
end

features.theta_delta = theta_delta;
features.emg_rms     = emg_rms;
features.motion_rms  = motion_rms;
features.time        = timeVec;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isREM = detectREM(features, params)

td = features.theta_delta;
emg = features.emg_rms;
mot = features.motion_rms;

% Adaptive thresholds
TH_td  = getThresh(td, params.thetaDelta_thresh, params.thetaDelta_prctile);
TH_emg = getThresh(emg, params.emg_thresh, params.emg_prctile);
TH_mot = getThresh(mot, params.motion_thresh, params.motion_prctile);

isREM = (td > TH_td) & ...
        (emg < TH_emg) & ...
        (mot < TH_mot);

% Smooth
isREM = medfilt1(double(isREM), 5) > 0;

end


function TH = getThresh(data, manualTH, prct)

if ~isempty(manualTH)
    TH = manualTH;
else
    TH = prctile(data, prct);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function REM_events = getREMtimestamps(isREM, timeVec, params)

d = diff([0; isREM; 0]);

start_idx = find(d == 1);
end_idx   = find(d == -1) - 1;

REM_events = [];

for i = 1:length(start_idx)

    tStart = timeVec(start_idx(i));
    tEnd   = timeVec(end_idx(i));

    if (tEnd - tStart) >= params.minREM_sec

        REM_events = [REM_events; ...
            sec2minsec(tStart), sec2minsec(tEnd)];
    end
end

end


function out = sec2minsec(t)

m = floor(t/60);
s = t - m*60;

out = [m s];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%