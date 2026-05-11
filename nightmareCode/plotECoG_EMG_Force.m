function plotECoG_EMG_Force(tdmsPath, tStart, tStop)
% plotECoG_EMG_Force
% Plots ECoG spectrogram, EMG, and force signals over a specified time window
%

%% Open TDMS file
tdmsDataStruct = openTDMS(tdmsPath);

%% Extract signals
ecog  = tdmsDataStruct.Analog_Data.ECoG;
emg   = tdmsDataStruct.Analog_Data.EMG;
force = tdmsDataStruct.Analog_Data.Force_Sensor;
fs = str2double(tdmsDataStruct.AnalogSamplingRate_Hz_);

t = (0:length(ecog)-1)/fs;

%% --- Compute spectrogram (Chronux) ---
spec = computeEcogSpectrogram(ecog, fs);

S = spec.S;
T = spec.T;
F = spec.F;

% Convert to dB + normalize
S_db = 10*log10(S + eps);
S_db = S_db - median(S_db(:));

%% --- Plot ---
figure;

% Spectrogram
subplot(3,1,1)
imagesc(T, F, S_db);
axis xy;
ylim([1 100]);

caxis([-30 30]); % adjust if needed
colormap turbo;
colorbar;

ylabel('Frequency (Hz)');
title('ECoG Spectrogram (Chronux)');

% EMG
subplot(3,1,2)
plot(t, emg, 'k');
ylabel('EMG');

% Force
subplot(3,1,3)
plot(t, force, 'b');
ylabel('Force');
xlabel('Time (s)');

% Link axes
linkaxes(findall(gcf,'Type','axes'),'x');

% Apply viewing window only
if nargin == 4
    xlim([tStart tStop]);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = computeEcogSpectrogram(ecog, Fs)

ecogDet = detrend(ecog,'constant');

params.tapers = [3 5];
params.Fs     = Fs;
params.fpass  = [1 100];
params.err    = 0;
params.trialave = 0;

movingwin = [5 0.2]; % [window length, step] = [5 sec, 200 ms]

[S,T,F] = mtspecgramc(ecogDet, movingwin, params);

out.S = S';   % freq x time (better for imagesc)
out.T = T;
out.F = F;
out.params = params;
out.movingwin = movingwin;

end