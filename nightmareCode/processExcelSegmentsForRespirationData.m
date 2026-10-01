function respirationStruct = processExcelSegmentsForRespirationData(excelPath)

T = readtable(excelPath,'Sheet',2,'ReadVariableNames',false);
if strcmp(T{1,1}{1},'file name')
    T = T(2:end,:);
end
nRec = height(T);

for r = 1:nRec
    
    tdmsPath = T{r,1}{1};
    segments = table2array(T(r,2:end));
    segments = segments(~isnan(segments));
    segments = reshape(segments,2,[])';
    [~,name,~] = fileparts(tdmsPath);
    for n = 1:length(segments(:,1))
        [respFreqCentroidVec,respirationSignal,respirationSignalSection] = respirationSpectrogramPlot_SF(tdmsPath,segments(n,1),segments(n,2));
        respFreqCentroid{n} = [linspace(0,100,length(respFreqCentroidVec))' , respFreqCentroidVec'];
        if n == 1
            respirationStruct.(['file_' name '_fullRaw']) = respirationSignal;
            respirationStruct.(['file_' name '_REMSegmentsRaw']) = {};
        end
        respirationStruct.(['file_' name '_REMSegmentsRaw']){n} = respirationSignalSection;
    end
end


end




function [freq_centroid, signalRaw, signalSection] = respirationSpectrogramPlot_SF(tdmsPath, segmentStart, segmentEnd)
% BANDPASS_PLOT filters a signal between 0.5 and 3 Hz,
% plots the power spectrum (DFT), and plots the spectrogram
% with frequency center of mass overlay.
%
tdmsDataStruct = openTDMS(tdmsPath);

Fs = str2double(tdmsDataStruct.CameraFrameratePerSecond);
signal = tdmsDataStruct.Digital_Data.Respiration_Sum;
signalRaw = signal;
timeVec = (0:length(signal))./Fs;
tMask = timeVec >= segmentStart & timeVec <= segmentEnd;
signal = signal(tMask);
signalSection = signal;

FsTemp = str2double(tdmsDataStruct.TemperatureSamplingRate_Hz_);
signalTemp = tdmsDataStruct.Serial_Data.Temperature;
timeVecTemp = (0:length(signalTemp))./FsTemp;
tMaskTemp = timeVecTemp >= segmentStart & timeVecTemp <= segmentEnd;
meanTemp = mean(signalTemp(tMaskTemp));

% --- Design bandpass filter ---
bpFilt = designfilt('bandpassiir', ...
    'FilterOrder', 4, ...
    'HalfPowerFrequency1', 0.5, ...
    'HalfPowerFrequency2', 3, ...
    'SampleRate', Fs);

% --- Apply zero-phase filter ---
filtered_signal = filtfilt(bpFilt, signal);

% --- Compute DFT ---
N = length(filtered_signal);
X = fft(filtered_signal);
f = (0:N-1)*(Fs/N);
powerX = abs(X).^2 / N;

% figure;
% subplot(4,1,1)
% plot((1:length(signal))./Fs, signal)
% xlabel('Time (s)')
% ylabel('ROI Pixel Sum')
% title('Respiration ROI (Front Chest) Pixel Sum');
% xlim([1 15])
% 
% subplot(4,1,2)
% plot((1:length(filtered_signal))./Fs, filtered_signal)
% xlabel('Time (s)')
% ylabel('Filtered ROI Pixel Sum')
% title('Band-Passed (0.5-3Hz) Respiration ROI (Front Chest) Pixel Sum');
% xlim([1 15])

% % --- Plot Power Spectrum ---
% subplot(4,1,3)
% plot(f, powerX);
% xlim([0 10]);
% xlabel('Frequency (Hz)');
% ylabel('Power');
% title('Power Spectrum of Band-Passed Signal');
% grid on;
% 
% % --- Plot Spectrogram ---
% subplot(4,1,4)

window  = hamming(round(2*Fs));      % 2-second window
overlap = round(0.9 * length(window));
nfft    = 1024;

% Compute spectrogram explicitly (so we can use the data)
[S,F,T] = spectrogram(filtered_signal, window, overlap, nfft, Fs);

% Power spectrogram
P = abs(S).^2;

% % Plot spectrogram
% imagesc(T, F, 10*log10(P));
% axis xy
% colormap jet;
% ax = gca;
% cb = colorbar(ax, 'eastoutside');
% ax.PositionConstraint = 'innerposition';
% ylabel(cb, 'Power')
% title('Spectrogram of Band-Passed Signal');
% xlabel('Time (s)')
% ylabel('Frequency (Hz)')
% ylim([0 5])
% xlim([1 15])

% --- CENTER OF MASS CALCULATION (Frequency Centroid) ---
% Weighted mean frequency at each time slice
freq_centroid = sum(F .* P, 1) ./ sum(P, 1);

% % --- Overlay center of mass ---
% hold on;
% plot(T, freq_centroid, 'w', 'LineWidth', 2);
% plot(T, freq_centroid, 'k--', 'LineWidth', 1); % outline for contrast
% hold off;
end