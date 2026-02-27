function results = compareUltrasoundAudio(wavOff, wavOn, tStartOff, tStopOff, tStartOn, tStopOn)

% COMPARE_RIG_STATES
% Ultrasonic comparison between Rig OFF and Rig ON
% Adds:
%   - ΔSPL for all significant ultrasonic peaks
%   - Spectrogram of ON - OFF difference

%% ---- Defaults ----
if nargin < 3, tStartOff = []; end
if nargin < 4, tStopOff  = []; end
if nargin < 5, tStartOn  = []; end
if nargin < 6, tStopOn   = []; end

%% ---- Load ----
[xOff, Fs1] = loadInput(wavOff);
[xOn,  Fs2] = loadInput(wavOn);

if Fs1 ~= Fs2
    error('Sample rates do not match.');
end

Fs = Fs1;
nyquist = Fs/2;
maxFreq = min(100e3, nyquist);

%% ---- Time Crop ----
xOff = cropSignal(xOff, Fs, tStartOff, tStopOff);
xOn  = cropSignal(xOn,  Fs, tStartOn,  tStopOn);

%% ---- Spectral Params ----
windowLength = 2048;
overlap = 0.75;
nfft = 4096;

window = hann(windowLength);
noverlap = round(windowLength * overlap);

[Soff,F,T] = spectrogram(xOff, window, noverlap, nfft, Fs);
[Son,~,~]  = spectrogram(xOn,  window, noverlap, nfft, Fs);

Poff = abs(Soff).^2;
Pon  = abs(Son).^2;

meanOff = mean(Poff, 2);
meanOn  = mean(Pon, 2);

meanOff_dB = 10*log10(meanOff + eps);
meanOn_dB  = 10*log10(meanOn  + eps);

delta_dB = meanOn_dB - meanOff_dB;

%% ---- Calibration ----
S_effective = 0.100;   % V/Pa (v1.2.0 medium gain approx)
Vref = 0.707;
Pref = 20e-6;

%% ---- 20–100 kHz Band ΔSPL ----
bandIdx = F >= 20e3 & F <= maxFreq;

bandOff = sqrt(mean(meanOff(bandIdx)));
bandOn  = sqrt(mean(meanOn(bandIdx)));

SPLoff = 20*log10((bandOff*Vref/S_effective)/Pref);
SPLon  = 20*log10((bandOn *Vref/S_effective)/Pref);

deltaSPL_band = SPLon - SPLoff;

%% ---- Detect Ultrasonic Peaks (ON spectrum) ----

ultraMask = F >= 20e3 & F <= maxFreq;
ultraIdx = find(ultraMask);

meanOnUltra_dB = meanOn_dB(ultraMask);

[pks, locs] = findpeaks(meanOnUltra_dB, ...
                        'MinPeakProminence', 2, ...
                        'MinPeakDistance', 20);

peakFreqs = F(ultraIdx(locs));
deltaSPL_peaks = zeros(size(peakFreqs));

%% ---- Compute ΔSPL Using Band RMS Around Each Peak ----

bw = 500;  % ±500 Hz integration window
for k = 1:length(peakFreqs)
    idx = F > (peakFreqs(k)-bw) & F < (peakFreqs(k)+bw);
    toneOff = sqrt(mean(meanOff(idx)));
    toneOn  = sqrt(mean(meanOn(idx)));
    SPLoff_t = 20*log10((toneOff*Vref/S_effective)/Pref);
    SPLon_t  = 20*log10((toneOn *Vref/S_effective)/Pref);
    deltaSPL_peaks(k) = SPLon_t - SPLoff_t;
end

%% ---- Overlay Plot ----
figure;
subplot(2,1,1)
plot(F/1000, meanOff_dB, 'LineWidth',1.5);
hold on;
plot(F/1000, meanOn_dB,  'LineWidth',1.5);
xlim([0 maxFreq/1000]);

xlabel('Frequency (kHz)');
ylabel('Power (dB)');
title('Time-Averaged Spectrum: Rig OFF vs ON');
legend({'Rig OFF','Rig ON'}, 'Location','best');
grid on;

%% ---- Annotate Peak ΔSPL (RMS Based) ----

if ~isempty(peakFreqs)
    for k = 1:length(peakFreqs)
        text(peakFreqs(k)/1000, ...
             interp1(F, meanOn_dB, peakFreqs(k)) + 2, ...
             sprintf('\\Delta %.1f dB', deltaSPL_peaks(k)), ...
             'HorizontalAlignment','center', ...
             'FontWeight','bold');
    end
end

%% ---- Add Summary Box (Upper Left) ----

summaryText = sprintf(['\\DeltaSPL 20–100 kHz: %.2f dB\n' ...
                       'Strongest Peak: %.2f kHz\n' ...
                       '\\Delta Peak (max): %.2f dB'], ...
                       deltaSPL_band, ...
                       peakFreqs(1)/1000, ...
                       max(deltaSPL_peaks));

text(0.02*maxFreq/1000, ...
     max(meanOn_dB)-1, ...
     summaryText, ...
     'VerticalAlignment','top', ...
     'BackgroundColor','white', ...
     'EdgeColor','black');

%% ---- Difference Spectrum Plot ----
subplot(2,1,2)
plot(F/1000, delta_dB, 'LineWidth',1.5);
xlim([0 maxFreq/1000]);
xlabel('Frequency (kHz)');
ylabel('\Delta Power (dB)');
title('Difference Spectrum (ON - OFF)');
grid on;

%% ---- Cleaned Difference Spectrogram ----

% Convert to dB first
Soff_dB = 10*log10(Poff + eps);
Son_dB  = 10*log10(Pon  + eps);

% Smooth across time (reduces flicker noise)
Soff_dB_s = movmean(Soff_dB, 5, 2);
Son_dB_s  = movmean(Son_dB,  5, 2);
diffSpec = Son_dB_s - Soff_dB_s;
figure;
imagesc(T, F/1000, diffSpec);
axis xy;
ylim([0 maxFreq/1000]);
xlabel('Time (s)');
ylabel('Frequency (kHz)');
title('Spectrogram Difference (ON - OFF)');
colorbar;
clim([-15 15]);   % tighter dynamic range = cleaner view
colormap turbo;

%% ---- Stable Difference Spectrogram ----

Soff_mean = mean(10*log10(Poff + eps), 2);
Son_mean  = mean(10*log10(Pon  + eps), 2);

diffSpec_stable = Son_mean - Soff_mean;

% Replicate across time for visualization
diffSpec_plot = repmat(diffSpec_stable, 1, length(T));

figure;
imagesc(T, F/1000, diffSpec_plot);
axis xy;
ylim([0 maxFreq/1000]);

xlabel('Time (s)');
ylabel('Frequency (kHz)');
title('Time-Averaged Spectral Difference (ON - OFF)');
colorbar;
clim([-10 10]);
colormap turbo;

%% ---- Console Output ----
fprintf('\n20–100 kHz Band ΔSPL: %.1f dB\n', deltaSPL_band);

for k = 1:length(peakFreqs)
    fprintf('Peak at %.2f kHz: ΔSPL = %.1f dB\n', ...
        peakFreqs(k)/1000, deltaSPL_peaks(k));
end

%% ---- Return ----
results.deltaSPL_20_100kHz = deltaSPL_band;
results.peakFrequencies_kHz = peakFreqs/1000;
results.deltaSPL_peaks = deltaSPL_peaks;

end

%% ---- Helper Functions ----
function x = cropSignal(x, Fs, tStart, tStop)
if isempty(tStart) || isempty(tStop), return; end
startIdx = max(1, round(tStart*Fs));
stopIdx  = min(length(x), round(tStop*Fs));
x = x(startIdx:stopIdx);
end

function [x, Fs] = loadInput(input)
[x, Fs] = audioread(input);
if size(x,2) > 1
    x = mean(x,2);
end
end