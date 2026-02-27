function procData = filterNeckEMGData(emg,sampleRate)

rawEMG = emg;
t = length(emg)/sampleRate;
rawEMG=filloutliers(rawEMG, 'center','mean', 'ThresholdFactor', 4);%this takes out any large artifacts
% process EMG data
% fpass = [300,(sampleRate/2)-1];   % Hz
fpass = [300,2999];   % Hz
trialDuration_sec = t(end);   % read this variable in from your data in seconds
analogSamplingRate = sampleRate;   % Hz - change if yours is different
% dsFs = 30;   % Hz - downsampled frequency
analogExpectedLength = trialDuration_sec*analogSamplingRate;
trimmedEMG = rawEMG(1:min(analogExpectedLength,length(rawEMG)));
[z,p,k] = butter(5,fpass./(analogSamplingRate/2));
[sos,g] = zp2sos(z,p,k);
filtEMG = filtfilt(sos,g,trimmedEMG - mean(trimmedEMG));
% kernelWidth = 0.5;
kernelWidth = 0.005;
smoothingKernel = gausswin(kernelWidth*analogSamplingRate)/sum(gausswin(kernelWidth*analogSamplingRate));
EMGPwr = log10(conv(filtEMG.^2,smoothingKernel,'same'));
% EMGPwr = conv(filtEMG.^2,smoothingKernel,'same');
% EMGPwr = filtEMG.^2;
% resampEMG = resample(EMGPwr,dsFs,analogSamplingRate);
% procEMG = resampEMG;   % save this as your final array
% procT = 0:1/dsFs:t(end);
% procT = procT(1:end-1);
procEMG = EMGPwr;
procT = linspace(0,t,length(trimmedEMG));
if length(procT) > length(procEMG)
    procT = procT(1:length(procEMG));
elseif length(procT) < length(procEMG)
    procEMG = procEMG(1:length(procT));
end
emgPlot = figure;
plot(procT,procEMG);
title('Select start and stop time (t) to use to determine mean baseline')
tVals = ginput(2);
close(emgPlot)
i = tVals(1,1) <= procT & procT <= tVals(2,1);
baseline = mean(procEMG(i));
procEMG = procEMG - (baseline - 1);
for n = length(procEMG)-100:length(procEMG)
    if procEMG(n) < 0
        procEMG(n) = procEMG(n-1);
    end
end
% procData = [procT',procEMG];
procData = procEMG;

% subplot(2,1,1)
% plot(emg(:,1),emg(:,2))
% subplot(2,1,2)
% plot(procData(:,1),procData(:,2))
end