x = load('C:\Workspace\210218_008.txt');
emg = x(:,[1 3]);
sampleRate = 10000;

close all
rawEMG = emg(:,2);
t = emg(:,1);
% process EMG data
% fpass = [300,(sampleRate/2)-1];   % Hz
fpass = [300,3000];   % Hz
trialDuration_sec = t(end);   % read this variable in from your data in seconds
analogSamplingRate = sampleRate;   % Hz - change if yours is different
dsFs = 30;   % Hz - downsampled frequency
analogExpectedLength = trialDuration_sec*analogSamplingRate;
trimmedEMG = rawEMG(1:min(analogExpectedLength,length(rawEMG)));
[z,p,k] = butter(3,fpass/(analogSamplingRate/2));
[sos,g] = zp2sos(z,p,k);
filtEMG = filtfilt(sos,g,trimmedEMG - mean(trimmedEMG));
% kernelWidth = 0.5;
kernelWidth = 0.01;
smoothingKernel = gausswin(kernelWidth*analogSamplingRate)/sum(gausswin(kernelWidth*analogSamplingRate));
% EMGPwr = log10(conv(filtEMG.^2,smoothingKernel,'same'));
EMGPwr = conv(filtEMG.^2,smoothingKernel,'same');
resampEMG = resample(EMGPwr,dsFs,analogSamplingRate);
procEMG = resampEMG;   % save this as your final array
procT = 0:1/dsFs:t(end);
% procT = procT(1:end-1);
semilogy(procT,procEMG)
title('Select start and stop time (t) to use to determine mean baseline')
tVals = ginput(2);
close
i = tVals(1,1) <= procT & procT <= tVals(2,1);
baseline = mean(procEMG(i));
procEMG = procEMG - (baseline - 1);
procData = [procT',procEMG];

emgSimpleFilt = bandpass(rawEMG,[300 3000],sampleRate);

subplot(3,1,1)
plot(t,rawEMG)
title('Raw EMG')

subplot(3,1,2)
plot(t,emgSimpleFilt)
title('Simple bandpass 300-3000 Hz')

subplot(3,1,3)
plot(procData(:,1),procData(:,2))
title('Filter Code at 300-3000 Hz w/ no log10 after smoothing')