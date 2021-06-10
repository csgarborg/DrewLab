function procData = filterEKGData(EKG,sampleRate,fileName)

close all

rawEKG = EKG(:,2);
t = EKG(:,1);
rawEKG=filloutliers(rawEKG, 'center','mean', 'ThresholdFactor', 4);%this takes out any large artifacts
% process EKG data
% fpass = [300,(sampleRate/2)-1];   % Hz
fpass = [70,200];   % Hz
trialDuration_sec = t(end);   % read this variable in from your data in seconds
analogSamplingRate = sampleRate;   % Hz - change if yours is different
dsFs = 30;   % Hz - downsampled frequency
analogExpectedLength = trialDuration_sec*analogSamplingRate;
trimmedEKG = rawEKG(1:min(analogExpectedLength,length(rawEKG)));
[z,p,k] = butter(5,fpass./(analogSamplingRate/2));
[sos,g] = zp2sos(z,p,k);
filtEKG = filtfilt(sos,g,trimmedEKG - mean(trimmedEKG));
% kernelWidth = 0.5;
kernelWidth = 0.005;
smoothingKernel = gausswin(kernelWidth*analogSamplingRate)/sum(gausswin(kernelWidth*analogSamplingRate));
% EKGPwr = log10(conv(filtEKG.^2,smoothingKernel,'same'));
EKGPwr = conv(filtEKG.^2,smoothingKernel,'same');
% EKGPwr = filtEKG.^2;
resampEKG = resample(EKGPwr,dsFs,analogSamplingRate);
procEKG = resampEKG;   % save this as your final array
procT = 0:1/dsFs:t(end);
% procT = procT(1:end-1);
if length(procT) > length(procEKG)
    procT = procT(1:length(procEKG));
elseif length(procT) < length(procEKG)
    procEKG = procEKG(1:length(procT));
end
plot(procT,procEKG)
title('Select start and stop time (t) to use to determine mean baseline')
tVals = ginput(2);
close
i = tVals(1,1) <= procT & procT <= tVals(2,1);
baseline = mean(procEKG(i));
procEKG = procEKG - (baseline - 1);
procData = [procT',procEKG];


h(1) = figure('Color','White');
subplot(211)
title('EKG')
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
plot(procData(:,1),procData(:,2))
subplot(212)
title('Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('Power')
dt = 1/sampleRate;
tapers = [2 3];
params.Fs = 1/dt;
params.tapers = tapers;
params.fpass = [0.025 100];
params.err = [2 0.05]; %use jack-knife resampling confidence intervals p = 0.05
% [Power, Hz, Error] = mtspectrumc(rawEKG',params);
[Power, Hz, Error] = mtspectrumc(EKGPwr',params);
semilogy(Hz,Power,'k')
[pks,locs] = findpeaks(Power,Hz);
hrSec = pks(5 <= locs & locs <= 15);
hrAmp = max(hrSec);
i = pks == hrAmp;
hold on
plot(locs(i),pks(i),'or')
text(locs(i)+5,pks(i),['HR = ' num2str(locs(i)) ' Hz'])
hold off
savefig(h,fileName);
close
end