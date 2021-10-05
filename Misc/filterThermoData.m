function procData = filterThermoData(Thermo,sampleRate)

close all

rawThermo = Thermo(:,2);
t = Thermo(:,1);
rawThermo=filloutliers(rawThermo, 'center','mean', 'ThresholdFactor', 4);%this takes out any large artifacts
% process EKG data
% fpass = [300,(sampleRate/2)-1];   % Hz
fpass = [2,40];   % Hz
trialDuration_sec = t(end);   % read this variable in from your data in seconds
analogSamplingRate = sampleRate;   % Hz - change if yours is different
dsFs = 30;   % Hz - downsampled frequency
analogExpectedLength = trialDuration_sec*analogSamplingRate;
trimmedThermo = rawThermo(1:min(analogExpectedLength,length(rawThermo)));
[z,p,k] = butter(5,fpass./(analogSamplingRate/2));
[sos,g] = zp2sos(z,p,k);
filtEKG = filtfilt(sos,g,trimmedThermo - mean(trimmedThermo));
% kernelWidth = 0.5;
kernelWidth = 0.005;
smoothingKernel = gausswin(kernelWidth*analogSamplingRate)/sum(gausswin(kernelWidth*analogSamplingRate));
% EKGPwr = log10(conv(filtEKG.^2,smoothingKernel,'same'));
EKGPwr = conv(filtEKG,smoothingKernel,'same');
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
% plot(procT,procEKG)
% title('Select start and stop time (t) to use to determine mean baseline')
% tVals = ginput(2);
% close
% i = tVals(1,1) <= procT & procT <= tVals(2,1);
% baseline = mean(procEKG(i));
% procEKG = procEKG - (baseline - 1);
procData = [procT',procEKG];

plot(procData(:,1),procData(:,2))
title('Select upper or lower limit')
limitVal = ginput(1);
close
for n = 1:size(procData,1)
    if procData(n,2) > abs(limitVal(2))
        procData(n,2) = abs(limitVal(2));
    elseif procData(n,2) < -abs(limitVal(2))
        procData(n,2) = -abs(limitVal(2));
    end
end
end