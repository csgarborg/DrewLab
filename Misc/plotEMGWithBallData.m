function plotEMGWithBallData(filename)
close all

ballData = load(filename);

t = ballData(:,1);
rot = ballData(:,2);
emg = ballData(:,3);

rawEMG = emg;
% process EMG data
fpass = [100,400];   % Hz
trialDuration_sec = t(end);   % read this variable in from your data in seconds
analogSamplingRate = 1000;   % Hz - change if yours is different
dsFs = 30;   % Hz - downsampled frequency
analogExpectedLength = trialDuration_sec*analogSamplingRate;
trimmedEMG = rawEMG(1:min(analogExpectedLength,length(rawEMG)));
[z,p,k] = butter(3,fpass/(analogSamplingRate/2));
[sos,g] = zp2sos(z,p,k);
filtEMG = filtfilt(sos,g,trimmedEMG - mean(trimmedEMG));
kernelWidth = 0.5;
smoothingKernel = gausswin(kernelWidth*analogSamplingRate)/sum(gausswin(kernelWidth*analogSamplingRate));
EMGPwr = log10(conv(filtEMG.^2,smoothingKernel,'same'));
resampEMG = resample(EMGPwr,dsFs,analogSamplingRate);
procEMG = resampEMG;   % save this as your final array
procT = 0:1/dsFs:t(end);
semilogy(procT,procEMG)
tVals = ginput(2);
close
i = tVals(1,1) <= procT & procT <= tVals(2,1);
baseline = mean(procEMG(i));
procEMG = procEMG - (baseline - 1);

h(1) = figure('Color','White');
x1 = subplot(3,1,1);
plot(t,emg,'k')
title('\fontsize{20pt}\bf{Raw Abdominal EMG Recording}')
xlabel('Time (s)')
ylabel('V')
grid on
axis([min(t) max(t) floor(min(emg)) ceil(max(emg))])

x2 = subplot(3,1,2);
semilogy(procT,procEMG,'k')
title('\fontsize{20pt}\bf{Filtered Abdominal EMG Recording}')
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
grid on
axis([min(procT) max(procT) floor(min(procEMG)) ceil(max(procEMG))])


x3 = subplot(3,1,3);
plot(t,convertBallVoltToMPS(rot),'k')
title('\fontsize{20pt}\bf{Ball Movement}')
xlabel('Time (s)')
ylabel('Movement')
grid on
axis([min(t) max(t) -1 ceil(max(convertBallVoltToMPS(rot)))])
end