function procData = smoothBallData(ball,sampleRate)

rawBall = ball(:,2);
t = ball(:,1);
dsFs = 30;   % Hz - downsampled frequency
trialDuration_sec = t(end);   % read this variable in from your data in seconds
analogSamplingRate = sampleRate;   % Hz - change if yours is different
analogExpectedLength = trialDuration_sec*analogSamplingRate;
trimmedEMG = rawBall(1:min(analogExpectedLength,length(rawBall)));
kernelWidth = 0.005;
smoothingKernel = gausswin(kernelWidth*analogSamplingRate)/sum(gausswin(kernelWidth*analogSamplingRate));
ballSmooth = conv(trimmedEMG,smoothingKernel,'same');
resampBall = resample(ballSmooth,dsFs,analogSamplingRate);
procBall = resampBall;   % save this as your final array
procT = 0:1/dsFs:t(end);
if length(procT) > length(procBall)
    procT = procT(1:length(procBall));
elseif length(procT) < length(procBall)
    procBall = procBall(1:length(procT));
end
procData = [procT',procBall];
end