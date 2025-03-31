function procData = smoothBallData(ball,sampleRate)
close all
for n = 2
    rawBall = ball(:,2);
    t = ball(:,1);
    fpass = [1,50];   % Hz
    dsFs = 30;   % Hz - downsampled frequency
    trialDuration_sec = t(end);   % read this variable in from your data in seconds
    analogSamplingRate = sampleRate;   % Hz - change if yours is different
    analogExpectedLength = trialDuration_sec*analogSamplingRate;
    trimmedBall = rawBall(1:min(analogExpectedLength,length(rawBall)));
    if n == 2
        [z,p,k] = butter(5,fpass./(analogSamplingRate/2));
        [sos,g] = zp2sos(z,p,k);
        filtBall = filtfilt(sos,g,trimmedBall - mean(trimmedBall));
%         filtBall = filtfilt(sos,g,trimmedBall);
    else
        filtBall = trimmedBall;
    end
    kernelWidth = 0.005;
    smoothingKernel = gausswin(kernelWidth*analogSamplingRate)/sum(gausswin(kernelWidth*analogSamplingRate));
    ballSmooth = conv(filtBall,smoothingKernel,'same');
    resampBall = resample(ballSmooth,dsFs,analogSamplingRate);
    procBall = resampBall;   % save this as your final array
    procT = 0:1/dsFs:t(end);
    if length(procT) > length(procBall)
        procT = procT(1:length(procBall));
    elseif length(procT) < length(procBall)
        procBall = procBall(1:length(procT));
    end
    procData = [procT',procBall];
%     plot(procT',abs(procBall))
%     plot(procT',procBall)
%     hold on
end
end