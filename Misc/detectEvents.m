%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:    detectEvents
%
% FUNCTION:         detectEvents
%
% DESCRIPTION:      detects running and EMG events from data set
%
% INPUT:
%
% VARIABLES:
%
% OUTPUT:
%
% FUNCTIONS USED:
%
% LIBARIES USED:
%
% NOTES:
%
% WRITTEN BY:       Spencer Garborg 2/22/21
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [motionEvents,stopMotionEvents,EMGEvents,stopEMGEvents,EMGNoMotionEvents] = detectEvents(procBallData,procEMGData,analogSampleRate,secondsPerFrame)

[motionEvents,stopMotionEvents,motionVelocityThresh] = detectMotionEvents(procBallData,analogSampleRate,secondsPerFrame);
[EMGEvents,stopEMGEvents,EMGThresh] = detectEMGEvents(procEMGData,secondsPerFrame);
EMGNoMotionEvents = detectEMGNoMotionEvents(procBallData,procEMGData,analogSampleRate,secondsPerFrame,motionVelocityThresh,EMGThresh);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subfunctions
function [motionEvents,stopMotionEvents,motionVelocityThresh] = detectMotionEvents(procBallData,analogSampleRate,secondsPerFrame)

motionEvents = [];
stopMotionEvents = [];
motionVelocityThresh = [];
analogSampleRate = 1/(procBallData(2,1) - procBallData(1,1));
if size(procBallData,1) <= 7.1*analogSampleRate
    return
else
    procBallData = abs(procBallData);
    plot(procBallData(:,1),procBallData(:,2));
    [~,motionVelocityThresh] = ginput(1);
    close;
    startFrame = 2*analogSampleRate + 1;
    stopFrame = size(procBallData,1) - (3*analogSampleRate) - 1;
    i = find(procBallData(:,2)>motionVelocityThresh);
    for n = 1:length(i)
        index = i(n);
        iStart = index - (2*analogSampleRate);
        iStop = index + (3*analogSampleRate);
        if index > startFrame && index < stopFrame && ~any(procBallData(iStart:index-1,2)>motionVelocityThresh)
            motionEvents(end+1,1:3) = [iStart/analogSampleRate index/analogSampleRate iStop/analogSampleRate];
            motionEvents(end,4:6) = motionEvents(end,1:3)./secondsPerFrame;
            if motionEvents(end,5) - round(motionEvents(end,5)) > 0
                motionEvents(end,4:6) = floor(motionEvents(end,4:6));
            elseif motionEvents(end,5) - round(motionEvents(end,5)) < 0
                motionEvents(end,4:6) = ceil(motionEvents(end,4:6));
            end
        end
    end
    i = find(procBallData(:,2)<motionVelocityThresh);
    for n = 1:length(i)
        index = i(n);
        iStart = index - (2*analogSampleRate);
        iStop = index + (3*analogSampleRate);
        if index > startFrame && index < stopFrame && ~any(procBallData(iStart:index-1,2)<motionVelocityThresh)
            stopMotionEvents(end+1,1:3) = [iStart/analogSampleRate index/analogSampleRate iStop/analogSampleRate];
            stopMotionEvents(end,4:6) = stopMotionEvents(end,1:3)./secondsPerFrame;
            if stopMotionEvents(end,5) - round(stopMotionEvents(end,5)) > 0
                stopMotionEvents(end,4:6) = floor(stopMotionEvents(end,4:6));
            elseif stopMotionEvents(end,5) - round(stopMotionEvents(end,5)) < 0
                stopMotionEvents(end,4:6) = ceil(stopMotionEvents(end,4:6));
            end
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [EMGEvents,stopEMGEvents,EMGThresh] = detectEMGEvents(procEMGData,secondsPerFrame)

EMGEvents = [];
stopEMGEvents = [];
EMGThresh = [];
EMGSampleRate = 1/(procEMGData(2,1) - procEMGData(1,1));
if size(procEMGData,1) <= 7.1*EMGSampleRate
    return
else
    plot(procEMGData(:,1),procEMGData(:,2));
    [~,EMGThresh] = ginput(1);
    close;
    startFrame = 2*EMGSampleRate + 1;
    stopFrame = size(procEMGData,1) - (3*EMGSampleRate) - 1;
    i = find(procEMGData(:,2)>EMGThresh);
    for n = 1:length(i)
        index = i(n);
        iStart = index - (2*EMGSampleRate);
        iStop = index + (3*EMGSampleRate);
        if index > startFrame && index < stopFrame && ~any(procEMGData(iStart:index-1,2)>EMGThresh)
            EMGEvents(end+1,1:3) = [iStart/EMGSampleRate index/EMGSampleRate iStop/EMGSampleRate];
            EMGEvents(end,4:6) = EMGEvents(end,1:3)./secondsPerFrame;
            if EMGEvents(end,5) - round(EMGEvents(end,5)) > 0
                EMGEvents(end,4:6) = floor(EMGEvents(end,4:6));
            elseif EMGEvents(end,5) - round(EMGEvents(end,5)) < 0
                EMGEvents(end,4:6) = ceil(EMGEvents(end,4:6));
            end
        end
    end
    i = find(procEMGData(:,2)<EMGThresh);
    for n = 1:length(i)
        index = i(n);
        iStart = index - (2*EMGSampleRate);
        iStop = index + (3*EMGSampleRate);
        if index > startFrame && index < stopFrame && ~any(procEMGData(iStart:index-1,2)<EMGThresh)
            stopEMGEvents(end+1,1:3) = [iStart/EMGSampleRate index/EMGSampleRate iStop/EMGSampleRate];
            stopEMGEvents(end,4:6) = stopEMGEvents(end,1:3)./secondsPerFrame;
            if stopEMGEvents(end,5) - round(stopEMGEvents(end,5)) > 0
                stopEMGEvents(end,4:6) = floor(stopEMGEvents(end,4:6));
            elseif stopEMGEvents(end,5) - round(stopEMGEvents(end,5)) < 0
                stopEMGEvents(end,4:6) = ceil(stopEMGEvents(end,4:6));
            end
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EMGNoMotionEvents = detectEMGNoMotionEvents(procBallData,procEMGData,analogSampleRate,secondsPerFrame,motionVelocityThresh,EMGThresh)

EMGNoMotionEvents = [];
EMGSampleRate = 1/(procEMGData(2,1) - procEMGData(1,1));
if size(procEMGData,1) <= 7.1*EMGSampleRate
    return
else
    startFrame = 2*analogSampleRate + 1;
    stopFrame = size(procEMGData,1) - (3*analogSampleRate) - 1;
    i = find(procEMGData(:,2)>EMGThresh);
    for n = 1:length(i)
        index = i(n);
        iStart = index - (2*EMGSampleRate);
        iStop = index + (3*EMGSampleRate);
        indexMotion = (analogSampleRate/EMGSampleRate) * index;
        iStartMotion = indexMotion - (2*analogSampleRate);
        iStopMotion = indexMotion + (3*analogSampleRate);
        if index > startFrame && index < stopFrame && ~any(procEMGData(iStart:index-1,2)>EMGThresh) && ~any(procBallData(iStartMotion:iStopMotion,2)>motionVelocityThresh)
            EMGNoMotionEvents(end+1,1:3) = [iStart/analogSampleRate index/analogSampleRate iStop/analogSampleRate];
            EMGNoMotionEvents(end,4:6) = EMGNoMotionEvents(end,1:3)./secondsPerFrame;
            if EMGNoMotionEvents(end,5) - round(EMGNoMotionEvents(end,5)) > 0
                EMGNoMotionEvents(end,4:6) = floor(EMGNoMotionEvents(end,4:6));
            elseif EMGNoMotionEvents(end,5) - round(EMGNoMotionEvents(end,5)) < 0
                EMGNoMotionEvents(end,4:6) = ceil(EMGNoMotionEvents(end,4:6));
            end
        end
    end
end
end