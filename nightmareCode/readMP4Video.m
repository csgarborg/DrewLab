function videoImages = readMP4Video(saveVideoTF,startTime,stopTime,videoPath)

% Time range (seconds)
if ~exist('stopTime','var')
    disp('Need start and stop time input')
    videoImages = [];
    return
end
tStart = startTime;
tEnd   = stopTime;

% Input video file
if ~exist('videoPath','var')
    [file, loc] = uigetfile('*.mp4');
    videoFile = fullfile(loc,file);
else
    videoFile = videoPath;
end

% Create VideoReader object
v = VideoReader(videoFile);

if saveVideoTF
    % Output folder for frames
    if ~exist(fullfile(loc,'extractedFrames'),'dir')
        mkdir(fullfile(loc,'extractedFrames'))
    end
    outputVideo = fullfile(loc,'extractedFrames',[file(1:end-4) '_trimmed.mp4']);

    % Create VideoWriter object
    vw = VideoWriter(outputVideo, 'MPEG-4');
    vw.FrameRate = v.FrameRate;   % Match original frame rate
    open(vw);
end

% Initialize frame counter
frameIdx = 1;

% Set the current time to the start time
v.CurrentTime = tStart;

% Read and write frames
frameCount = 1;
videoImages = {};
while hasFrame(v) && v.CurrentTime <= tEnd
    frame = readFrame(v);
    if saveVideoTF
        writeVideo(vw, frame);
    end
    videoImages{frameCount} = frame;
    frameCount = frameCount + 1;
end

% Close the video writer
if saveVideoTF
    close(vw);
end

% fprintf('Extracted %d frames from %.1f to %.1f seconds.\n', ...
%     frameIdx-1, tStart, tEnd);
end