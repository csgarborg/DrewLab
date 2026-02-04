function mp4ToGif(tStart, tEnd, inputVideoFile, outputGifFile)
% mp4ToGif  Convert MP4 segment to a compressed grayscale GIF

% Settings (unchanged)
numColors = 64;      % Number of grayscale levels
scale     = 1;    % 1 = original size
targetFPS = 20;      % Reduce if source FPS is high

v = VideoReader(inputVideoFile);

tStart = max(tStart, 0);
tEnd   = min(tEnd, v.Duration);
assert(tStart < tEnd, 'Invalid time range.');

% Frame skipping
skip = max(1, round(v.FrameRate / targetFPS));
gifDelay = skip / v.FrameRate;

v.CurrentTime = tStart;

% Read first frame to build grayscale colormap
frame = readFrame(v);
if scale ~= 1
    frame = imresize(frame, scale);
end

grayFrame = rgb2gray(frame);

% Fixed grayscale colormap
cmap = gray(numColors);

% Convert to indexed image using fixed colormap
imgIndexed = gray2ind(grayFrame, numColors);

imwrite(imgIndexed, cmap, outputGifFile, ...
    'gif', 'LoopCount', inf, 'DelayTime', gifDelay);

frameCount = 1;

while hasFrame(v) && v.CurrentTime <= tEnd
    frame = readFrame(v);
    frameCount = frameCount + 1;

    if mod(frameCount, skip) ~= 0
        continue
    end

    if scale ~= 1
        frame = imresize(frame, scale);
    end

    grayFrame = rgb2gray(frame);
    imgIndexed = gray2ind(grayFrame, numColors);

    imwrite(imgIndexed, cmap, outputGifFile, ...
        'gif', 'WriteMode', 'append', 'DelayTime', gifDelay);
end

fprintf('Grayscale GIF written to "%s"\n', outputGifFile);
end