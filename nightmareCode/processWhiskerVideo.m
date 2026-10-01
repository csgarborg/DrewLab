function processWhiskerVideo(videoPath, saveVideoFlag, startFrame, stopFrame)
%PROCESS_WHISKER_VIDEO Load an mp4, draw a (rotatable) box on the first
%frame, then loop through every frame processing only the pixels inside
%that box.
%
%   process_whisker_video(videoPath) draws the box and loops through
%   frames without saving anything.
%
%   process_whisker_video(videoPath, true) additionally saves the boxed
%   pixels (after processing) as a new mp4 in the same folder as the
%   original, named "<originalname>_whiskers.mp4".

%   process_whisker_video(videoPath, saveVideoFlag, startFrame, stopFrame)
%   only processes frames startFrame:stopFrame (1-indexed, inclusive)
%   instead of the whole video. Omit or pass [] to use the defaults
%   (startFrame = 1, stopFrame = last frame).
%
%   ROTATION: the drawn box is a drawrectangle ROI with 'Rotatable' set
%   to true, so after drawing you'll see a small rotation handle above
%   the box - drag it to angle the box to match the whisker orientation.
%   Double-click inside the box to confirm and continue.

if nargin < 2
    saveVideoFlag = false;
end
if nargin < 3 || isempty(startFrame)
    startFrame = 1;
end


vr = VideoReader(videoPath);
 
if nargin < 4 || isempty(stopFrame)
    TotalFrameNum = vr.NumFrames;
    stopFrame = TotalFrameNum;
end

firstFrame = readFrame(vr);

fig = figure('Name', 'Draw whisker box');
imshow(firstFrame);
title('Drag to draw box, use rotation handle to angle it, double-click to confirm');

roi = drawrectangle('Rotatable', true, 'Color', 'r');
wait(roi);  % pauses execution until user double-clicks the ROI

% Unrotated position [x y w h] + rotation angle (deg) about the box center
pos    = roi.Position;
angle  = roi.RotationAngle;
boxW   = pos(3);
boxH   = pos(4);
center = [pos(1) + (boxW/2), pos(2) + (boxH/2)];

close(fig);

% Set up the output video writer if requested
if saveVideoFlag
    [folder, name, ~] = fileparts(videoPath);
    outPath = fullfile(folder, [name '_whiskers.mp4']);
    vw = VideoWriter(outPath, 'MPEG-4');
    vw.FrameRate = vr.FrameRate;
    open(vw);
end

% Seek to startFrame. CurrentTime seeking is approximate (nearest
% keyframe/timestamp), so this can occasionally land off by a frame.
if startFrame == 1
    vr.CurrentTime = 0;
else
    vr.CurrentTime = (startFrame - 1) / vr.FrameRate;
end

nFramesToProcess = stopFrame - startFrame + 1;
wb = waitbar(0, sprintf('Processing frame %d/%d...', startFrame, stopFrame), ...
    'Name', 'Processing video');
waitbarUpdateEvery = 30;  % throttle UI updates so they don't slow the loop

firstFrame = true;
frameIdx = startFrame - 1;
while hasFrame(vr) && frameIdx < stopFrame
    frameIdx = frameIdx + 1;
    frame = readFrame(vr);

    if firstFrame
        [tform, xyPadCoords, localCenter] = getTform(frame, center, boxW, boxH, angle);
        firstFrame = false;
    end
    boxPixels = extractRotatedBox(frame, center, boxW, boxH, angle, tform, xyPadCoords, localCenter);

    % ----------------  RADON TRANSFORM GOES HERE  ----------------

    % boxPixels is the (boxH x boxW) image patch for this frame, already de-rotated so it's axis-aligned regardless of the box angle.

    processedPixels = boxPixels;

    % -------------------------------------------------------------

    if saveVideoFlag
        writeVideo(vw, processedPixels);
    end

    nDone = frameIdx - startFrame + 1;
    if mod(nDone, waitbarUpdateEvery) == 0 || frameIdx == stopFrame
        waitbar(nDone / nFramesToProcess, wb, ...
            sprintf('Processing frame %d/%d...', frameIdx, stopFrame));
    end
end

close(wb);

if saveVideoFlag
    close(vw);
    fprintf('Saved processed box video to: %s\n', outPath);
end

fprintf('Processed %d frames.\n', frameIdx);

end

% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tform, xyPadCoords, localCenter] = getTform(frame, center, boxW, boxH, angleDeg)
% For speed, pad enough to cover the box at any angle (its circumscribing circle), plus a couple pixels of margin for interpolation at the edges.
halfDiag = 0.5 * sqrt(boxW^2 + boxH^2);
pad = ceil(halfDiag) + 2;

[H, W, ~] = size(frame);
xMin = max(1, floor(center(1) - pad));
xMax = min(W, ceil(center(1) + pad));
yMin = max(1, floor(center(2) - pad));
yMax = min(H, ceil(center(2) + pad));
xyPadCoords = [xMin xMax yMin yMax];

% Box center's coordinates relative to the cropped subFrame.
localCenter = [center(1) - xMin + 1, center(2) - yMin + 1];

% % Rotate only the small subFrame about localCenter to undo the box angle.
% t1 = [1 0 0; 0 1 0; -localCenter(1) -localCenter(2) 1];
% theta = deg2rad(angleDeg);   % undo the ROI's rotation
% R  = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
% t2 = [1 0 0; 0 1 0; localCenter(1) localCenter(2) 1];
% tform = affine2d(t1 * R * t2);

% Rotate only the small subFrame about localCenter to undo the box angle using newer affinetform2d
t1 = [1 0 -localCenter(1); 0 1 -localCenter(2); 0 0 1];
theta = deg2rad(angleDeg);   % undo the ROI's rotation
R  = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
t2 = [1 0 localCenter(1); 0 1 localCenter(2); 0 0 1];
tform = affinetform2d(t2 * R * t1);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function patch = extractRotatedBox(frame, center, boxW, boxH, angleDeg, tform, xyPadCoords, localCenter)
% Return the boxW x boxH patch of FRAME centered at CENTER and de-rotated by -angleDeg so it's axis-aligned. For speed on large frames with a small box, this only crops+rotates a small local neighborhood around the box instead of warping the entire frame.

% Fast path: no rotation needed, just crop. Computationally free.
if angleDeg == 0
    x1 = round(center(1) - boxW/2);
    y1 = round(center(2) - boxH/2);
    patch = imcrop(frame, [x1, y1, boxW - 1, boxH - 1]);
    return
end

subFrame = frame(xyPadCoords(3):xyPadCoords(4), xyPadCoords(1):xyPadCoords(2), :);

outputView  = imref2d(size(subFrame));
rotatedSub  = imwarp(subFrame, tform, 'OutputView', outputView);

x1 = round(localCenter(1) - boxW/2);
y1 = round(localCenter(2) - boxH/2);
patch = imcrop(rotatedSub, [x1, y1, boxW - 1, boxH - 1]);

end