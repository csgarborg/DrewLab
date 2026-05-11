function detectLipFlutterFromExcel(excelPath)

%% Find video files to process

T = readtable(excelPath, 'ReadVariableNames', false);
folderPaths = T{:,1};

parpool; % start parallel pool

for i = 1:length(folderPaths)
    
    folder = folderPaths{i};
    
    if ~isfolder(folder)
        warning('Invalid folder: %s', folder);
        continue;
    end
    
    videoFiles = dir(fullfile(folder, '*.mp4'));
    
    if isempty(videoFiles)
        continue;
    end
    
    fprintf('\n=== Folder: %s ===\n', folder);
    
    % --- Select ROI once per folder ---
    firstVideo = fullfile(folder, videoFiles(1).name);
    [roi1_pos, roi2_pos] = selectROIs(firstVideo);
    
    % --- Parallel video processing ---
    parfor j = 1:length(videoFiles)
        
        videoPath = fullfile(folder, videoFiles(j).name);
        
        try
            processVideoBatched(videoPath, roi1_pos, roi2_pos);
        catch ME
            fprintf('Error: %s\n', videoPath);
            disp(ME.message);
        end
    end
end

end


%% Subfunctions

function processVideoBatched(videoFile, roi1_pos, roi2_pos)

[folder, name, ~] = fileparts(videoFile);
savePath = fullfile(folder, [name '_lipsDiff.mat']);

% ---- SKIP if already processed ----
if exist(savePath, 'file')
    fprintf('Skipping (exists): %s\n', name);
    return;
end

v = VideoReader(videoFile);

% ROI coords
x1 = roi1_pos(1); y1 = roi1_pos(2); w1 = roi1_pos(3); h1 = roi1_pos(4);
x2 = roi2_pos(1); y2 = roi2_pos(2); w2 = roi2_pos(3); h2 = roi2_pos(4);

batchSize = 200;

motionLeft = [];
motionRight = [];
timeVec = [];

idx = 1;

while hasFrame(v)
    
    frames_r1 = [];
    frames_r2 = [];
    t_batch = [];
    
    count = 0;
    
    % ---- LOAD BATCH ----
    while hasFrame(v) && count < batchSize
        
        frame = rgb2gray(readFrame(v));
        
        % Crop immediately
        r1 = frame(y1:y1+h1, x1:x1+w1);
        r2 = frame(y2:y2+h2, x2:x2+w2);
        
        frames_r1(:,:,count+1) = single(r1);
        frames_r2(:,:,count+1) = single(r2);
        
        t_batch(count+1) = v.CurrentTime;
        
        count = count + 1;
    end
    
    if count < 2
        break;
    end
    
    % ---- VECTORIZED DIFF ----
    diff_r1 = abs(diff(frames_r1, 1, 3));
    diff_r2 = abs(diff(frames_r2, 1, 3));
    
    % Mean over space (separate signals)
    motion_r1 = squeeze(mean(mean(diff_r1,1),2));
    motion_r2 = squeeze(mean(mean(diff_r2,1),2));
    
    n = length(motion_r1);
    
    motionLeft(idx:idx+n-1) = motion_r1;
    motionRight(idx:idx+n-1) = motion_r2;
    
    timeVec(idx:idx+n-1) = t_batch(2:end); % align with diff
    
    idx = idx + n;
end

% ---- SAVE ----
lipData.motionValLeft = motionLeft;
lipData.motionValRight = motionRight;
lipData.timeVec = timeVec;
lipData.motionCombined = motionLeft + motionRight;
lipData.motionLeft_dFF = normalizePixelDiff(motionLeft);
lipData.motionRight_dFF = normalizePixelDiff(motionRight);

save(savePath, 'lipData');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eventTimes = detectFlutterEvents(signal, timeVec)

signal = smoothdata(signal, 'gaussian', 5);

thresh = 0.15;
above = signal > thresh;

d = diff([0 above 0]);
startIdx = find(d == 1);
endIdx = find(d == -1) - 1;

eventTimes = [];

minDuration = 0.15;

for i = 1:length(startIdx)
    
    duration = timeVec(endIdx(i)) - timeVec(startIdx(i));
    
    if duration > minDuration
        eventTimes(end+1) = timeVec(startIdx(i));
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str = secToMinSec(t)
m = floor(t/60);
s = t - m*60;
str = sprintf('%02d:%05.2f', m, s);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [roi1_pos, roi2_pos] = selectROIs(videoFile)

v = VideoReader(videoFile);
frame = readFrame(v);
grayFrame = rgb2gray(frame);

figure; imshow(grayFrame);
title('Draw LEFT lip ROI, double-click inside when done');
roi1 = drawrectangle;
wait(roi1);

title('Draw RIGHT lip ROI, double-click inside when done');
roi2 = drawrectangle;
wait(roi2);

roi1_pos = round(roi1.Position);
roi2_pos = round(roi2.Position);

close;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function motion_dFF = normalizePixelDiff(motionSignal)

motionSignal = motionSignal(:); % ensure column

baselineIdx = 1:min(450, length(motionSignal)); % safety check

F0 = mean(motionSignal(baselineIdx));

if F0 == 0
    F0 = eps;
end

motion_dFF = (motionSignal - F0) ./ F0;

end