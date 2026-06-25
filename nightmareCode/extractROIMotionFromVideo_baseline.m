function roiMotion = extractROIMotionFromVideo_baseline(videoFile,excelPath,roiPath,segments)

%% ==============================
%% 1) Load video / cached
%% ==============================
[vidPath, name, ~] = fileparts(videoFile);
[~, roiName, ~] = fileparts(roiPath);
[~, excelName, ~] = fileparts(excelPath);

if size(segments,2) == 3
    % wake event
    saveFile = fullfile(vidPath, [name '_' excelName '_wakeEvents_' roiName '_roiMotion_baseline.mat']);
    baselineSegmentStart = segments(:,1);
    segments = segments(:,[2 3]);
else
    saveFile = fullfile(vidPath, [name '_' excelName '_' roiName '_roiMotion_baseline.mat']);
    baselineSegmentStart = segments(:,1);
    segments = segments(:,[1 2]);
end

if exist(saveFile,'file')
    load(saveFile,'roiMotion')
    disp(['ROI data loaded: ' saveFile])
    return
end

v = VideoReader(videoFile);
fprintf('Video loaded: %s\n', videoFile);

%% ==============================
%% 2) Load ROIs
%% ==============================
if exist("roiPath","var") && ~isempty(roiPath)
    load(roiPath);
    roiMasks  = roiStruct.roiMasks;
    roiLabels = roiStruct.roiLabels;
else
    error('ROI path required for this version.')
end

nROIs = length(roiMasks);

%% ==============================
%% 3) Process segments
%% ==============================
if ~exist("segments","var") || isempty(segments)
    segments = [0 v.Duration];
end

nSeg = size(segments,1);
motionCells = cell(nSeg,1);
timeCells   = cell(nSeg,1);

baselineWindow = 5; % seconds

for s = 1:nSeg

    startTime = segments(s,1);
    stopTime  = segments(s,2);

    fprintf('\n==============================\n');
    fprintf('Segment %d / %d\n', s, nSeg);
    fprintf('REM segment: %.2f to %.2f sec\n', ...
        startTime, stopTime);

    %% =========================================
    %% 3A) Determine baseline window
    %% =========================================

    startTimeBaseline = baselineSegmentStart(s);

    if startTimeBaseline >= baselineWindow
        % Preferred:
        % use 5 sec BEFORE REM onset
        baselineStart = startTimeBaseline - baselineWindow;
        baselineStop  = startTimeBaseline;

        fprintf('Using PRE-REM baseline: %.2f to %.2f sec\n', ...
            baselineStart, baselineStop);

    else
        % Fallback:
        % too close to file start → use first 5 sec after segment starts
        baselineStart = startTimeBaseline;
        baselineStop  = min(startTimeBaseline + baselineWindow, stopTime);

        fprintf(['Segment too close to file start.\n' ...
                 'Using POST-start fallback baseline: %.2f to %.2f sec\n'], ...
                 baselineStart, baselineStop);
    end

    %% =========================================
    %% 3B) Build baseline image
    %% =========================================

    v.CurrentTime = baselineStart;

    baselineSum   = [];
    baselineCount = 0;

    while hasFrame(v) && v.CurrentTime <= baselineStop

        frame = readFrame(v);
        gray = double(rgb2gray(frame));

        % Skip black frames
        if mean(gray(:)) < 5
            continue
        end

        if isempty(baselineSum)
            baselineSum = zeros(size(gray));
        end

        baselineSum = baselineSum + gray;
        baselineCount = baselineCount + 1;
    end

    if baselineCount == 0
        warning('No valid baseline frames found for segment %d', s);
        continue
    end

    baselineFrame = baselineSum / baselineCount;

    fprintf('Baseline built using %d frames\n', baselineCount);

    %% =========================================
    %% 3C) Process REM motion vs baseline
    %% =========================================

    v.CurrentTime = startTime;

    segData = [];
    segTime = [];

    h = waitbar(0, ...
        sprintf('Segment %d / %d — 0%%', s, nSeg));

    while hasFrame(v) && v.CurrentTime <= stopTime

        frame = readFrame(v);
        gray = double(rgb2gray(frame));

        % Skip black frames
        if mean(gray(:)) < 5
            continue
        end

        % Difference from baseline
        diffFrame = abs(gray - baselineFrame);

        roiVals = zeros(nROIs,1,'single');

        for r = 1:nROIs
            roiVals(r) = sum(diffFrame(roiMasks{r}), 'all');
        end

        segData(:,end+1) = roiVals; %#ok<AGROW>
        segTime(end+1)   = v.CurrentTime; %#ok<AGROW>

        % Progress bar
        frac = (v.CurrentTime - startTime) / ...
               (stopTime - startTime);
        frac = max(0, min(1, frac));

        waitbar(frac, h, ...
            sprintf('Segment %d / %d — %2.0f%%', ...
            s, nSeg, frac*100));
    end

    motionCells{s} = segData;
    timeCells{s}   = segTime;

    close(h);
end

%% ==============================
%% 4) Concatenate segments
%% ==============================
motionMat = cat(2, motionCells{:});
t         = cat(2, timeCells{:});

%% ==============================
%% 5) Save output struct
%% ==============================
roiMotion = struct;

for r = 1:nROIs
    roiMotion.(strrep(roiLabels{r},' ','_')) = motionMat(r,:);
end

roiMotion.time = t;

%% ==============================
%% 6) Save
%% ==============================
save(saveFile, 'roiMotion', '-v7.3');

fprintf('\nSaved: %s\n', saveFile);

end