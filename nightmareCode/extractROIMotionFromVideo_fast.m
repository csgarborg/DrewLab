function roiMotion = extractROIMotionFromVideo_fast(videoFile,excelPath,roiPath,segments)

%% ==============================
%% 1) Load video / cached
%% ==============================
[vidPath, name, ~] = fileparts(videoFile);
[~, roiName, ~] = fileparts(roiPath);
[~, excelName, ~] = fileparts(excelPath);
saveFile = fullfile(vidPath, [name '_' excelName '_' roiName '_roiMotion.mat']);

if exist(saveFile,'file')
    load(saveFile,'roiMotion')
    disp(['ROI data loaded: ' saveFile])
    return
end

v = VideoReader(videoFile);
fprintf('Video loaded: %s\n', videoFile);

%% ==============================
%% 2) Draw / load ROIs
%% ==============================
if exist("roiPath","var") && ~isempty(roiPath)
    load(roiPath);
    roiMasks = roiStruct.roiMasks;
    roiLabels = roiStruct.roiLabels;
else
    frame = readFrame(v);
    figure; imshow(frame);
    title('Draw ROIs → Double-click → Press Enter when done');

    roiMasks = {};
    roiLabels = {};
    roiCount = 0;

    while true
        roi = drawrectangle('Color','r');
        if isempty(roi)
            break;
        end

        roiCount = roiCount + 1;

        label = input(sprintf('Enter label for ROI %d: ', roiCount), 's');
        roiMasks{roiCount} = createMask(roi);
        roiLabels{roiCount} = matlab.lang.makeValidName(label);

        choice = input('Add another ROI? (y/n): ','s');
        if lower(choice) ~= 'y'
            break;
        end
    end
    close
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

for s = 1:nSeg
    
    startTime = segments(s,1);
    stopTime  = segments(s,2);
    
    v.CurrentTime = startTime;
    prevFrame = [];
    
    segData = [];
    segTime = [];
    
    h = waitbar(0, sprintf('Segment %d / %d — 0%%', s, nSeg));
    
    while hasFrame(v) && v.CurrentTime <= stopTime
        
        frame = readFrame(v);
        gray = double(rgb2gray(frame));
        
        % Black frame detection
        if mean(gray(:)) < 5
            prevFrame = [];
            continue
        end
        
        if isempty(prevFrame)
            prevFrame = gray;
            continue
        end
        
        diffFrame = abs(gray - prevFrame);
        
        roiVals = zeros(nROIs,1,'single');
        
        % Vectorized ROI extraction
        for r = 1:nROIs
            roiVals(r) = sum(diffFrame(roiMasks{r}), 'all');
        end
        
        segData(:,end+1) = roiVals; %#ok<AGROW> (small per segment, OK)
        segTime(end+1)   = v.CurrentTime; %#ok<AGROW>
        
        prevFrame = gray;

        % Compute fractional progress for this segment
        frac = (v.CurrentTime - startTime) / (stopTime - startTime);
        frac = max(0, min(1, frac)); % clamp to [0,1]

        % Update waitbar with percent and segment info
        waitbar(frac, h, sprintf('Segment %d / %d — %2.0f%%', s, nSeg, frac*100));
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
%% 5) Save as struct
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

fprintf('Saved: %s\n', saveFile);

end