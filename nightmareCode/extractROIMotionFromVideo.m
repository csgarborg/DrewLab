function extractROIMotionFromVideo(videoFile,roiPath)

%% ==============================
%% 1) Load video
%% ==============================
[vidPath, name, ~] = fileparts(videoFile);
v = VideoReader(videoFile);

fprintf('Loaded video: %s\n', videoFile);

%% ==============================
%% 2) Select ROIs interactively
%% ==============================
if exist("roiPath","var")
    load(roiPath);
    roiMasks = roiStruct.roiMasks; % Assuming roiMasks are stored in roiStruct
    roiLabels = roiStruct.roiLabels; % Assuming roiLabels are stored in roiStruct
    roiList = roiStruct.roiList;

else
    frame = readFrame(v);
    figure; imshow(frame);
    title('Draw ROIs. Double-click inside ROI when done. Press Enter when finished');

    roiList = {};
    roiMasks = {};
    roiLabels = {};

    roiCount = 0;

    while true
        roi = drawrectangle('Color','r');
        if isempty(roi)
            break;
        end

        roiCount = roiCount + 1;

        % Label input
        label = input(sprintf('Enter label for ROI %d: ', roiCount), 's');

        % Create mask
        mask = createMask(roi);

        roiList{roiCount} = roi;
        roiMasks{roiCount} = mask;
        roiLabels{roiCount} = label;

        choice = input('Add another ROI? (y/n): ','s');
        if lower(choice) ~= 'y'
            break;
        end
    end

    close
end

nROIs = length(roiMasks);
fprintf('Total ROIs: %d\n', nROIs);

%% ==============================
%% 3) Reset video + prep
%% ==============================
v.CurrentTime = 0;

prevFrame = [];
motionVals = cell(nROIs,1);
for i = 1:nROIs
    motionVals{i} = [];
end

frameIdx = 0;

%% ==============================
%% 4) Loop through frames
%% ==============================
while hasFrame(v)
    frame = readFrame(v);
    frameIdx = frameIdx + 1;
    
    gray = rgb2gray(frame);
    gray = double(gray);
    
    % Detect black frame (trial separator)
    if mean(gray(:)) < 5
        prevFrame = [];
        continue
    end
    
    if isempty(prevFrame)
        prevFrame = gray;
        continue
    end
    
    diffFrame = abs(gray - prevFrame);
    
    % Compute ROI motion
    for r = 1:nROIs
        roiDiff = diffFrame .* roiMasks{r};
        motionVals{r}(end+1) = sum(roiDiff(:));
    end

    if length(motionVals{1}) > 2
        if motionVals{1}(end) - motionVals{1}(end-1) > (motionVals{1}(end-1) - motionVals{1}(end-2))*5
            subplot(1,2,1)
            imshow(gray)
            subplot(1,2,2)
            imshow(prevFrame)
        end
    end
    
    prevFrame = gray;
end

fprintf('Finished processing %d frames\n', frameIdx);

%% ==============================
%% 5) Package output
%% ==============================
roiMotion = struct;

for r = 1:nROIs
    roiMotion(r).label = roiLabels{r};
    roiMotion(r).motion = motionVals{r};
end

%% ==============================
%% 6) Save file
%% ==============================
saveFile = fullfile(vidPath, [name '_roiMotion.mat']);
save(saveFile, 'roiMotion');

fprintf('Saved: %s\n', saveFile);

end