function roiStruct = createROIFile(videoFile,roiPath)

%% ==============================
%% 1) Load video / cached
%% ==============================

if exist(roiPath,'file')
    disp(['ROI selections exist: ' roiPath])
    return
end

v = VideoReader(videoFile);
fprintf('Video loaded: %s\n', videoFile);

%% ==============================
%% 2) Draw / load ROIs
%% ==============================
frame = readFrame(v);
figure; imshow(frame);

roiList = {};
roiMasks = {};
roiLabels = {};

roiCount = 0;

roiLabels = {'Eye top', 'Eye bottom', 'Nose', 'Right ear bottom', 'Right ear top', 'Left ear bottom', 'Left ear top', 'Right lip bottom front', 'Right lip top front', 'Left lip bottom front', 'Left lip top front', 'Lower jaw bottom', 'Right whisker rear', 'Right whisker front', 'Left whisker rear', 'Left whisker front'};
while true
    if roiCount < length(roiLabels)
        title(['Draw ' roiLabels{roiCount+1} ' ROI. Double-click inside ROI when done. Press Enter when finished.']);
    else
        title('Draw next ROI. Double-click inside ROI when done. Press Enter when finished.');
    end
    roi = drawrectangle('Color','r');
    if isempty(roi)
        break;
    end

    roiCount = roiCount + 1;

    % Label input
    if roiCount <= length(roiLabels)
        label = roiLabels{roiCount};
    else
        label = input(sprintf('Enter label for ROI %d: ', roiCount), 's');
    end

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

roiStruct.roiList = roiList;
roiStruct.roiMasks = roiMasks;
roiStruct.roiLabels = roiLabels;
save(roiPath,'roiStruct')