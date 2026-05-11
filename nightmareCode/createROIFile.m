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

roiStruct.roiList = roiList;
roiStruct.roiMasks = roiMasks;
roiStruct.roiLabels = roiLabels;
save(roiPath,'roiStruct')