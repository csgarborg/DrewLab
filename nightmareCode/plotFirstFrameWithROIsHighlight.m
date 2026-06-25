function plotFirstFrameWithROIsHighlight(videoFile,roiPath,highlightROIs)

if nargin < 3
    highlightROIs = {};
end

%% ==============================
%% Read first frame
%% ==============================
v = VideoReader(videoFile);
frame = readFrame(v);

figure
imshow(frame)
title('ROI Locations')
hold on

%% ==============================
%% Load ROIs
%% ==============================
load(roiPath,'roiStruct')

roiMasks  = roiStruct.roiMasks;
roiLabels = roiStruct.roiLabels;

%% ==============================
%% Overlay ROI outlines
%% ==============================
for r = 1:length(roiMasks)

    % Highlight selected ROIs
    if ismember(roiLabels{r},highlightROIs)

        roiColor = 'g';
        textColor = 'g';
        lineWidth = 1.5;

    else

        roiColor = 'r';
        textColor = 'r';
        lineWidth = 1.5;

    end

    % Draw ROI boundary
    B = bwboundaries(roiMasks{r});

    for k = 1:length(B)

        boundary = B{k};

        plot(boundary(:,2),...
             boundary(:,1),...
             'Color',roiColor,...
             'LineWidth',lineWidth);

    end

    % % Label at ROI center
    % [y,x] = find(roiMasks{r});
    % 
    % text(mean(x),...
    %      mean(y)-50,...
    %      roiLabels{r},...
    %      'Color',textColor,...
    %      'FontSize',6,...
    %      'FontWeight','bold',...
    %      'HorizontalAlignment','center');

end

hold off

end