%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:    motionTrackingCalibration2P
%
% FUNCTION:         motionTrackingCalibration2P(tifFileName,redoFilterTF,tifFrameBounds)
%
% DESCRIPTION:      Processes a tif file and extracts movement data based
%                   on user inputs
%
% INPUT:
%
% VARIABLES:
%
% OUTPUT:
%
% FUNCTIONS USED:
%
% LIBARIES USED:
%
% NOTES:
%
% WRITTEN BY:       Spencer Garborg 1/23/19
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function motionTrackingCalibration2P(tifFileName,redoFilterTF,micronJumpVal,medFiltTF,tifFrameBounds)
%% Initialization
close all;

% Input video file which needs to be stabilized.
if exist('tifFrameBounds','var')
    aviFileName = lowpassImageFilter2P(tifFileName,redoFilterTF,medFiltTF,tifFrameBounds);
else
    aviFileName = lowpassImageFilter2P(tifFileName,redoFilterTF,medFiltTF);
    tifFrameBounds = [1 length(imfinfo(tifFileName))];
end

% Get file name
[tokens,~] = regexpi(tifFileName,'\\([^\\]*).TIF','tokens','match');
fileName = tokens{1}{1};

% Create a System object(TM) to read video from a multimedia file. We set the
% output to be of intensity only video.
hVideoSource = vision.VideoFileReader(aviFileName, ...
                                      'ImageColorSpace', 'Intensity',...
                                      'VideoOutputDataType', 'double');

%% Create template
% Create a template matcher System object to compute the location of the
% best match of the target in the video frame. We use this location to find
% translation between successive video frames.
hTM = vision.TemplateMatcher('ROIInputPort', true, ...
                            'BestMatchNeighborhoodOutputPort', true);
                        
%% Create motion tracking output display
% Create a System object to display the original video and the stabilized
% video.
hVideoOut = vision.VideoPlayer('Name', 'Video Stabilization');
hVideoOut.Position(1) = round(0.4*hVideoOut.Position(1));
hVideoOut.Position(2) = round(.5*(hVideoOut.Position(2)));
hVideoOut.Position(3:4) = [1050 550];

%% Initialize variables for processing loop
% Get target window
if medFiltTF
    initialImage = medfilt2(imread(tifFileName, tifFrameBounds(1)));
else
    initialImage = imread(tifFileName, tifFrameBounds(1));
end
imshow(initialImage);
title('Select upper left, then lower right target corners and press enter');
[inputCoordTargetX,inputCoordTargetY] = getpts(gcf);
close(gcf);
inputCoordTargetX = round(inputCoordTargetX);
inputCoordTargetY = round(inputCoordTargetY);
pos.template_orig = [inputCoordTargetX(1) inputCoordTargetY(1)]; % [x y] upper left corner
pos.template_size = [inputCoordTargetX(2) inputCoordTargetY(2)] - [inputCoordTargetX(1) inputCoordTargetY(1)];   % [width height]

% Get search window
imshow(initialImage);
title('Select upper left search area corner (target box pictured) and press enter');
rectangle('Position',[pos.template_orig(1) pos.template_orig(2) pos.template_size(1) pos.template_size(2)],'EdgeColor','w');
[inputCoordSearchX,inputCoordSearchY] = getpts(gcf);
close(gcf);
inputCoordSearchX = round(inputCoordSearchX);
inputCoordSearchY = round(inputCoordSearchY);
pos.search_border = [abs(inputCoordSearchX(1) - inputCoordTargetX(1)),abs(inputCoordSearchY(1) - inputCoordTargetY(1))];   % max horizontal and vertical displacement

% Calculate important parameters
pos.template_center = floor((pos.template_size-1)/2);
pos.template_center_pos = (pos.template_orig + pos.template_center - 1);
fileInfo = info(hVideoSource);
W = fileInfo.VideoSize(1); % Width of video in pixels
H = fileInfo.VideoSize(2); % Height of video in pixels
sz = fileInfo.VideoSize;
% pos.search_border = [W-10 H-10];
BorderCols = [1:pos.search_border(1)+4 W-pos.search_border(1)+4:W];
BorderRows = [1:pos.search_border(2)+4 H-pos.search_border(2)+4:H];
TargetRowIndices = ...
  pos.template_orig(2)-1:pos.template_orig(2)+pos.template_size(2)-2;
TargetColIndices = ...
  pos.template_orig(1)-1:pos.template_orig(1)+pos.template_size(1)-2;
SearchRegion = pos.template_orig - pos.search_border - 1;
Offset = [0 0];
Target = zeros(18,22);
firstTime = true;
n = 1;
MoveDist = [];
TargetPosition = [0,0];

%% Stream Processing Loop
% This is the main processing loop which uses the objects we instantiated
% above to stabilize the input video.
while ~isDone(hVideoSource)
    input = hVideoSource();

    % Find location of Target in the input video frame
    if firstTime
      Idx = int32(pos.template_center_pos);
      MotionVector = [0 0];
      firstTime = false;
    else
      IdxPrev = Idx;
% IdxPrev = int32(pos.template_center_pos);

      ROI = [SearchRegion, pos.template_size+2*pos.search_border];
      Idx = hTM(input,Target,ROI);
      
      MotionVector = double(Idx-IdxPrev);
    end

%     [Offset, SearchRegion] = updatesearch(sz, MotionVector, ...
%         SearchRegion, Offset, pos);
    [Offset] = updatesearch(sz, MotionVector, ...
        SearchRegion, Offset, pos);

    % Translate video frame to offset the camera motion
    Stabilized = imtranslate(input, Offset, 'linear');
    
    if n == 1
        Target = Stabilized(round(TargetRowIndices), round(TargetColIndices));
        n = 0;
    end

    % Add black border for display
    Stabilized(:, BorderCols) = 0;
    Stabilized(BorderRows, :) = 0;

    TargetRect = [pos.template_orig-Offset, pos.template_size];
    SearchRegionRect = [SearchRegion, pos.template_size + 2*pos.search_border];

    % Draw rectangles on input to show target and search region
    input = insertShape(input, 'Rectangle', [TargetRect; SearchRegionRect],...
                        'Color', 'white');
    % Display the offset (displacement) values on the input image
    txt = sprintf('(%+05.1f,%+05.1f)', Offset);
    input = insertText(input(:,:,1),[1 1],txt,'FontSize',16, ...
                    'TextColor', 'white', 'BoxOpacity', 0);
    % Display video
    hVideoOut([input(:,:,1) Stabilized]);
    
    % Add pixel motion to data
    MotionVector(2) = -MotionVector(2);
    MoveDist = [MoveDist;MotionVector];
    TargetPosition = [TargetPosition;TargetPosition(end,:)+MotionVector];
end

%% Release
% Here you call the release method on the objects to close any open files
% and devices.
release(hVideoSource);

%% Output data

% Get binary ball data to compare to frame movement data
% ballDataID = fopen([tifFileName(1:end-3) 'bin']);
% ballData = fread(ballDataID);
% fclose(ballDataID);

subplot(2,1,1)
plot(1:size(MoveDist,1),MoveDist(:,1),'r')
title('Object Movement Between Frames')
xlabel('Frame')
ylabel('X Movement (Pixels)')
grid on
subplot(2,1,2)
plot(1:size(MoveDist,1),MoveDist(:,2),'b')
title('Object Movement Between Frames')
xlabel('Frame')
ylabel('Y Movement (Pixels)')
grid on
% subplot(3,1,3)
% plot(1:size(ballData,1),ballData,'.k')
frameSelection = [];
cont = true;
while cont
    figure(2)
    subplot(2,1,1)
    plot(1:size(TargetPosition,1),TargetPosition(:,1),'r')
    hold on
    if length(frameSelection) == 1
        line([frameSelection(1) frameSelection(1)],[min(TargetPosition(:,1)) max(TargetPosition(:,1))],'Color','g','LineStyle','--');
    end
    for n = 1:2:length(frameSelection)
        line([frameSelection(n) frameSelection(n)],[min(TargetPosition(:,1)) max(TargetPosition(:,1))],'Color','g','LineStyle','--');
    end
    for n = 2:2:length(frameSelection)
        line([frameSelection(n) frameSelection(n)],[min(TargetPosition(:,1)) max(TargetPosition(:,1))],'Color','r','LineStyle','--');
    end
    hold off
    title('Object Position per Frame')
    xlabel('Frame')
    ylabel('X Position (Pixels)')
    axis([-10 ceil(size(TargetPosition,1)/10)*10 floor(min(TargetPosition(:,1))/10)*10 ceil(max(TargetPosition(:,1))/10)*10])
    grid on
    subplot(2,1,2)
    plot(1:size(TargetPosition,1),TargetPosition(:,2),'b')
    hold on
    if length(frameSelection) == 1
        line([frameSelection(1) frameSelection(1)],[min(TargetPosition(:,2)) max(TargetPosition(:,2))],'Color','g','LineStyle','--');
    end
    for n = 1:2:length(frameSelection)
        line([frameSelection(n) frameSelection(n)],[min(TargetPosition(:,2)) max(TargetPosition(:,2))],'Color','g','LineStyle','--');
    end
    for n = 2:2:length(frameSelection)
        line([frameSelection(n) frameSelection(n)],[min(TargetPosition(:,2)) max(TargetPosition(:,2))],'Color','r','LineStyle','--');
    end
    hold off
    title('Object Position per Frame')
    xlabel('Frame')
    ylabel('Y Position (Pixels)')
    axis([-10 ceil(size(TargetPosition,1)/10)*10 floor(min(TargetPosition(:,2))/10)*10 ceil(max(TargetPosition(:,2))/10)*10])
    grid on
    selectedCoord = ginput(1);
    frameSelection(end+1) = round(selectedCoord(1));
    
    if frameSelection(end) < 0
        frameSelection = frameSelection(1:end-1);
        cont = false;
    elseif length(frameSelection) > 1 && frameSelection(end) <= frameSelection(end-1)
        frameSelection = frameSelection(1:end-2);
    end
    % subplot(3,1,3)
    % plot(1:size(ballData,1),ballData,'.k')
end

if exist('C:\Workspace\Code\DrewLab\calibrationValues.mat')
    load('C:\Workspace\Code\DrewLab\calibrationValues.mat');
end

avgPixelPos = [];
if length(frameSelection) < 4
    disp('Could not generate calibration values because too few points were selected')
else
    if length(frameSelection)/2 ~= floor(length(frameSelection)/2)
        frameSelection = frameSelection(1:end-1);
    end
    for n = 1:length(frameSelection)/2
        avgPixelPos(n,:) = [mean(TargetPosition([frameSelection(n*2-1),frameSelection(n*2)],1)), mean(TargetPosition([frameSelection(n*2-1),frameSelection(n*2)],2))];
    end
    for n = 1:size(avgPixelPos,1)-1
        pixelDiff(n) = sqrt(((avgPixelPos(n+1,1)-avgPixelPos(n,1))^2) + ((avgPixelPos(n+1,2)-avgPixelPos(n,2))^2));
    end
    avgPixelDiff = mean(pixelDiff);
    
    calibrationValues.(['file_' fileName]).pixelDiff = pixelDiff;
    calibrationValues.(['file_' fileName]).avgPixelDiff = avgPixelDiff;
    calibrationValues.(['file_' fileName]).micronJumpVal = micronJumpVal;
    calibrationValues.(['file_' fileName]).pixelsPerMicron = avgPixelDiff/micronJumpVal;
    calibrationValues.(['file_' fileName]).micronsPerPixel = micronJumpVal/avgPixelDiff;
    calibrationValues.(['file_' fileName]).frames = tifFrameBounds;
end
save('C:\Workspace\Code\DrewLab\calibrationValues.mat','calibrationValues');

figure(3)
k = convhull(TargetPosition(:,1),TargetPosition(:,2));
plot(TargetPosition(k,1),TargetPosition(k,2),'b',TargetPosition(:,1),TargetPosition(:,2),'k');
maxVal = ceil(max(max(abs(TargetPosition)))/10)*10;
axis equal square
axis([-maxVal maxVal -maxVal maxVal])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title('Position of Target Object')
xlabel('Pixels')
ylabel('Pixels')
end