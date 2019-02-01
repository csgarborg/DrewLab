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
function motionTrackingCalibration2P(tifFileName,redoFilterTF,micronJumpVal,medFiltTF,saveOutputVideoTF,framesPerSecond,tifFrameBounds)
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
firstTimeTarget = true;
n = 1;
moveDist = [];
velocity = [];
targetPosition = [0,0];
secondsPerFrame = 1/framesPerSecond;
tifLength = tifFrameBounds(2)-tifFrameBounds(1)+1;
imageStack = cell(1,tifLength);

%% Stream Processing Loop
% This is the main processing loop which uses the objects we instantiated
% above to stabilize the input video.
while ~isDone(hVideoSource)
    input = hVideoSource();

    % Find location of Target in the input video frame
    if firstTime
      Idx = int32(pos.template_center_pos);
      motionVector = [0 0];
      firstTime = false;
    else
      IdxPrev = Idx;
% IdxPrev = int32(pos.template_center_pos);

      ROI = [SearchRegion, pos.template_size+2*pos.search_border];
      Idx = hTM(input,Target,ROI);
      
      motionVector = double(Idx-IdxPrev);
    end

%     [Offset, SearchRegion] = updatesearch(sz, motionVector, ...
%         SearchRegion, Offset, pos);
    [Offset] = updatesearch(sz, motionVector, ...
        SearchRegion, Offset, pos);

    % Translate video frame to offset the camera motion
    Stabilized = imtranslate(input, Offset, 'linear');
    
    if firstTimeTarget
        Target = Stabilized(round(TargetRowIndices), round(TargetColIndices));
        firstTimeTarget = false;
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
    
    % Save output video to variable
    if saveOutputVideoTF
        imageStack{1,n} = [input(:,:,1) Stabilized];
        n = n + 1;
    end
    
    % Add pixel motion to data
    motionVector(2) = -motionVector(2);
    moveDist = [moveDist;motionVector];
    velocity = [velocity;sqrt((motionVector(1))^2+(motionVector(2))^2)/secondsPerFrame];
    targetPosition = [targetPosition;targetPosition(end,:)+motionVector];
end

%% Release
% Here you call the release method on the objects to close any open files
% and devices.
release(hVideoSource);

%% Output data

% Save output video to AVI file
if saveOutputVideoTF
    [HOut,WOut] = size(imageStack{1,1});
    imageStackMat = cell2mat(imageStack);
    imageStackMat = double(reshape(imageStackMat,HOut,WOut,tifLength));
    aviFileNameOutput = [tifFileName(1:end-4) '_output.avi'];
    if exist(aviFileNameOutput,'file')
        delete(aviFileNameOutput)
    end
    aviObject = VideoWriter(aviFileNameOutput,'Uncompressed AVI');
    aviObject.FrameRate = 10;
    open(aviObject);
    f = waitbar(0,'Creating output AVI file');
    for k = 1:tifLength
        waitbar(round(k/tifLength,2),f,'Creating output AVI file');
        %     imagesc(filteredData(:,:,k));
        %     filteredFrame = getframe(fig,[0.05 0.05 0.9 0.9]);
        writeVideo(aviObject,imageStackMat(:,:,k));
    end
    close(f)
    close(aviObject);
end

% Get binary ball data to compare to frame movement data
secondsBounds = [tifFrameBounds(1)*secondsPerFrame tifFrameBounds(2)*secondsPerFrame];
ballData = load([tifFileName(1:end-3) 'txt']);
ballDataIndex = secondsBounds(1)<=ballData(:,1) & ballData(:,1)<= secondsBounds(2);
ballData = [ballData(ballDataIndex,1) ballData(ballDataIndex,2)];

figure('Color','White')
subplot(3,1,1)
plot(1:size(moveDist,1),moveDist(:,1),'r')
title('Object Movement Between Frames')
xlabel('Frame')
ylabel('X Movement (Pixels)')
grid on
axis([1 size(moveDist,1) floor(min(moveDist(:,1))/10)*10 ceil(max(moveDist(:,1))/10)*10])
subplot(3,1,2)
plot(1:size(moveDist,1),moveDist(:,2),'b')
title('Object Movement Between Frames')
xlabel('Frame')
ylabel('Y Movement (Pixels)')
grid on
axis([1 size(moveDist,1) floor(min(moveDist(:,2))/10)*10 ceil(max(moveDist(:,2))/10)*10])
subplot(3,1,3)
plot(ballData(:,1),ballData(:,2),'k')
title('Ball Movement')
xlabel('Time (s)')
ylabel('Movement')
grid on
axis([min(ballData(:,1)) max(ballData(:,1)) -1 1])

figure('Color','White')
subplot(2,1,1)
plot(1:length(velocity),velocity,'r')
title('Object Velocity Between Frames')
xlabel('Frame')
ylabel('Velocity (Pixels/s)')
grid on
axis([1 size(velocity,1) floor(min(velocity(:,1))/10)*10 ceil(max(velocity(:,1))/10)*10])
subplot(2,1,2)
plot(ballData(:,1),ballData(:,2),'k')
title('Ball Movement')
xlabel('Time (s)')
ylabel('Movement')
grid on
axis([min(ballData(:,1)) max(ballData(:,1)) -1 1])

frameSelection = [];
cont = true;
figure('Color','White','Name','Select frame starts and stops, select behind line to erase, press enter when finished','NumberTitle','off');
while cont
    subplot(3,1,1)
    plot(1:size(targetPosition,1),targetPosition(:,1),'r')
    hold on
    if length(frameSelection) == 1
        line([frameSelection(1) frameSelection(1)],[min(targetPosition(:,1)) max(targetPosition(:,1))],'Color','g','LineStyle','--');
    end
    for n = 1:2:length(frameSelection)
        line([frameSelection(n) frameSelection(n)],[min(targetPosition(:,1)) max(targetPosition(:,1))],'Color','g','LineStyle','--');
    end
    for n = 2:2:length(frameSelection)
        line([frameSelection(n) frameSelection(n)],[min(targetPosition(:,1)) max(targetPosition(:,1))],'Color','r','LineStyle','--');
    end
    hold off
    title('Object Position per Frame')
    xlabel('Frame')
    ylabel('X Position (Pixels)')
    axis([1 ceil(size(targetPosition,1)/10)*10 floor(min(targetPosition(:,1))/10)*10 ceil(max(targetPosition(:,1))/10)*10])
    grid on
    
    subplot(3,1,2)
    plot(1:size(targetPosition,1),targetPosition(:,2),'b')
    hold on
    if length(frameSelection) == 1
        line([frameSelection(1) frameSelection(1)],[min(targetPosition(:,2)) max(targetPosition(:,2))],'Color','g','LineStyle','--');
    end
    for n = 1:2:length(frameSelection)
        line([frameSelection(n) frameSelection(n)],[min(targetPosition(:,2)) max(targetPosition(:,2))],'Color','g','LineStyle','--');
    end
    for n = 2:2:length(frameSelection)
        line([frameSelection(n) frameSelection(n)],[min(targetPosition(:,2)) max(targetPosition(:,2))],'Color','r','LineStyle','--');
    end
    hold off
    title('Object Position per Frame')
    xlabel('Frame')
    ylabel('Y Position (Pixels)')
    axis([1 ceil(size(targetPosition,1)/10)*10 floor(min(targetPosition(:,2))/10)*10 ceil(max(targetPosition(:,2))/10)*10])
    grid on
    
    subplot(3,1,3)
    plot(ballData(:,1),ballData(:,2),'k')
    title('Ball Movement')
    xlabel('Time (s)')
    ylabel('Movement')
    grid on
    axis([min(ballData(:,1)) max(ballData(:,1)) -1 1])
    
    selectedCoord = ginput(1);
    
    if isempty(selectedCoord)
        cont = false;
    elseif ~isempty(frameSelection) && round(selectedCoord(1)) <= frameSelection(end)
        frameSelection = frameSelection(1:end-1);
    else
        frameSelection(end+1) = round(selectedCoord(1));
    end
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
        avgPixelPos(n,:) = [mean(targetPosition([frameSelection(n*2-1),frameSelection(n*2)],1)), mean(targetPosition([frameSelection(n*2-1),frameSelection(n*2)],2))];
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
    
    save('C:\Workspace\Code\DrewLab\calibrationValues.mat','calibrationValues');
end

figure('Color','White')
k = convhull(targetPosition(:,1),targetPosition(:,2));
plot(targetPosition(k,1),targetPosition(k,2),'b',targetPosition(:,1),targetPosition(:,2),'k');
maxVal = ceil(max(max(abs(targetPosition)))/10)*10;
axis equal square
axis([-maxVal maxVal -maxVal maxVal])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title('Position of Target Object')
xlabel('Pixels')
ylabel('Pixels')
end