%% Video Stabilization
% This example shows how to remove the effect of camera motion from a video stream.

%   Copyright 2006-2014 The MathWorks, Inc.

clc;
close all;
clear;

%% Introduction
% In this example we first define the target to track. In this case, it is the
% back of a car and the license plate. We also establish a dynamic search
% region, whose position is determined by the last known target location.
% We then search for the target only within this search region, which
% reduces the number of computations required to find the target. In each
% subsequent video frame, we determine how much the target has moved
% relative to the previous frame. We use this information to remove
% unwanted translational camera motions and generate a stabilized video. 

%% Initialization
% Create a System object(TM) to read video from a multimedia file. We set the
% output to be of intensity only video.

% Input video file which needs to be stabilized.
tifFileName = 'F:\19-01-18_PaperExp\190118_006.TIF';
redoFilterTF = false;

aviFileName = lowpassImageFilter2P(tifFileName,redoFilterTF,[275 775]);

hVideoSource = vision.VideoFileReader(aviFileName, ...
                                      'ImageColorSpace', 'Intensity',...
                                      'VideoOutputDataType', 'double');

%%
% Create a template matcher System object to compute the location of the
% best match of the target in the video frame. We use this location to find
% translation between successive video frames.
hTM = vision.TemplateMatcher('ROIInputPort', true, ...
                            'BestMatchNeighborhoodOutputPort', true);
                        
%%
% Create a System object to display the original video and the stabilized
% video.
hVideoOut = vision.VideoPlayer('Name', 'Video Stabilization');
hVideoOut.Position(1) = round(0.4*hVideoOut.Position(1));
hVideoOut.Position(2) = round(.5*(hVideoOut.Position(2)));
hVideoOut.Position(3:4) = [1050 550];

%% Here we initialize some variables used in the processing loop.
% Get target window
imshow(imread(tifFileName, 275));
[inputCoordTargetX,inputCoordTargetY] = getpts(gcf);
close(gcf);
pos.template_orig = [inputCoordTargetX(1) inputCoordTargetY(1)]; % [x y] upper left corner
pos.template_size = [inputCoordTargetX(2) inputCoordTargetY(2)] - [inputCoordTargetX(1) inputCoordTargetY(1)];   % [width height]

% Get search window
imshow(imread(tifFileName, 275));
rectangle('Position',[pos.template_orig(1) pos.template_orig(2) pos.template_size(1) pos.template_size(2)],'EdgeColor','w');
[inputCoordSearchX,inputCoordSearchY] = getpts(gcf);
close(gcf);
pos.search_border = [abs(inputCoordSearchX(1) - inputCoordTargetX(1)),abs(inputCoordSearchY(1) - inputCoordTargetY(1))];   % max horizontal and vertical displacement

% Calculate important parameters
pos.template_center = floor((pos.template_size-1)/2);
pos.template_center_pos = (pos.template_orig + pos.template_center - 1);
fileInfo = info(hVideoSource);
W = fileInfo.VideoSize(1); % Width in pixels
H = fileInfo.VideoSize(2); % Height in pixels
sz = fileInfo.VideoSize;
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

%% Stream Processing Loop
% This is the main processing loop which uses the objects we instantiated
% above to stabilize the input video.
n = 1;
MoveDist = [];
TargetPosition = [0,0];
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
        Target = Stabilized(TargetRowIndices, TargetColIndices);
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
ylabel('Pixels (x)')
grid on
subplot(2,1,2)
plot(1:size(MoveDist,1),MoveDist(:,2),'b')
title('Object Movement Between Frames')
xlabel('Frame')
ylabel('Pixels (y)')
grid on
% subplot(3,1,3)
% plot(1:size(ballData,1),ballData,'.k')

figure(2)
subplot(2,1,1)
plot(1:size(TargetPosition,1),TargetPosition(:,1),'r')
title('Object Position per Frame')
xlabel('Frame')
ylabel('Position (Pixels) (x)')
grid on
subplot(2,1,2)
plot(1:size(TargetPosition,1),TargetPosition(:,2),'b')
title('Object Position per Frame')
xlabel('Frame')
ylabel('Position (Pixels) (y)')
grid on
% subplot(3,1,3)
% plot(1:size(ballData,1),ballData,'.k')

figure(3)
plot(TargetPosition(:,1),TargetPosition(:,2));

%% Conclusion
% Using the Computer Vision System Toolbox(TM) functionality from
% MATLAB(R) command line it is easy to implement complex systems like video
% stabilization.

%% Appendix
% The following helper function is used in this example.
%
% * <matlab:edit('updatesearch.m') updatesearch.m>
