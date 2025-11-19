%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:    motionTracking2P
%
% FUNCTION:         motionTracking2P(tifFileName,calibrationFileString,redoFilterTF,updateSearchTF,medFiltTF,saveOutputVideoTF,threeStepTF,compiledTifTF,targetAvgNum,framesPerSecond,objMag,digMag,turnabout,commentString,tifFrameBounds)
%
% DESCRIPTION:      Processes a tif file and extracts movement data based
%                   on user inputs
%
% INPUT:            tif must be extracted from MDF file and have the same name as ball data text file (also extracted from MDF file as ANSI encoding (8-bit))
%
% VARIABLES:        tifFileName: (string) complete path to tif data file
%                   calibrationFileString: (string) the filename of the tif file that was processed through the calibration code to be used for this file (ex. 190315_001)
%                   redoFilterTF: (logical) select false when data have already been filtered and output AVI created but want to reprocess without redoing the filter
%                   updateSearchTF: (logical) false if the searched ROI should remain in the same selected place for tracking or true if should update every frame based on target location
%                   medFiltTF: (logical) true if spatial median filter should be applied to the image data
%                   saveOutputVideoTF: (logical) true if tracking output video should be saved to AVI file
%                   threeStepTF: (logical) true if speed of tracking code is more important than accuracy - not every point of the cross-correlation is calculated and position of target is narrowed down quickly
%                   compiledTifTF: (logical) true if tif file is larger than 4 GB and had to be stitched together in ImageJ from multiple seperate tif files
%                   targetAvgNum: (double) number of targets to average together at the beginning of the tif to get final target used to match and track object
%                   framesPerSecond: (double) frames per second value
%                   objMag: (double) objective magnification value
%                   digMag: (double) digital magnification value (in MSCAN)
%                   turnabout: (double) turnabout value in MSCAN (affects calibration of image)
%                   commentString: (string, NOT REQUIRED) include any notes about this run to be displayed on the output plots
%                   tifFrameBounds: (1x2 double array, NOT REQUIRED) start and stop frames to process of the tif, default is 2nd frame to the end to avoid occasional distortion in first image
%
% OUTPUT:
%
% FUNCTIONS USED:
%
% LIBARIES USED:
%
% NOTES:
%
% WRITTEN BY:       Spencer Garborg 1/22/19
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function motionTrackingStationaryAccuracy2P(tifFileName,redoFilterTF,updateSearchTF,medFiltTF,saveOutputVideoTF,threeStepTF,compiledTifTF,binSize,targetAvgNum,framesPerSecond,objMag,digMag,turnabout,commentString,tifFrameBounds)
%% Initialization
close all;

% Input video file which needs to be stabilized.
if compiledTifTF
    [imStack,tifLength] = imread_big(tifFileName);
else
    tifLength = length(imfinfo(tifFileName));
end

if exist('tifFrameBounds','var')
    aviFileName = imageFilter2P(tifFileName,redoFilterTF,medFiltTF,compiledTifTF,tifFrameBounds);
else
    aviFileName = imageFilter2P(tifFileName,redoFilterTF,medFiltTF,compiledTifTF);
    tifFrameBounds = [2 tifLength];
end

% % Get calibration values
% if exist('C:\Workspace\Code\DrewLab\calibrationValues.mat','file')
%     load('C:\Workspace\Code\DrewLab\calibrationValues.mat');
% else
%     disp('No calibration file found')
%     return
% end
% fns = fieldnames(calibrationValues);
% if ~isempty(fns)
%     if ~any(strcmp(['file_' calibrationFileString],fns))
%         disp('No matching calibration values found')
%         return
%     end
%     calibrationFileString = ['file_' calibrationFileString];
% else
%     disp('No calibration values found');
%     return
% end

% Get file name
[tokens,~] = regexpi(tifFileName,'\\([^\\]*).TIF','tokens','match');
fileName = tokens{1}{1};

% Create a System object(TM) to read video from a multimedia file. We set the
% output to be of intensity only video.
hVideoSource = vision.VideoFileReader(aviFileName, 'ImageColorSpace', 'Intensity', 'VideoOutputDataType', 'double');

%% Create template
% Create a template matcher System object to compute the location of the
% best match of the target in the video frame. We use this location to find
% translation between successive video frames.
if threeStepTF
    hTM = vision.TemplateMatcher('ROIInputPort', true, 'BestMatchNeighborhoodOutputPort', true, 'SearchMethod', 'Three-step');
else
    hTM = vision.TemplateMatcher('ROIInputPort', true, 'BestMatchNeighborhoodOutputPort', true);
end
                        
%% Create motion tracking output display
% Create a System object to display the original video and the stabilized
% video.
hVideoOut = vision.VideoPlayer('Name', 'Video Stabilization');
hVideoOut.Position(1) = round(0.4*hVideoOut.Position(1));
hVideoOut.Position(2) = round(.5*(hVideoOut.Position(2)));
hVideoOut.Position(3:4) = [1050 550];

%% Initialize variables for processing loop
% Get target window
if compiledTifTF
    if medFiltTF
        initialImage = medfilt2(imStack(:,:,tifFrameBounds(1)));
    else
        initialImage = imStack(:,:,tifFrameBounds(1));
    end
else
    if medFiltTF
        initialImage = medfilt2(imread(tifFileName,tifFrameBounds(1)));
    else
        initialImage = imread(tifFileName,tifFrameBounds(1));
    end
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
targetNum = 0;
targetSum = zeros(length(TargetRowIndices),length(TargetColIndices));
moveDist = [];
velocity = [];
% surfaceCalibFitX = calibrationValues.(calibrationFileString).surfaceCalibFitX;
% surfaceCalibFitY = calibrationValues.(calibrationFileString).surfaceCalibFitY;
tifLength = tifFrameBounds(2)-tifFrameBounds(1)+1;
imageStack = cell(1,tifLength);
midlineX = round(W/2);
midlineY = round(H/2);
secondsPerFrame = 1/framesPerSecond;
if ~exist('commentString','var') || isempty(commentString)
    commentString = 'No comments';
end

%% Stream Processing Loop
% This is the main processing loop which uses the objects we instantiated
% above to stabilize the input video.
while ~isDone(hVideoSource)
    input = hVideoSource();

    % Find location of Target in the input video frame
    if firstTime
      Idx = int32(pos.template_center_pos);
      motionVector = [0 0];
    else
      IdxPrev = Idx;
% IdxPrev = int32(pos.template_center_pos);

      ROI = [SearchRegion, pos.template_size+2*pos.search_border];
      Idx = hTM(input,Target,ROI);
      
      motionVector = double(Idx-IdxPrev);
    end

    if updateSearchTF
        [Offset, SearchRegion] = updatesearch(sz, motionVector, ...
            SearchRegion, Offset, pos);
    else
        [Offset] = updatesearch(sz, motionVector, ...
            SearchRegion, Offset, pos);
    end

    % Translate video frame to offset the camera motion
    Stabilized = imtranslate(input, Offset, 'linear');
    
    targetNum = targetNum + 1;
    if targetNum <= targetAvgNum
        targetSum = targetSum + Stabilized(round(TargetRowIndices), round(TargetColIndices));
        Target = targetSum ./ targetNum;
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
    txt = sprintf('(%+05.1f pixels,%+05.1f pixels)', Offset);
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
    if firstTime
        targetPositionPixel = [double(Idx(1)), H-double(Idx(2))+1];
        targetPosition = [0 0];
        firstTime = false;
    else
        motionVector(2) = -motionVector(2);
        targetPositionPixel = [targetPositionPixel;targetPositionPixel(end,:)+motionVector];
    end
end

% movementLength = size(targetPositionPixel,1);
% f = waitbar(0,'Calculating position from calibration');
% for n = 2:movementLength
%     % Get x movement
%     micronDistTraveledX = getDistFromImCentX(targetPositionPixel(n-1,1),targetPositionPixel(n-1,2),targetPositionPixel(n,1),targetPositionPixel(n,2),midlineX,surfaceCalibFitX);
%     
%     % Get y movement
%     micronDistTraveledY = getDistFromImCentY(targetPositionPixel(n-1,1),targetPositionPixel(n-1,2),targetPositionPixel(n,1),targetPositionPixel(n,2),midlineY,surfaceCalibFitY);
%     
%     moveDist(n-1,:) = [micronDistTraveledX,micronDistTraveledY];
%     velocity(n-1,1) = sqrt((micronDistTraveledX)^2+(micronDistTraveledY)^2)/secondsPerFrame;
%     targetPosition(n,:) = targetPosition(n-1,:)+[micronDistTraveledX,micronDistTraveledY];
%     waitbar(round((n-1)/(movementLength-1),2),f,'Calculating position from calibration');
% end
% close(f)

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

% Write data values to .mat file structure
movementData.fileName = fileName;
movementData.moveDist = moveDist;
movementData.velocity = velocity;
movementData.targetPosition = targetPosition;
movementData.frames = tifFrameBounds;
movementData.imageSize = sz;
movementData.medFiltTF = medFiltTF;
movementData.pos = pos;
movementData.ballData = ballData;
movementData.secondsPerFrame = secondsPerFrame;
movementData.objMag = objMag;
movementData.digMag = digMag;
movementData.turnabout = turnabout;
movementData.commentString = commentString;

matFileName = [tifFileName(1:end-4) '_processed.mat'];
save(matFileName,'movementData');

% Process data
binVec = 1:binSize:size(targetPositionPixel,1);
for n = 1:length(binVec)-1
    averageTargetPixelPosition(n,:) = mean(targetPositionPixel(binVec(n):binVec(n+1)-1,:));
end

% Plot data
subtitle = [num2str(1/secondsPerFrame) ' Frames/s, ' num2str(secondsPerFrame*(diff(tifFrameBounds)+1)) ' Seconds, ' num2str(objMag*digMag) 'x Magnification (' num2str(objMag) 'x Objective, ' num2str(digMag) 'x Digital), Turnabout = ' num2str(turnabout)];

h(1) = figure('Color','White');
k = convhull(targetPositionPixel(:,1),targetPositionPixel(:,2));
plot(targetPositionPixel(k,1),targetPositionPixel(k,2),'b',targetPositionPixel(:,1),targetPositionPixel(:,2),'k');
maxVal = ceil(max(max(abs(targetPositionPixel)))/10)*10;
axis equal square
% axis([-maxVal maxVal -maxVal maxVal])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title(['\fontsize{20pt}\bf{Position of Target Object}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixels')
ylabel('Pixels')

h(2) = figure('Color','White');
hist(targetPositionPixel(:,1),500);
title(['\fontsize{20pt}\bf{Position of Target Object Histogram (X)}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixels')
ylabel('Number of Data Points')

h(3) = figure('Color','White');
hist(targetPositionPixel(:,2),500);
title(['\fontsize{20pt}\bf{Position of Target Object Histogram (Y)}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixels')
ylabel('Number of Data Points')

h(4) = figure('Color','White');
k = convhull(averageTargetPixelPosition(:,1),averageTargetPixelPosition(:,2));
plot(averageTargetPixelPosition(k,1),averageTargetPixelPosition(k,2),'b',averageTargetPixelPosition(:,1),averageTargetPixelPosition(:,2),'k.');
maxVal = ceil(max(max(abs(averageTargetPixelPosition)))/10)*10;
axis equal square
% axis([-maxVal maxVal -maxVal maxVal])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title(['\fontsize{20pt}\bf{Position of Target Object - Binned}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixels')
ylabel('Pixels')

% Save figures to single .fig file
savefig(h,[matFileName(1:end-14) '_outputPlots.fig']);
end

