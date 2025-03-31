%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:    motionTracking2P3DRef
%
% FUNCTION:         motionTracking2P3DRef(tifFileName,calibrationFileString,redoFilterTF,updateSearchTF,medFiltTF,saveOutputVideoTF,threeStepTF,compiledTifTF,targetAvgNum,framesPerSecond,objMag,digMag,turnabout,commentString,tifFrameBounds)
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
%                   layer: (double) layer 1 or 2 designation
%                   hemisphere: (double) 1 = right, 2 = left
%                   commentString: (string, NOT REQUIRED) include any notes about this run to be displayed on the output plots
%                   tifFrameBounds: (1x2 double array, NOT REQUIRED) start and stop frames to process of the tif, default is 2nd frame to the end to avoid occasional distortion in first image
%
% OUTPUT:
%
% FUNCTIONS USED:
%
% LIBARIES USED:
%
% NOTES:            max(3dVar,[],1) = xz, max(3dVar,[],2) = yz, max(3dVar,[],3) = xy
%
% WRITTEN BY:       Spencer Garborg 6/4/20
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function motionTracking2P3DRef(tifFileName,referenceStackTifFileName,calibrationFileString,updateSearchTF,medFiltTF,saveOutputVideoTF,threeStepTF,compiledTifTF,targetAvgNum,framesPerSecond,objMag,digMag,turnabout,layers,micronsPerLayerRef,hemisphere,commentString,tifFrameBounds)
%% Initialization
close all;

% Input video file which needs to be stabilized.
% f = waitbar(0,'Compiling and Filtering Images');
if compiledTifTF
    [imStack,tifLength] = imread_big(tifFileName);
    if ~exist('tifFrameBounds','var')
        tifFrameBounds = [2 tifLength];
    end
else
    tifLength = length(imfinfo(tifFileName));
    if ~exist('tifFrameBounds','var')
        tifFrameBounds = [2 tifLength];
    end
%     imStack = [];
%     for n = 1:tifLength
%         waitbar(round((n-tifFrameBounds(1))/tifLength,2),f,'Compiling and Filtering Images');
%         imStack(:,:,end+1) = imread(tifFileName, n);
%     end
end
tifLengthRef = length(imfinfo(referenceStackTifFileName));
% zMicronVals = [0:tifLengthRef-1] * micronsPerLayerRef;
referenceStack = [];
for i = 1:tifLengthRef
    referenceStack(:,:,end+1) = im2double(imread(referenceStackTifFileName, i));
end
% close(f);

% if exist('tifFrameBounds','var')
%     aviFileName = imageFilter2P(tifFileName,redoFilterTF,medFiltTF,compiledTifTF,tifFrameBounds);
% else
%     aviFileName = imageFilter2P(tifFileName,redoFilterTF,medFiltTF,compiledTifTF);
%     tifFrameBounds = [2 tifLength];
% end

% Get calibration values
if exist('C:\Workspace\Code\DrewLab\calibrationValues.mat','file')
    load('C:\Workspace\Code\DrewLab\calibrationValues.mat');
else
    disp('No calibration file found')
    return
end
fns = fieldnames(calibrationValues);
if ~isempty(fns)
    if ~any(strcmp(['file_' calibrationFileString],fns))
        disp('No matching calibration values found')
        return
    end
    calibrationFileString = ['file_' calibrationFileString];
else
    disp('No calibration values found');
    return
end

% Get file name
[tokens,~] = regexpi(tifFileName,'\\([^\\]*).TIF','tokens','match');
fileName = tokens{1}{1};

% Create a System object(TM) to read video from a multimedia file. We set the
% output to be of intensity only video.
% hVideoSource = vision.VideoFileReader(aviFileName, 'ImageColorSpace', 'Intensity', 'VideoOutputDataType', 'double');

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
hVideoOut = vision.VideoPlayer('Name', 'Video Stabilization Single Target');
hVideoOut.Position(1) = round(.3*hVideoOut.Position(1));
hVideoOut.Position(2) = round(.5*(hVideoOut.Position(2)));
hVideoOut.Position(3:4) = [1050 550];

%% Initialize variables for processing loop
if medFiltTF
    initialImage = im2double(medfilt2(imread(tifFileName, tifFrameBounds(1))));
else
    initialImage = im2double(imread(tifFileName, tifFrameBounds(1)));
end

initialImageRef = [];
if medFiltTF
    for n = 1:length(imfinfo(referenceStackTifFileName))
        initialImageRef(:,:,n) = im2double(medfilt2(imread(referenceStackTifFileName, n)));
    end
else
    for n = 1:length(imfinfo(referenceStackTifFileName))
        initialImageRef(:,:,n) = im2double(imread(referenceStackTifFileName, n));
    end
end

% Get target window
imshow(initialImage);
title('Select upper left, then lower right target corners and press enter');
[inputCoordTargetX,inputCoordTargetY] = getpts(gcf);
close(gcf);
inputCoordTargetX = round(inputCoordTargetX);
inputCoordTargetY = round(inputCoordTargetY);
pos.template_orig = [inputCoordTargetX(1) inputCoordTargetY(1)]; % [x y] upper left corner
pos.template_size = [inputCoordTargetX(2) inputCoordTargetY(2)] - [inputCoordTargetX(1) inputCoordTargetY(1)];   % [width height]

% Find match for first target window
TargetRowIndicesMovie = pos.template_orig(2)-1:pos.template_orig(2)+pos.template_size(2)-2;
TargetColIndicesMovie = pos.template_orig(1)-1:pos.template_orig(1)+pos.template_size(1)-2;
Target = initialImage(round(TargetRowIndicesMovie), round(TargetColIndicesMovie));
for i = 1:size(initialImageRef,3)
    corrVals(:,:,i) = normxcorr2(Target,initialImageRef(:,:,i));
end
[~,idx] = max(corrVals(:));
[r,c,p] = ind2sub(size(corrVals),idx);
TargetRowIndices = r-length(TargetRowIndicesMovie)+1:r;
TargetColIndices = c-length(TargetColIndicesMovie)+1:c;
pos.template_orig = [c-length(TargetColIndicesMovie)+2 r-length(TargetRowIndicesMovie)+2]; % [x y] upper left corner
pos.template_size = [c+2 r+2] - [c-length(TargetColIndicesMovie)+2 r-length(TargetRowIndicesMovie)+2];   % [width height]

% Get search window
imshow(initialImageRef(:,:,p));
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
W = size(initialImage,2); % Width of video in pixels
H = size(initialImage,1); % Height of video in pixels
sz = [W,H];
SearchRegion = pos.template_orig - pos.search_border - 1;
BorderCols = [1:pos.search_border(1)+4 W-pos.search_border(1)+4:W];
BorderRows = [1:pos.search_border(2)+4 H-pos.search_border(2)+4:H];
Offset = [0 0];
Target = zeros(length(TargetRowIndices),length(TargetColIndices));
firstTime = true;
n = 1;
targetNum = 0;
targetSum = zeros(length(TargetRowIndices),length(TargetColIndices));
moveDist = [];
velocity = [];
surfaceCalibFitX = calibrationValues.(calibrationFileString).surfaceCalibFitX;
surfaceCalibFitY = calibrationValues.(calibrationFileString).surfaceCalibFitY;
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
for i = tifFrameBounds(1):tifFrameBounds(2)
    if medFiltTF
        inputMovie = im2double(medfilt2(imread(tifFileName, i)));
    else
        inputMovie = im2double(imread(tifFileName, i));
    end
    
    Target = inputMovie(round(TargetRowIndicesMovie), round(TargetColIndicesMovie));
    for j = 1:size(initialImageRef,3)
        corrVals(:,:,j) = normxcorr2(Target,initialImageRef(:,:,j));
    end
    [~,idx] = max(corrVals(:));
    [r,c,p] = ind2sub(size(corrVals),idx);
    
    input = initialImageRef(:,:,p);
    
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
%     Stabilized = imtranslate(input, Offset, 'linear');
    
%     targetNum = targetNum + 1;
%     if targetNum <= targetAvgNum
%         targetSum = targetSum + Stabilized(round(TargetRowIndices), round(TargetColIndices));
%         Target = targetSum ./ targetNum;
%         %         imshow(Target);
%     end

    % Add black border for display
    Stabilized(:, BorderCols) = 0;
    Stabilized(BorderRows, :) = 0;

    TargetRect = [pos.template_orig-Offset, pos.template_size];
    SearchRegionRect = [SearchRegion, pos.template_size + 2*pos.search_border];
    TargetRectMovie = [TargetColIndicesMovie(1) TargetRowIndicesMovie(1) length(TargetColIndicesMovie) length(TargetRowIndicesMovie)];

    % Draw rectangles on input to show target and search region
    input = insertShape(input, 'Rectangle', [TargetRect; SearchRegionRect],...
                        'Color', 'white');
    inputMovie = insertShape(inputMovie, 'Rectangle', TargetRectMovie,...
                        'Color', 'white');
    % Display the offset (displacement) values on the input image
    if firstTime
        Offset3 = [Offset 0];
    else
        Offset3 = [Offset p - targetPositionPixel(1,3)];
    end
    txt = sprintf('(%+05.1f pixels,%+05.1f pixels,%+05.1f pixels)', Offset3);
    input = insertText(input(:,:,1),[1 1],txt,'FontSize',16, ...
                    'TextColor', 'white', 'BoxOpacity', 0);
    
    % Save output video to variable
    if saveOutputVideoTF
        imageStack{1,n} = [input(:,:,1) Stabilized];      
        n = n + 1;
    end
    
    % Add pixel motion to data
    if firstTime
        targetPositionPixel = [double(Idx(1)), H-double(Idx(2))+1, p];
        targetPosition = [0 0 0];
        firstTime = false;
    else
        motionVector(2) = -motionVector(2);
        motionVector(3) = targetPositionPixel(end,3) - p;
%         targetPositionPixel = [targetPositionPixel;targetPositionPixel(end,:)+motionVector];
        targetPositionPixel = [targetPositionPixel;targetPositionPixel(end,:)-motionVector];
    end
    
    %%%%%%%%%%% Display video %%%%%%%%%%%%%
    hVideoOut([input(:,:,1) inputMovie(:,:,1)]);
end

% Perform movement calculations
movementLength = size(targetPositionPixel,1);
f = waitbar(0,'Calculating position from calibration');
for n = 2:movementLength
    %%%%%%%%%%%% XY %%%%%%%%%%%%%%
    % Get x movement
    micronDistTraveledX = getDistFromImCentX(targetPositionPixel(n-1,1),targetPositionPixel(n-1,2),targetPositionPixel(n,1),targetPositionPixel(n,2),midlineX,surfaceCalibFitX);
    
    % Get y movement
    micronDistTraveledY = getDistFromImCentY(targetPositionPixel(n-1,1),targetPositionPixel(n-1,2),targetPositionPixel(n,1),targetPositionPixel(n,2),midlineY,surfaceCalibFitY);
    
    % Get z movement
    micronDistTraveledZ = (targetPositionPixel(n,3) - targetPositionPixel(n-1,3)) * micronsPerLayerRef;
    
    moveDist(n-1,:) = [micronDistTraveledX,micronDistTraveledY,micronDistTraveledZ];
    velocity(n-1,1) = sqrt((micronDistTraveledX)^2+(micronDistTraveledY)^2)/(layers*secondsPerFrame);
    targetPosition(n,:) = targetPosition(n-1,:)+[micronDistTraveledX,micronDistTraveledY,micronDistTraveledZ];
    
    waitbar(round((n-1)/(movementLength-1),2),f,'Calculating position from calibration');
end
close(f)
% meanPosX = mean(targetPosition(:,1));
% meanPosY = mean(targetPosition(:,2));
% targetPosition(:,1) = targetPosition(:,1) - meanPosX;
% targetPosition(:,2) = targetPosition(:,2) - meanPosY;
% 
% meanPosX = mean(targetPositionXZ(:,1));
% meanPosY = mean(targetPositionXZ(:,2));
% targetPositionXZ(:,1) = targetPositionXZ(:,1) - meanPosX;
% targetPositionXZ(:,2) = targetPositionXZ(:,2) - meanPosY;
% 
% meanPosX = mean(targetPositionYZ(:,1));
% meanPosY = mean(targetPositionYZ(:,2));
% targetPositionYZ(:,1) = targetPositionYZ(:,1) - meanPosX;
% targetPositionYZ(:,2) = targetPositionYZ(:,2) - meanPosY;

%% Release
% Here you call the release method on the objects to close any open files
% and devices.
% release(hVideoSource);

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
movementData.calibrationFileString = calibrationFileString;
movementData.frames = tifFrameBounds;
movementData.imageSize = sz;
movementData.medFiltTF = medFiltTF;
movementData.pos = pos;
movementData.ballData = ballData;
movementData.secondsPerFrame = secondsPerFrame*layers;
movementData.objMag = objMag;
movementData.digMag = digMag;
movementData.turnabout = turnabout;
movementData.commentString = commentString;
movementData.layers = layers;
movementData.micronsPerLayerRef = micronsPerLayerRef;
movementData.hemisphere = hemisphere;

matFileName = [tifFileName(1:end-4) '3DRef.mat'];
save(matFileName,'movementData');

% Plot data
plotMotionTracking3DRef(matFileName);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subfunctions
function micronDistTraveledX = getDistFromImCentX(prevPixelLocX,prevPixelLocY,currPixelLocX,currPixelLocY,midlinePixelVal,calibSurf)

% Get distance between midline and previous x value at certain y value from surface calibration
prevCalibVec = [];
if prevPixelLocX == midlinePixelVal % Dist is zero
    prevMicronDistFromMid = 0;
elseif prevPixelLocX < midlinePixelVal % Dist is negative from center of image
    for x = prevPixelLocX:.5:midlinePixelVal
        prevCalibVec(end+1,1) = calibSurf(x,prevPixelLocY);
    end
    prevMicronDistFromMid = -trapz(prevPixelLocX:.5:midlinePixelVal,prevCalibVec);
else % Dist is positive from center of image
    for x = midlinePixelVal:.5:prevPixelLocX
        prevCalibVec(end+1,1) = calibSurf(x,prevPixelLocY);
    end
    prevMicronDistFromMid = trapz(midlinePixelVal:.5:prevPixelLocX,prevCalibVec);
end

% Get distance between midline and current x value at certain y value from surface calibration
currCalibVec = [];
if currPixelLocX == midlinePixelVal % Dist is zero
    currMicronDistFromMid = 0;
elseif currPixelLocX < midlinePixelVal % Dist is negative from center of image
    for x = currPixelLocX:.5:midlinePixelVal
        currCalibVec(end+1,1) = calibSurf(x,currPixelLocY);
    end
    currMicronDistFromMid = -trapz(currPixelLocX:.5:midlinePixelVal,currCalibVec);
else % Dist is positive from center of image
    for x = midlinePixelVal:.5:currPixelLocX
        currCalibVec(end+1,1) = calibSurf(x,currPixelLocY);
    end
    currMicronDistFromMid = trapz(midlinePixelVal:.5:currPixelLocX,currCalibVec);
end

% Distance travelled in x is difference between two distances from midline
micronDistTraveledX = currMicronDistFromMid - prevMicronDistFromMid;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function micronDistTraveledY = getDistFromImCentY(prevPixelLocX,prevPixelLocY,currPixelLocX,currPixelLocY,midlinePixelVal,calibSurf)

% Get distance between midline and previous x value at certain y value from surface calibration
prevCalibVec = [];
if prevPixelLocY == midlinePixelVal % Dist is zero
    prevMicronDistFromMid = 0;
elseif prevPixelLocY < midlinePixelVal % Dist is negative from center of image
    for y = prevPixelLocY:.5:midlinePixelVal
        prevCalibVec(end+1,1) = calibSurf(prevPixelLocX,y);
    end
    prevMicronDistFromMid = -trapz(prevPixelLocY:.5:midlinePixelVal,prevCalibVec);
else % Dist is positive from center of image
    for y = midlinePixelVal:.5:prevPixelLocY
        prevCalibVec(end+1,1) = calibSurf(prevPixelLocX,y);
    end
    prevMicronDistFromMid = trapz(midlinePixelVal:.5:prevPixelLocY,prevCalibVec);
end

% Get distance between midline and current x value at certain y value from surface calibration
currCalibVec = [];
if currPixelLocY == midlinePixelVal % Dist is zero
    currMicronDistFromMid = 0;
elseif currPixelLocY < midlinePixelVal % Dist is negative from center of image
    for y = currPixelLocY:.5:midlinePixelVal
        currCalibVec(end+1,1) = calibSurf(currPixelLocX,y);
    end
    currMicronDistFromMid = -trapz(currPixelLocY:.5:midlinePixelVal,currCalibVec);
else % Dist is positive from center of image
    for y = midlinePixelVal:.5:currPixelLocY
        currCalibVec(end+1,1) = calibSurf(currPixelLocX,y);
    end
    currMicronDistFromMid = trapz(midlinePixelVal:.5:currPixelLocY,currCalibVec);
end

% Distance travelled in x is difference between two distances from midline
micronDistTraveledY = currMicronDistFromMid - prevMicronDistFromMid;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function micronDistTraveledX = getDistFromImCentX_XZ(prevPixelLocX,currPixelLocX,midlinePixelVal,meanXCalibVec)
% Get distance between midline and previous x value at certain y value from surface calibration
if prevPixelLocX == midlinePixelVal % Dist is zero
    prevMicronDistFromMid = 0;
elseif prevPixelLocX < midlinePixelVal % Dist is negative from center of image
    prevMicronDistFromMid = -trapz(prevPixelLocX:.5:midlinePixelVal,meanXCalibVec((2*prevPixelLocX)-1:(2*midlinePixelVal)-1));
else % Dist is positive from center of image
    prevMicronDistFromMid = trapz(midlinePixelVal:.5:prevPixelLocX,meanXCalibVec((2*midlinePixelVal)-1:(2*prevPixelLocX)-1));
end

% Get distance between midline and current x value at certain y value from surface calibration
if currPixelLocX == midlinePixelVal % Dist is zero
    currMicronDistFromMid = 0;
elseif currPixelLocX < midlinePixelVal % Dist is negative from center of image
    currMicronDistFromMid = -trapz(currPixelLocX:.5:midlinePixelVal,meanXCalibVec((2*currPixelLocX)-1:(2*midlinePixelVal)-1));
else % Dist is positive from center of image
    currMicronDistFromMid = trapz(midlinePixelVal:.5:currPixelLocX,meanXCalibVec((2*midlinePixelVal)-1:(2*currPixelLocX)-1));
end

% Distance travelled in x is difference between two distances from midline
micronDistTraveledX = currMicronDistFromMid - prevMicronDistFromMid;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function micronDistTraveledY = getDistFromImCentY_YZ(prevPixelLocY,currPixelLocY,midlinePixelVal,meanYCalibVec)
% Get distance between midline and previous x value at certain y value from surface calibration
if prevPixelLocY == midlinePixelVal % Dist is zero
    prevMicronDistFromMid = 0;
elseif prevPixelLocY < midlinePixelVal % Dist is negative from center of image
    prevMicronDistFromMid = -trapz(prevPixelLocY:.5:midlinePixelVal,meanYCalibVec((2*prevPixelLocY)-1:(2*midlinePixelVal)-1));
else % Dist is positive from center of image
    prevMicronDistFromMid = trapz(midlinePixelVal:.5:prevPixelLocY,meanYCalibVec((2*midlinePixelVal)-1:(2*prevPixelLocY)-1));
end

% Get distance between midline and current x value at certain y value from surface calibration
if currPixelLocY == midlinePixelVal % Dist is zero
    currMicronDistFromMid = 0;
elseif currPixelLocY < midlinePixelVal % Dist is negative from center of image
    currMicronDistFromMid = -trapz(currPixelLocY:.5:midlinePixelVal,meanYCalibVec((2*currPixelLocY)-1:(2*midlinePixelVal)-1));
else % Dist is positive from center of image
    currMicronDistFromMid = trapz(midlinePixelVal:.5:currPixelLocY,meanYCalibVec((2*midlinePixelVal)-1:(2*currPixelLocY)-1));
end

% Distance travelled in x is difference between two distances from midline
micronDistTraveledY = currMicronDistFromMid - prevMicronDistFromMid;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [meanXCalibVec,meanYCalibVec] = getMeanCalibVecs(width,height,calibSurfX,calibSurfY)

meanXCalibVec = [];
for n = 1:.5:width
    calibValSum = 0;
    numVals = 0;
    for i = ceil(.25*height):2:floor(height-(.25*height))
        calibValSum = calibValSum + calibSurfX(n,i);
        numVals = numVals + 1;
    end
    meanXCalibVec(end+1) = calibValSum/numVals;
end

meanYCalibVec = [];
for n = 1:.5:height
    calibValSum = 0;
    numVals = 0;
    for i = ceil(.25*width):2:floor(width-(.25*width))
        calibValSum = calibValSum + calibSurfY(n,i);
        numVals = numVals + 1;
    end
    meanYCalibVec(end+1) = calibValSum/numVals;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mip = MIPXY(Mat3D)

mip = max(Mat3D,[],3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mip = MIPXZ(Mat3D)

maxInt3D = max(Mat3D,[],1);
H = size(maxInt3D,3);
for n = 1:H
    mip(n,:) = maxInt3D(:,:,H+1-n);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mip = MIPYZ(Mat3D)

maxInt3D = max(Mat3D,[],2);
H = size(maxInt3D,3);
for n = 1:H
    mip(n,:) = maxInt3D(:,:,H+1-n)';
end
end