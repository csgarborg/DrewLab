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
function motionTracking2P(tifFileName,calibrationFileString,updateSearchTF,medFiltTF,saveOutputVideoTF,threeStepTF,compiledTifTF,tempMedFiltTF,targetAvgNum,framesPerSecond,analogSampleRate,objMag,digMag,turnabout,hemisphere,commentString,tifFrameBounds)
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
hVideoOut = vision.VideoPlayer('Name', 'Video Stabilization');
hVideoOut.Position(1) = round(0.4*hVideoOut.Position(1));
hVideoOut.Position(2) = round(.5*(hVideoOut.Position(2)));
hVideoOut.Position(3:4) = [1050 550];

%% Initialize variables for processing loop
% Get target window
if compiledTifTF
    if medFiltTF
        initialImage = im2double(medfilt2(imStack(:,:,tifFrameBounds(1))));
    else
        initialImage = im2double(imStack(:,:,tifFrameBounds(1)));
    end
else
    if medFiltTF
        initialImage = im2double(medfilt2(imread(tifFileName, tifFrameBounds(1))));
    else
        initialImage = im2double(imread(tifFileName, tifFrameBounds(1)));
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
W = size(initialImage,2); % Width of video in pixels
H = size(initialImage,1); % Height of video in pixels
sz = [W,H];
BorderCols = [1:pos.search_border(1)+4 W-pos.search_border(1)+4:W];
BorderRows = [1:pos.search_border(2)+4 H-pos.search_border(2)+4:H];
TargetRowIndices = ...
  pos.template_orig(2)-1:pos.template_orig(2)+pos.template_size(2)-2;
TargetColIndices = ...
  pos.template_orig(1)-1:pos.template_orig(1)+pos.template_size(1)-2;
SearchRegion = pos.template_orig - pos.search_border - 1;
Offset = [0 0];
Target = zeros(length(TargetRowIndices),length(TargetColIndices));
firstTime = true;
nImg = 1;
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
if tempMedFiltTF
    f = waitbar(0,'Loading frames for temporal median filtering');    
    inputStack = [];
    for n = tifFrameBounds(1):tifFrameBounds(2)
        if compiledTifTF
            if medFiltTF
                inputStack(:,:,end+1) = im2double(medfilt2(imStack(:,:,n)));
            else
                inputStack(:,:,end+1) = im2double(imStack(:,:,n));
            end
        else
            if medFiltTF
                inputStack(:,:,end+1) = im2double(medfilt2(imread(tifFileName, n)));
            else
                inputStack(:,:,end+1) = im2double(imread(tifFileName, n));
            end
        end
        waitbar(round((n-tifFrameBounds(1)+1)/(tifFrameBounds(2)-tifFrameBounds(1)+1),2),f,'Loading frames for temporal median filtering');
    end
    close(f)
    f = waitbar(0,'Applying temporal median filter');
    for i = 1:size(inputStack,1)
        for j = 1:size(inputStack,2)
            inputStack(i,j,:) = medfilt1(inputStack(i,j,:),3);
            waitbar(round((((i-1)*size(inputStack,2))+j)/(size(inputStack,1)*size(inputStack,2)),2),f,'Applying temporal median filter');
        end
    end
    close(f)
end

%% Stream Processing Loop
% This is the main processing loop which uses the objects we instantiated
% above to stabilize the input video.
for i = tifFrameBounds(1)+2:tifFrameBounds(2)-2
    if tempMedFiltTF
        input = inputStack(:,:,i-tifFrameBounds(1)+1);
    else
        if compiledTifTF
            if medFiltTF
                input = im2double(medfilt2(imStack(:,:,i)));
            else
                input = im2double(imStack(:,:,i));
            end
        else
            if medFiltTF
                input = im2double(medfilt2(imread(tifFileName, i)));
            else
                input = im2double(imread(tifFileName, i));
            end
        end
    end

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
        imshow(Target);
    end
    
    % Add black border for display
    Stabilized(:, BorderCols) = 0;
    Stabilized(BorderRows, :) = 0;
    
    TargetRect = [pos.template_orig-Offset, pos.template_size];
    SearchRegionRect = [SearchRegion, pos.template_size + 2*pos.search_border];
    
    % Draw rectangles on input to show target and search region
%         input = insertShape(input, 'Rectangle', [TargetRect; SearchRegionRect],...
%                             'Color', 'White');
    input = insertShape(input, 'Rectangle', SearchRegionRect,...
        'Color', [58 58 58], 'LineWidth', 2);
    input = insertShape(input, 'Rectangle', TargetRect,...
        'Color', 'white');
    % Display the offset (displacement) values on the input image
%     txt = sprintf('(%+05.1f pixels,%+05.1f pixels)', Offset);
%     input = insertText(input(:,:,1),[1 1],txt,'FontSize',16, ...
%         'TextColor', 'white', 'BoxOpacity', 0);
    % Display video
    hVideoOut([input(:,:,1) Stabilized]);
    if firstTime
        savedImage = input(:,:,1);
    end
    
    % Save output video to variable
    if saveOutputVideoTF
%         imageStack{1,n} = [input(:,:,1) Stabilized];
        imageStack{1,nImg} = input(:,:,1);
        nImg = nImg + 1;
    end
    
    % Add pixel motion to data
    if firstTime
        targetPositionPixel = [double(Idx(1)), 512-double(Idx(2))+1];
        targetPosition = [0 0];
        firstTime = false;
    else
        motionVector(2) = -motionVector(2);
        targetPositionPixel = [targetPositionPixel;targetPositionPixel(end,:)+motionVector];
    end
end

movementLength = size(targetPositionPixel,1);
uniquePositions = unique(targetPositionPixel,'rows');
lookupTableX = zeros(512,512);
lookupTableY = zeros(512,512);
f = waitbar(0,'Calculating position look-up table from calibration');
for n = 1:size(uniquePositions,1)
    lookupTableX(uniquePositions(n,1),uniquePositions(n,2)) = getDistFromImCentX(uniquePositions(n,1),uniquePositions(n,2),midlineX,surfaceCalibFitX);
    lookupTableY(uniquePositions(n,1),uniquePositions(n,2)) = getDistFromImCentY(uniquePositions(n,1),uniquePositions(n,2),512-midlineY,surfaceCalibFitY);
    waitbar(round((n)/(size(uniquePositions,1)),2),f,'Calculating position look-up table from calibration');
end
close(f)
f = waitbar(0,'Calculating position from look-up table');
for n = 2:movementLength
    % Get x movement
    micronDistTraveledX = lookupTableX(targetPositionPixel(n,1),targetPositionPixel(n,2)) - lookupTableX(targetPositionPixel(n-1,1),targetPositionPixel(n-1,2));
    
    % Get y movement
    micronDistTraveledY = lookupTableY(targetPositionPixel(n,1),targetPositionPixel(n,2)) - lookupTableY(targetPositionPixel(n-1,1),targetPositionPixel(n-1,2));
    
    moveDist(n-1,:) = [micronDistTraveledX,micronDistTraveledY];
    velocity(n-1,1) = sqrt((micronDistTraveledX)^2+(micronDistTraveledY)^2)/secondsPerFrame;
    targetPosition(n,:) = targetPosition(n-1,:)+[micronDistTraveledX,micronDistTraveledY];
    waitbar(round((n-1)/(movementLength-1),2),f,'Calculating position from look-up table');
end
close(f)
% meanPosX = mean(targetPosition(:,1));
% meanPosY = mean(targetPosition(:,2));
% targetPosition(:,1) = targetPosition(:,1) - meanPosX;
% targetPosition(:,2) = targetPosition(:,2) - meanPosY;

%% Release
% Here you call the release method on the objects to close any open files
% and devices.
% release(hVideoSource);

%% Output data

% Save output video to AVI file
if saveOutputVideoTF
%     [HOut,WOut] = size(imageStack{1,1});
%     imageStackMat = cell2mat(imageStack);
%     imageStackMat = double(reshape(imageStackMat,HOut,WOut,tifLength));
%     aviFileNameOutput = [tifFileName(1:end-4) '_output.avi'];
%     if exist(aviFileNameOutput,'file')
%         delete(aviFileNameOutput)
%     end
%     aviObject = VideoWriter(aviFileNameOutput,'Uncompressed AVI');
%     aviObject.FrameRate = 10;
%     open(aviObject);
%     f = waitbar(0,'Creating output AVI file');
%     for k = 1:tifLength
%         waitbar(round(k/tifLength,2),f,'Creating output AVI file');
%         %     imagesc(filteredData(:,:,k));
%         %     filteredFrame = getframe(fig,[0.05 0.05 0.9 0.9]);
%         writeVideo(aviObject,imageStackMat(:,:,k));
%     end
%     close(f)
%     close(aviObject);
    tiffFileNameOutput = [tifFileName(1:end-4) '_output.tif'];
    imwrite(imageStack{1,1},tiffFileNameOutput)
    for n = 1:length(imageStack)
        if isempty(imageStack{1,n})
            continue
        end
        imwrite(imageStack{1,n},tiffFileNameOutput,'WriteMode','append')
    end
end

if exist([tifFileName(1:end-4) '_processed_1.mat'],'file')
    load([tifFileName(1:end-4) '_processed_1.mat']);
    procBallData = movementData.ballData;
    procEMGData = movementData.emgData;
    procEKGData = movementData.ekgData;
    motionEvents = movementData.motionEvents;
    EMGEvents = movementData.EMGEvents;
    EMGNoMotionEvents = movementData.EMGNoMotionEvents;
    clear movementData
elseif ~exist([tifFileName(1:end-4) '.txt'],'file')
    procBallData(:,1) = 0:1/30:(tifFrameBounds(2)-tifFrameBounds(1))*secondsPerFrame';
    procBallData(:,2) = zeros(length(procBallData),1);
    procEMGData = [procBallData(:,1) zeros(size(procBallData,1),1)];
    procEKGData = [procBallData(:,1) zeros(size(procBallData,1),1)];
    motionEvents = [];
    EMGEvents = [];
    EMGNoMotionEvents = [];
else
    % Get binary ball and EMG data to compare to frame movement data
    secondsBounds = [tifFrameBounds(1)*secondsPerFrame tifFrameBounds(2)*secondsPerFrame];
    ballData = load([tifFileName(1:end-4) '.txt']);
    ballDataIndex = secondsBounds(1)<=ballData(:,1) & ballData(:,1)<= secondsBounds(2);
    if size(ballData,2) > 2
        ballDataOnly = [ballData(ballDataIndex,1) ballData(ballDataIndex,2)];
        emgDataOnly = [ballData(ballDataIndex,1) ballData(ballDataIndex,3)];
        
        % procBallData = filterEMGData(ballDataOnly,analogSampleRate);
        procEMGData = filterEMGData(emgDataOnly,analogSampleRate);
        procEKGData = filterEKGData(emgDataOnly,analogSampleRate,[tifFileName(1:end-4) '_EKG.fig']);
        procBallData = smoothBallData(ballDataOnly,analogSampleRate);
    else
        procBallData = smoothBallData([ballData(ballDataIndex,1) ballData(ballDataIndex,2)],analogSampleRate);
        procEMGData = [procBallData(:,1) zeros(size(procBallData,1),1)];
        procEKGData = [procBallData(:,1) zeros(size(procBallData,1),1)];
    end
    
    % Detect motion and events
    [motionEvents,EMGEvents,EMGNoMotionEvents] = detectEvents(procBallData,procEMGData,analogSampleRate,secondsPerFrame);
end

% Write data values to .mat file structure
movementData.fileName = fileName;
movementData.moveDistRaw = moveDist;
movementData.velocityRaw = velocity;
movementData.targetPositionRaw = targetPosition;
movementData.moveDist = medfilt1(moveDist,6);
movementData.velocity = medfilt1(velocity,6);
movementData.targetPosition = medfilt1(targetPosition,6);
movementData.calibrationFileString = calibrationFileString;
movementData.frames = tifFrameBounds;
movementData.imageSize = sz;
movementData.medFiltTF = medFiltTF;
movementData.pos = pos;
movementData.ballData = procBallData;
movementData.emgData = procEMGData;
movementData.ekgData = procEKGData;
movementData.motionEvents = motionEvents;
movementData.EMGEvents = EMGEvents;
movementData.EMGNoMotionEvents = EMGNoMotionEvents;
movementData.secondsPerFrame = secondsPerFrame;
movementData.objMag = objMag;
movementData.digMag = digMag;
movementData.turnabout = turnabout;
movementData.commentString = commentString;
movementData.hemisphere = hemisphere;
movementData.savedImage = savedImage;

for n = 1:100
    matFileName = [tifFileName(1:end-4) '_processed_' num2str(n) '.mat'];
    if ~exist(matFileName,'file')
        save(matFileName,'movementData');
        break
    end
end

% Plot data
plotMotionTracking(matFileName);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subfunctions
function micronDistFromMid = getDistFromImCentX(pixelLocX,pixelLocY,midlinePixelVal,calibSurf)

% Get distance between midline and x value at certain y value from surface calibration
calibVec = [];
if pixelLocX == midlinePixelVal % Dist is zero
    micronDistFromMid = 0;
elseif pixelLocX < midlinePixelVal % Dist is negative from center of image
    for x = pixelLocX:.5:midlinePixelVal
        calibVec(end+1,1) = calibSurf(x,pixelLocY);
    end
    micronDistFromMid = -trapz(pixelLocX:.5:midlinePixelVal,calibVec);
else % Dist is positive from center of image
    for x = midlinePixelVal:.5:pixelLocX
        calibVec(end+1,1) = calibSurf(x,pixelLocY);
    end
    micronDistFromMid = trapz(midlinePixelVal:.5:pixelLocX,calibVec);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function micronDistFromMid = getDistFromImCentY(pixelLocX,pixelLocY,midlinePixelVal,calibSurf)

% Get distance between midline and previous x value at certain y value from surface calibration
prevCalibVec = [];
if pixelLocY == midlinePixelVal % Dist is zero
    micronDistFromMid = 0;
elseif pixelLocY < midlinePixelVal % Dist is negative from center of image
    for y = pixelLocY:.5:midlinePixelVal
        prevCalibVec(end+1,1) = calibSurf(pixelLocX,y);
    end
    micronDistFromMid = -trapz(pixelLocY:.5:midlinePixelVal,prevCalibVec);
else % Dist is positive from center of image
    for y = midlinePixelVal:.5:pixelLocY
        prevCalibVec(end+1,1) = calibSurf(pixelLocX,y);
    end
    micronDistFromMid = trapz(midlinePixelVal:.5:pixelLocY,prevCalibVec);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function procData = filterEMGData(emg,sampleRate)
% 
% close all
% rawEMG = emg(:,2);
% t = emg(:,1);
% % process EMG data
% % fpass = [300,(sampleRate/2)-1];   % Hz
% fpass = [300,3000];   % Hz
% trialDuration_sec = t(end);   % read this variable in from your data in seconds
% analogSamplingRate = sampleRate;   % Hz - change if yours is different
% dsFs = 30;   % Hz - downsampled frequency
% analogExpectedLength = trialDuration_sec*analogSamplingRate;
% trimmedEMG = rawEMG(1:min(analogExpectedLength,length(rawEMG)));
% [z,p,k] = butter(3,fpass/(analogSamplingRate/2));
% [sos,g] = zp2sos(z,p,k);
% filtEMG = filtfilt(sos,g,trimmedEMG - mean(trimmedEMG));
% % kernelWidth = 0.5;
% kernelWidth = 0.005;
% smoothingKernel = gausswin(kernelWidth*analogSamplingRate)/sum(gausswin(kernelWidth*analogSamplingRate));
% % EMGPwr = log10(conv(filtEMG.^2,smoothingKernel,'same'));
% EMGPwr = conv(filtEMG.^2,smoothingKernel,'same');
% resampEMG = resample(EMGPwr,dsFs,analogSamplingRate);
% procEMG = resampEMG;   % save this as your final array
% procT = 0:1/dsFs:t(end);
% % procT = procT(1:end-1);
% if length(procT) > length(procEMG)
%     procT = procT(1:length(procEMG));
% elseif length(procT) < length(procEMG)
%     procEMG = procEMG(1:length(procT));
% end
% plot(procT,procEMG)
% title('Select start and stop time (t) to use to determine mean baseline')
% tVals = ginput(2);
% close
% i = tVals(1,1) <= procT & procT <= tVals(2,1);
% baseline = mean(procEMG(i));
% procEMG = procEMG - (baseline - 1);
% procData = [procT',procEMG];
% end