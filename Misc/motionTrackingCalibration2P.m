%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:    motionTrackingCalibration2P
%
% FUNCTION:         motionTrackingCalibration2P(tifFileName,redoFilterTF,updateSearchTF,medFiltTF,saveOutputVideoTF,threeStepTF,compiledTifTF,targetAvgNum,micronJumpVal,framesPerSecond,objMag,digMag,turnabout,commentString,tifFrameBounds)
%
% DESCRIPTION:      Processes a tif file and extracts movement data based
%                   on user inputs
%
% INPUT:            tif must be extracted from MDF file and have the same name as ball data text file (also extracted from MDF file as ANSI encoding (8-bit))
%
% VARIABLES:        tifFileName: (string) complete path to tif data file
%                   redoFilterTF: (logical) select false when data have already been filtered and output AVI created but want to reprocess without redoing the filter
%                   updateSearchTF: (logical) false if the searched ROI should remain in the same selected place for tracking or true if should update every frame based on target location
%                   medFiltTF: (logical) true if spatial median filter should be applied to the image data
%                   saveOutputVideoTF: (logical) true if tracking output video should be saved to AVI file
%                   threeStepTF: (logical) true if speed of tracking code is more important than accuracy - not every point of the cross-correlation is calculated and position of target is narrowed down quickly
%                   compiledTifTF: (logical) true if tif file is larger than 4 GB and had to be stitched together in ImageJ from multiple seperate tif files
%                   targetAvgNum: (double) number of targets to average together at the beginning of the tif to get final target used to match and track object
%                   micronJumpVal: (double) size of steps in microns used in the calibration
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
% WRITTEN BY:       Spencer Garborg 1/23/19
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function motionTrackingCalibration2P(tifFileName,redoFilterTF,updateSearchTF,medFiltTF,saveOutputVideoTF,threeStepTF,compiledTifTF,targetAvgNum,micronJumpVal,framesPerSecond,objMag,digMag,turnabout,commentString,tifFrameBounds)
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
targetNum = 0;
targetSum = zeros(length(TargetRowIndices),length(TargetColIndices));
targetMiddleRow = round(length(TargetRowIndices)/2);
moveDist = [];
velocity = [];
targetArea = [];
fwhmTarget = [];
fwhmTargetPixels = [];
secondsPerFrame = 1/framesPerSecond;
tifLength = tifFrameBounds(2)-tifFrameBounds(1)+1;
imageStack = cell(1,tifLength);
targetStack = cell(1,tifLength);
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
        [Offset, SearchRegion] = updatesearch(sz, motionVector, SearchRegion, Offset, pos);
    else
        [Offset] = updatesearch(sz, motionVector, SearchRegion, Offset, pos);
    end

    % Translate video frame to offset the camera motion
    Stabilized = imtranslate(input, Offset, 'linear');
    
%     targetStack{1,n} = Stabilized(round(TargetRowIndices), round(TargetColIndices));
%     targetStack{1,n} = vertcat(zeros(20,size(targetStack{1,n},2)),targetStack{1,n});
%     if n == 1
%         origIm = targetStack{1,1};
%     end
%     targetStack{1,n} = horzcat(origIm,targetStack{1,n});
%     txt = sprintf('%+05.1fp(x),%+05.1fp(y)', double(Idx));
%     targetStack{1,n} = insertText(targetStack{1,n},[1 1],txt,'FontSize',8, 'TextColor', 'white', 'BoxOpacity', 0);
%     imshow(targetStack{1,n});
%     n = n + 1;

    try
        currTargetImageProcessed = imfill(bwareaopen(imbinarize(Stabilized(round(TargetRowIndices), round(TargetColIndices))),30),'holes');
        [~,currTargetImageProcessed] = bwboundaries(currTargetImageProcessed,'noholes');
        areaStruct = regionprops(currTargetImageProcessed,'Area');
        targetArea(end+1,1) = areaStruct(1).Area;
    catch
        targetArea(end+1,1) = -1;
    end
    
    fwhmTargetPixels(:,:,end+1) = mat2gray(Stabilized(TargetRowIndices(targetMiddleRow-2:targetMiddleRow+2), round(TargetColIndices)));
    
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
    txt = sprintf('(%+05.1f,%+05.1f)', Offset);
    input = insertText(input(:,:,1),[1 1],txt,'FontSize',16, 'TextColor', 'white', 'BoxOpacity', 0);
    
    % Display video
    hVideoOut([input(:,:,1) Stabilized]);
    
    % Save output video to variable
    if saveOutputVideoTF
        imageStack{1,n} = [input(:,:,1) Stabilized];
        n = n + 1;
    end
    
    % Add pixel motion to data
    if firstTime
        targetPosition = [double(Idx(1)), H-double(Idx(2))+1];
        firstTime = false;
    else
        motionVector(2) = -motionVector(2);
        moveDist = [moveDist;motionVector];
        velocity = [velocity;sqrt((motionVector(1))^2+(motionVector(2))^2)/secondsPerFrame];
        targetPosition = [targetPosition;targetPosition(end,:)+motionVector];
    end
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

% Generate plots
targetPosition(:,1) = medfilt1(targetPosition(:,1));
targetPosition(:,2) = medfilt1(targetPosition(:,2));
subtitle = [num2str(framesPerSecond) ' Frames/s, ' num2str(secondsPerFrame*(diff(tifFrameBounds)+1)) ' Seconds, ' num2str(micronJumpVal) ' ' char(181) 'm Jumps, ' num2str(objMag*digMag) 'x Magnification (' num2str(objMag) 'x Objective, ' num2str(digMag) 'x Digital), Turnabout = ' num2str(turnabout)];
h(1) = figure('Color','White');
subplot(3,1,1)
plot(1:size(moveDist,1),moveDist(:,1),'r')
title(['\fontsize{20pt}\bf{Object Movement Between Frames}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Frame')
ylabel('X Movement (Pixels)')
grid on
axis([1 size(moveDist,1) floor(min([moveDist(:,1);moveDist(:,2)])/10)*10 ceil(max([moveDist(:,1);moveDist(:,2)])/10)*10])
subplot(3,1,2)
plot(1:size(moveDist,1),moveDist(:,2),'b')
xlabel('Frame')
ylabel('Y Movement (Pixels)')
grid on
axis([1 size(moveDist,1) floor(min([moveDist(:,1);moveDist(:,2)])/10)*10 ceil(max([moveDist(:,1);moveDist(:,2)])/10)*10])
subplot(3,1,3)
plot(ballData(:,1),ballData(:,2),'k')
title('\fontsize{20pt}\bf{Ball Movement}')
xlabel('Time (s)')
ylabel('Movement')
grid on
axis([min(ballData(:,1)) max(ballData(:,1)) -1 ceil(max(ballData(:,2)))])

h(2) = figure('Color','White');
subplot(2,1,1)
plot(1:length(velocity),velocity,'r')
title(['\fontsize{20pt}\bf{Object Velocity Between Frames}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Frame')
ylabel('Velocity (Pixels/s)')
grid on
axis([1 size(velocity,1) 0 ceil(max(velocity(:,1))/10)*10])
subplot(2,1,2)
plot(ballData(:,1),ballData(:,2),'k')
title('\fontsize{20pt}\bf{Ball Movement}')
xlabel('Time (s)')
ylabel('Movement')
grid on
axis([min(ballData(:,1)) max(ballData(:,1)) -1 ceil(max(ballData(:,2)))])


frameSelection = [];
lowerBoundX = 1;
upperBoundX = 500;
cont = true;
h(3) = figure('Color','White','Name','Select frame starts and stops, select behind line to erase, right and left arrow key to move plot, press enter when finished','NumberTitle','off');
while cont
    if lowerBoundX <= length(targetPosition) && length(targetPosition) < upperBoundX
        targetPositionSegmentX = [targetPosition(lowerBoundX:length(targetPosition),1); zeros(upperBoundX-length(targetPosition),1)];
        targetPositionSegmentY = [targetPosition(lowerBoundX:length(targetPosition),2); zeros(upperBoundX-length(targetPosition),1)];
    elseif lowerBoundX > length(targetPosition)
        lowerBoundX = lowerBoundX - 450;
        upperBoundX = upperBoundX - 450;
    else
        targetPositionSegmentX = targetPosition(lowerBoundX:upperBoundX,1);
        targetPositionSegmentY = targetPosition(lowerBoundX:upperBoundX,2);
    end
    
    subplot(3,1,1)
    plot(lowerBoundX:upperBoundX,targetPositionSegmentX,'r')
    hold on
    %     if length(frameSelection) == 1
    %         line([frameSelection(1) frameSelection(1)],[min(targetPosition(:,1)) max(targetPosition(:,1))],'Color','g','LineStyle','--');
    %     end
    for n = 1:2:length(frameSelection)
        if lowerBoundX <= frameSelection(n) &&  frameSelection(n) <= upperBoundX
            line([frameSelection(n) frameSelection(n)],[min(targetPosition(:,1)) max(targetPosition(:,1))],'Color','g','LineStyle','--');
        end
    end
    for n = 2:2:length(frameSelection)
        if lowerBoundX <= frameSelection(n) &&  frameSelection(n) <= upperBoundX
            line([frameSelection(n) frameSelection(n)],[min(targetPosition(:,1)) max(targetPosition(:,1))],'Color','r','LineStyle','--');
        end
    end 
    hold off
    title(['\fontsize{20pt}\bf{Object Position per Frame}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
    xlabel('Frame')
    ylabel('X Position (Pixels)')
    axis([lowerBoundX upperBoundX floor(min([targetPosition(:,1);targetPosition(:,2)])/10)*10 ceil(max([targetPosition(:,1);targetPosition(:,2)])/10)*10])
    grid on
    
    subplot(3,1,2)
    plot(lowerBoundX:upperBoundX,targetPositionSegmentY,'b')
    hold on
    %     if length(frameSelection) == 1
    %         line([frameSelection(1) frameSelection(1)],[min(targetPosition(:,2)) max(targetPosition(:,2))],'Color','g','LineStyle','--');
    %     end
    for n = 1:2:length(frameSelection)
        if lowerBoundX <= frameSelection(n) &&  frameSelection(n) <= upperBoundX
            line([frameSelection(n) frameSelection(n)],[min(targetPosition(:,2)) max(targetPosition(:,2))],'Color','g','LineStyle','--');
        end
    end
    for n = 2:2:length(frameSelection)
        if lowerBoundX <= frameSelection(n) &&  frameSelection(n) <= upperBoundX
            line([frameSelection(n) frameSelection(n)],[min(targetPosition(:,2)) max(targetPosition(:,2))],'Color','r','LineStyle','--');
        end
    end
    hold off
    xlabel('Frame')
    ylabel('Y Position (Pixels)')
    axis([lowerBoundX upperBoundX floor(min([targetPosition(:,1);targetPosition(:,2)])/10)*10 ceil(max([targetPosition(:,1);targetPosition(:,2)])/10)*10])
    grid on
    
    secondsBoundsSelect = [(tifFrameBounds(1)+lowerBoundX-1)*secondsPerFrame (tifFrameBounds(1)+upperBoundX)*secondsPerFrame];
    ballDataSegmentIdx = secondsBoundsSelect(1) <= ballData(:,1) & ballData(:,1) <= secondsBoundsSelect(2);
    ballDataSegment = [ballData(ballDataSegmentIdx,1) ballData(ballDataSegmentIdx,2)];
    
    subplot(3,1,3)
    plot(ballDataSegment(:,1),ballDataSegment(:,2),'k')
    title('\fontsize{20pt}\bf{Ball Movement}')
    xlabel('Time (s)')
    ylabel('Movement')
    grid on
    axis([secondsBoundsSelect(1) secondsBoundsSelect(2) -1 ceil(max(ballData(:,2)))])
    
    [selectedCoordX,~,arrowKey] = ginput(1);
    
    if isempty(arrowKey)
        cont = false;
    elseif arrowKey == 1
        if ~isempty(frameSelection) && round(selectedCoordX) <= frameSelection(end)
            frameSelection = frameSelection(1:end-1);
        else
            frameSelection(end+1) = round(selectedCoordX);
        end
    elseif arrowKey == 28
        lowerBoundX = lowerBoundX - 450;
        upperBoundX = upperBoundX - 450;
    elseif arrowKey == 29
        lowerBoundX = lowerBoundX + 450;
        upperBoundX = upperBoundX + 450;
    end
end

avgPixelPos = [];
pixelDiffX = [];
pixelDiffY = [];
pixelDiffPosX = [];
pixelDiffPosY = [];
if length(frameSelection) < 4
    disp('Could not generate calibration values because too few points were selected')
    return
else
    if length(frameSelection)/2 ~= floor(length(frameSelection)/2)
        frameSelection = frameSelection(1:end-1);
    end
    for n = 1:length(frameSelection)/2
        avgPixelPos(n,:) = [mean(targetPosition(frameSelection(n*2-1):frameSelection(n*2),1)), mean(targetPosition(frameSelection(n*2-1):frameSelection(n*2),2))];
        targetAreaFrames = targetArea(frameSelection(n*2-1):frameSelection(n*2));
        targetAreaFrames(targetAreaFrames == -1) = [];
        if isempty(targetAreaFrames)
            avgTargetArea(n,1) = -1;
        else
            avgTargetArea(n,1) = mean(targetAreaFrames);
        end
        targetFwhmFrames = fwhmTargetPixels(:,:,frameSelection(n*2-1):frameSelection(n*2));
        avgTargetFwhm(n,1) = fwhm(1:size(targetFwhmFrames,2),mean(mean(targetFwhmFrames),3));
    end
    for n = 1:size(avgPixelPos,1)-1
        if abs(avgPixelPos(n+1,1)-avgPixelPos(n,1)) > abs(avgPixelPos(n+1,2)-avgPixelPos(n,2))
            pixelDiffX(end+1,1) = sqrt(((avgPixelPos(n+1,1)-avgPixelPos(n,1))^2) + ((avgPixelPos(n+1,2)-avgPixelPos(n,2))^2));
            pixelDiffPosX(end+1,:) = [mean([avgPixelPos(n+1,1), avgPixelPos(n,1)]), mean([avgPixelPos(n+1,2), avgPixelPos(n,2)])];
        else
            pixelDiffY(end+1,1) = sqrt(((avgPixelPos(n+1,1)-avgPixelPos(n,1))^2) + ((avgPixelPos(n+1,2)-avgPixelPos(n,2))^2));
            pixelDiffPosY(end+1,:) = [mean([avgPixelPos(n+1,1), avgPixelPos(n,1)]), mean([avgPixelPos(n+1,2), avgPixelPos(n,2)])];
        end
    end
end

h(4) = figure('Color','White');
k = convhull(targetPosition(:,1),targetPosition(:,2));
plot(targetPosition(k,1),targetPosition(k,2),'b',targetPosition(:,1),targetPosition(:,2),'k');
hold on
scatter([pixelDiffPosX(:,1); pixelDiffPosY(:,1)],[pixelDiffPosX(:,2); pixelDiffPosY(:,2)],50,[pixelDiffX; pixelDiffY],'filled');
colorbar;
hold off
axis equal
axis([1 W 1 H])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title(['\fontsize{20pt}\bf{Position of Target Object}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixels')
ylabel('Pixels')

h(5) = figure('Color','White');
hist3(targetPosition,'CdataMode','auto','Nbins',[40 40]);
colorbar;
view(2);
axis equal
axis([1 W 1 H])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title(['\fontsize{20pt}\bf{Position of Target Object Histogram}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixels (X)')
ylabel('Pixels (Y)')

h(6) = figure('Color','White');
hist(targetPosition(:,1),500);
title(['\fontsize{20pt}\bf{Position of Target Object Histogram (X)}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixels (X)')
ylabel('Number of Data Points')

h(7) = figure('Color','White');
hist(targetPosition(:,2),500);
title(['\fontsize{20pt}\bf{Position of Target Object Histogram (Y)}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixels (Y)')
ylabel('Number of Data Points')

h(8) = figure('Color','White');
plot(pixelDiffPosX(:,1),pixelDiffX,'b*');
axis([1 W 0 ceil(max(pixelDiffX))])
title(['\fontsize{20pt}\bf{Pixel Change in X for Specific \mum Jumps vs X Location in Image}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixel Position (X)')
ylabel('Pixel Change from Jump')

h(9) = figure('Color','White');
plot(pixelDiffPosY(:,2),pixelDiffY,'b*');
axis([1 H 0 ceil(max(pixelDiffY))])
title(['\fontsize{20pt}\bf{Pixel Change in Y for Specific \mum Jumps vs Y Location in Image}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixel Position (Y)')
ylabel('Pixel Change from Jump')

h(10) = figure('Color','White');
surfaceCalibFitX = fit(pixelDiffPosX,micronJumpVal./pixelDiffX,'poly55');
plot(surfaceCalibFitX,pixelDiffPosX,micronJumpVal./pixelDiffX);
axis([1 W 1 H 0 ceil(max(micronJumpVal./pixelDiffX))])
title(['\fontsize{20pt}\bf{Surface Calibration Fit (X)}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixels (X)')
ylabel('Pixels (Y)')
zlabel('\mum/Pixel')

h(11) = figure('Color','White');
surfaceCalibFitY = fit(pixelDiffPosY,micronJumpVal./pixelDiffY,'poly55');
plot(surfaceCalibFitY,pixelDiffPosY,micronJumpVal./pixelDiffY);
axis([1 W 1 H 0 ceil(max(micronJumpVal./pixelDiffY))])
title(['\fontsize{20pt}\bf{Surface Calibration Fit (Y)}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixels (X)')
ylabel('Pixels (Y)')
zlabel('\mum/Pixel')

h(12) = figure('Color','White');
subplot(2,1,1)
plot(1:length(targetArea),targetArea,'k');
axis([1 length(targetArea) -1 ceil(max(targetArea))])
title(['\fontsize{20pt}\bf{Tracked Target Area Calculations}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Frame')
ylabel('Area (Pixels)')
subplot(2,1,2)
plot(1:length(targetArea),medfilt1(lowpass(targetArea,1,framesPerSecond),5),'k');
axis([1 length(targetArea) -1 ceil(max(targetArea))])
title('LP and Median Filtered')
xlabel('Frame')
ylabel('Area (Pixels)')

h(13) = figure('Color','White');
plot(avgPixelPos(:,1),avgTargetArea,'b*');
axis([1 W 0 ceil(max(avgTargetArea))])
title(['\fontsize{20pt}\bf{Tracked Target Area vs X Location in Image}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixel Position (X)')
ylabel('Area (Pixels)')

h(14) = figure('Color','White');
plot(avgPixelPos(:,2),avgTargetArea,'b*');
axis([1 H 0 ceil(max(avgTargetArea))])
title(['\fontsize{20pt}\bf{Tracked Target Area vs Y Location in Image}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixel Position (Y)')
ylabel('Area (Pixels)')

h(15) = figure('Color','White');
scatter3(avgPixelPos(:,1),avgPixelPos(:,2),avgTargetArea,'b*');
axis([1 W 1 H 0 ceil(max(avgTargetArea))])
title(['\fontsize{20pt}\bf{Plot of Average Target Areas Throughout Image}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixels (X)')
ylabel('Pixels (Y)')
zlabel('Area (Pixels)')

h(16) = figure('Color','White');
plot(avgPixelPos(:,1),avgTargetFwhm,'b*');
axis([1 W 0 ceil(max(avgTargetFwhm))])
title(['\fontsize{20pt}\bf{Tracked Target FWHM vs X Location in Image}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixel Position (X)')
ylabel('Width (Pixels)')

h(17) = figure('Color','White');
plot(avgPixelPos(:,2),avgTargetFwhm,'b*');
axis([1 H 0 ceil(max(avgTargetFwhm))])
title(['\fontsize{20pt}\bf{Tracked Target FWHM vs Y Location in Image}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixel Position (Y)')
ylabel('Width (Pixels)')

h(18) = figure('Color','White');
scatter3(avgPixelPos(:,1),avgPixelPos(:,2),avgTargetFwhm,'b*');
axis([1 W 1 H 0 ceil(max(avgTargetFwhm))])
title(['\fontsize{20pt}\bf{Plot of Average Target FWHM Throughout Image}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixels (X)')
ylabel('Pixels (Y)')
zlabel('Width (Pixels)')

% h(5) = figure('Color','White');
% [uxy,~,idx] = unique([targetPosition(:,1),targetPosition(:,2)],'rows');
% szscale = histc(idx,unique(idx));
% [x1,y1] = meshgrid(-maxVal:.05:maxVal,-maxVal:.05:maxVal);
% z1 = griddata(uxy(:,1),uxy(:,2),szscale,x1,y1);
% pcolor(x1,y1,log(sqrt(z1)))
% shading interp
% axis equal square
% axis([-maxVal maxVal -maxVal maxVal])
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% title(['\fontsize{20pt}\bf{Position of Target Object}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
% xlabel('Pixels')
% ylabel('Pixels')

% Save figures to single .fig file
savefig(h,[tifFileName(1:end-4) '_calibrationPlots.fig']);

% Save calibration values
if exist('C:\Workspace\Code\DrewLab\calibrationValues.mat')
    load('C:\Workspace\Code\DrewLab\calibrationValues.mat');
end

if length(frameSelection) >= 4
    calibrationValues.(['file_' fileName]).pixelDiffX = pixelDiffX;
    calibrationValues.(['file_' fileName]).pixelDiffY = pixelDiffY;
    calibrationValues.(['file_' fileName]).avgPixelDiffX = mean(pixelDiffX);
    calibrationValues.(['file_' fileName]).avgPixelDiffY = mean(pixelDiffY);
    calibrationValues.(['file_' fileName]).micronJumpVal = micronJumpVal;
    calibrationValues.(['file_' fileName]).frames = tifFrameBounds;
    calibrationValues.(['file_' fileName]).pixelDiffPosX = pixelDiffPosX;
    calibrationValues.(['file_' fileName]).pixelDiffPosY = pixelDiffPosY;
    calibrationValues.(['file_' fileName]).avgTargetArea = avgTargetArea;
    calibrationValues.(['file_' fileName]).avgTargetFwhm = avgTargetFwhm;
    calibrationValues.(['file_' fileName]).surfaceCalibFitX = surfaceCalibFitX;
    calibrationValues.(['file_' fileName]).surfaceCalibFitY = surfaceCalibFitY;
    
    save('C:\Workspace\Code\DrewLab\calibrationValues.mat','calibrationValues');
end
end