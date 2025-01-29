%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:    motionTracking2P3D
%
% FUNCTION:         motionTracking2P3D(tifFileName,calibrationFileString,redoFilterTF,updateSearchTF,medFiltTF,saveOutputVideoTF,threeStepTF,compiledTifTF,targetAvgNum,framesPerSecond,objMag,digMag,turnabout,commentString,tifFrameBounds)
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
% WRITTEN BY:       Spencer Garborg 5/20/20
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function motionTracking2P3D(tifFileName,calibrationFileString,updateSearchTF,medFiltTF,saveOutputVideoTF,threeStepTF,compiledTifTF,targetAvgNum,framesPerSecond,objMag,digMag,turnabout,layers,micronsPerLayer,hemisphere,commentString,tifFrameBounds)
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
    hTMXY = vision.TemplateMatcher('ROIInputPort', true, 'BestMatchNeighborhoodOutputPort', true, 'SearchMethod', 'Three-step');
else
    hTMXY = vision.TemplateMatcher('ROIInputPort', true, 'BestMatchNeighborhoodOutputPort', true);
end

if threeStepTF
    hTMXZ = vision.TemplateMatcher('ROIInputPort', true, 'BestMatchNeighborhoodOutputPort', true, 'SearchMethod', 'Three-step');
else
    hTMXZ = vision.TemplateMatcher('ROIInputPort', true, 'BestMatchNeighborhoodOutputPort', true);
end

if threeStepTF
    hTMYZ = vision.TemplateMatcher('ROIInputPort', true, 'BestMatchNeighborhoodOutputPort', true, 'SearchMethod', 'Three-step');
else
    hTMYZ = vision.TemplateMatcher('ROIInputPort', true, 'BestMatchNeighborhoodOutputPort', true);
end
                        
%% Create motion tracking output display
% Create a System object to display the original video and the stabilized
% video.
% hVideoOutXY = vision.VideoPlayer('Name', 'Video Stabilization');
% hVideoOutXY.Position(1) = round(0.01*hVideoOutXY.Position(1));
% hVideoOutXY.Position(2) = round(1.1*(hVideoOutXY.Position(2)));
% hVideoOutXY.Position(3:4) = [950 550];
% 
% hVideoOutXZ = vision.VideoPlayer('Name', 'Video Stabilization');
% hVideoOutXZ.Position(1) = round(1.25*hVideoOutXY.Position(1));
% hVideoOutXZ.Position(2) = round(1.1*(hVideoOutXY.Position(2)));
% hVideoOutXZ.Position(3:4) = [950 550];
% 
% hVideoOutYZ = vision.VideoPlayer('Name', 'Video Stabilization');
% hVideoOutYZ.Position(1) = round(.4*hVideoOutXY.Position(1));
% hVideoOutYZ.Position(2) = round(.1*(hVideoOutXY.Position(2)));
% hVideoOutYZ.Position(3:4) = [950 550];

hVideoOutXY = vision.VideoPlayer('Name', 'Video Stabilization XY Plane');
hVideoOutXY.Position(1) = round(.5*hVideoOutXY.Position(1));
hVideoOutXY.Position(2) = round(1.1*(hVideoOutXY.Position(2)));
hVideoOutXY.Position(3:4) = [1050 550];

hVideoOutXZ = vision.VideoPlayer('Name', 'Video Stabilization XZ Plane');
hVideoOutXZ.Position(1) = round(.5*hVideoOutXZ.Position(1));
hVideoOutXZ.Position(2) = round(.6*(hVideoOutXZ.Position(2)));
hVideoOutXZ.Position(3:4) = [1050 100];

hVideoOutYZ = vision.VideoPlayer('Name', 'Video Stabilization YZ Plane');
hVideoOutYZ.Position(1) = round(.5*hVideoOutYZ.Position(1));
hVideoOutYZ.Position(2) = round(.1*(hVideoOutYZ.Position(2)));
hVideoOutYZ.Position(3:4) = [1050 100];

%% Initialize variables for processing loop
initialImage = [];
if medFiltTF
    for n = 1:layers
        initialImage(:,:,n) = im2double(medfilt2(imread(tifFileName, tifFrameBounds(1)-1+n)));
    end
else
    for n = 1:layers
        initialImage(:,:,n) = im2double(imread(tifFileName, tifFrameBounds(1)-1+n));
    end
end
%%%%%%%%%%%% XY %%%%%%%%%%%%%%
% Get target window
initialImageXY = MIPXY(initialImage);
imshow(initialImageXY);
title('Select upper left, then lower right target corners and press enter');
[inputCoordTargetX,inputCoordTargetY] = getpts(gcf);
close(gcf);
inputCoordTargetX = round(inputCoordTargetX);
inputCoordTargetY = round(inputCoordTargetY);
posXY.template_orig = [inputCoordTargetX(1) inputCoordTargetY(1)]; % [x y] upper left corner
posXY.template_size = [inputCoordTargetX(2) inputCoordTargetY(2)] - [inputCoordTargetX(1) inputCoordTargetY(1)];   % [width height]

% Get search window
imshow(initialImageXY);
title('Select upper left search area corner (target box pictured) and press enter');
rectangle('Position',[posXY.template_orig(1) posXY.template_orig(2) posXY.template_size(1) posXY.template_size(2)],'EdgeColor','w');
[inputCoordSearchX,inputCoordSearchY] = getpts(gcf);
close(gcf);
inputCoordSearchX = round(inputCoordSearchX);
inputCoordSearchY = round(inputCoordSearchY);
posXY.search_border = [abs(inputCoordSearchX(1) - inputCoordTargetX(1)),abs(inputCoordSearchY(1) - inputCoordTargetY(1))];   % max horizontal and vertical displacement

% Calculate important parameters
posXY.template_center = floor((posXY.template_size-1)/2);
posXY.template_center_pos = (posXY.template_orig + posXY.template_center - 1);
W = size(initialImageXY,2); % Width of video in pixels
H = size(initialImageXY,1); % Height of video in pixels
sz = [W,H];
BorderCols = [1:posXY.search_border(1)+4 W-posXY.search_border(1)+4:W];
BorderRows = [1:posXY.search_border(2)+4 H-posXY.search_border(2)+4:H];
TargetRowIndices = ...
  posXY.template_orig(2)-1:posXY.template_orig(2)+posXY.template_size(2)-2;
TargetColIndices = ...
  posXY.template_orig(1)-1:posXY.template_orig(1)+posXY.template_size(1)-2;
SearchRegion = posXY.template_orig - posXY.search_border - 1;
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

%%%%%%%%%%%% XZ %%%%%%%%%%%%%%
% Get target window
initialImageXZ = MIPXZ(initialImage);
imshow(initialImageXZ);
title('Select upper left, then lower right target corners and press enter');
[inputCoordTargetX,inputCoordTargetY] = getpts(gcf);
close(gcf);
inputCoordTargetX = round(inputCoordTargetX);
inputCoordTargetY = round(inputCoordTargetY);
posXZ.template_orig = [inputCoordTargetX(1) inputCoordTargetY(1)]; % [x y] upper left corner
posXZ.template_size = [inputCoordTargetX(2) inputCoordTargetY(2)] - [inputCoordTargetX(1) inputCoordTargetY(1)];   % [width height]

% Get search window
imshow(initialImageXZ);
title('Select upper left search area corner (target box pictured) and press enter');
rectangle('Position',[posXZ.template_orig(1) posXZ.template_orig(2) posXZ.template_size(1) posXZ.template_size(2)],'EdgeColor','w');
[inputCoordSearchX,inputCoordSearchY] = getpts(gcf);
close(gcf);
inputCoordSearchX = round(inputCoordSearchX);
inputCoordSearchY = round(inputCoordSearchY);
posXZ.search_border = [abs(inputCoordSearchX(1) - inputCoordTargetX(1)),abs(inputCoordSearchY(1) - inputCoordTargetY(1))];   % max horizontal and vertical displacement

% Calculate important parameters
posXZ.template_center = floor((posXZ.template_size-1)/2);
posXZ.template_center_pos = (posXZ.template_orig + posXZ.template_center - 1);
W = size(initialImageXZ,2); % Width of video in pixels
HXZ = size(initialImageXZ,1); % Height of video in pixels
szXZ = [W,HXZ];
BorderColsXZ = [1:posXZ.search_border(1)+4 W-posXZ.search_border(1)+4:W];
BorderRowsXZ = [1:posXZ.search_border(2)+4 HXZ-posXZ.search_border(2)+4:HXZ];
TargetRowIndicesXZ = ...
  posXZ.template_orig(2)-1:posXZ.template_orig(2)+posXZ.template_size(2)-2;
TargetColIndicesXZ = ...
  posXZ.template_orig(1)-1:posXZ.template_orig(1)+posXZ.template_size(1)-2;
SearchRegionXZ = posXZ.template_orig - posXZ.search_border - 1;
OffsetXZ = [0 0];
TargetXZ = zeros(length(TargetRowIndicesXZ),length(TargetColIndicesXZ));
firstTimeXZ = true;
targetNumXZ = 0;
targetSumXZ = zeros(length(TargetRowIndicesXZ),length(TargetColIndicesXZ));
moveDistXZ = [];
velocityXZ = [];
% surfaceCalibFitX = calibrationValues.(calibrationFileString).surfaceCalibFitX;
% surfaceCalibFitY = calibrationValues.(calibrationFileString).surfaceCalibFitY;
midlineXXZ = round(W/2);
midlineYXZ = round(HXZ/2);

%%%%%%%%%%%% YZ %%%%%%%%%%%%%%
% Get target window
initialImageYZ = MIPYZ(initialImage);
imshow(initialImageYZ);
title('Select upper left, then lower right target corners and press enter');
[inputCoordTargetX,inputCoordTargetY] = getpts(gcf);
close(gcf);
inputCoordTargetX = round(inputCoordTargetX);
inputCoordTargetY = round(inputCoordTargetY);
posYZ.template_orig = [inputCoordTargetX(1) inputCoordTargetY(1)]; % [x y] upper left corner
posYZ.template_size = [inputCoordTargetX(2) inputCoordTargetY(2)] - [inputCoordTargetX(1) inputCoordTargetY(1)];   % [width height]

% Get search window
imshow(initialImageYZ);
title('Select upper left search area corner (target box pictured) and press enter');
rectangle('Position',[posYZ.template_orig(1) posYZ.template_orig(2) posYZ.template_size(1) posYZ.template_size(2)],'EdgeColor','w');
[inputCoordSearchX,inputCoordSearchY] = getpts(gcf);
close(gcf);
inputCoordSearchX = round(inputCoordSearchX);
inputCoordSearchY = round(inputCoordSearchY);
posYZ.search_border = [abs(inputCoordSearchX(1) - inputCoordTargetX(1)),abs(inputCoordSearchY(1) - inputCoordTargetY(1))];   % max horizontal and vertical displacement

% Calculate important parameters
posYZ.template_center = floor((posYZ.template_size-1)/2);
posYZ.template_center_pos = (posYZ.template_orig + posYZ.template_center - 1);
W = size(initialImageYZ,2); % Width of video in pixels
HYZ = size(initialImageYZ,1); % Height of video in pixels
szYZ = [W,HYZ];
BorderColsYZ = [1:posYZ.search_border(1)+4 W-posYZ.search_border(1)+4:W];
BorderRowsYZ = [1:posYZ.search_border(2)+4 HYZ-posYZ.search_border(2)+4:HYZ];
TargetRowIndicesYZ = ...
  posYZ.template_orig(2)-1:posYZ.template_orig(2)+posYZ.template_size(2)-2;
TargetColIndicesYZ = ...
  posYZ.template_orig(1)-1:posYZ.template_orig(1)+posYZ.template_size(1)-2;
SearchRegionYZ = posYZ.template_orig - posYZ.search_border - 1;
OffsetYZ = [0 0];
TargetYZ = zeros(length(TargetRowIndicesYZ),length(TargetColIndicesYZ));
firstTimeYZ = true;
targetNumYZ = 0;
targetSumYZ = zeros(length(TargetRowIndicesYZ),length(TargetColIndicesYZ));
moveDistYZ = [];
velocityYZ = [];
% surfaceCalibFitX = calibrationValues.(calibrationFileString).surfaceCalibFitX;
% surfaceCalibFitY = calibrationValues.(calibrationFileString).surfaceCalibFitY;
midlineXYZ = round(W/2);
midlineYYZ = round(HYZ/2);

%% Stream Processing Loop
% This is the main processing loop which uses the objects we instantiated
% above to stabilize the input video.
for i = tifFrameBounds(1):layers:tifFrameBounds(2)
    if medFiltTF
        for j = 1:layers
            input(:,:,j) = im2double(medfilt2(imread(tifFileName, i-1+j)));
        end
    else
        for j = 1:layers
            input(:,:,j) = im2double(imread(tifFileName, i-1+j));
        end
    end
    inputXY = MIPXY(input);
    inputXZ = MIPXZ(input);
    inputYZ = MIPYZ(input);
    %%%%%%%%%%%% XY %%%%%%%%%%%%%%
    % Find location of Target in the input video frame
    if firstTime
      Idx = int32(posXY.template_center_pos);
      motionVector = [0 0];
    else
      IdxPrev = Idx;
% IdxPrev = int32(pos.template_center_pos);

      ROI = [SearchRegion, posXY.template_size+2*posXY.search_border];
      Idx = hTMXY(inputXY,Target,ROI);
      
      motionVector = double(Idx-IdxPrev);
    end

    if updateSearchTF
        [Offset, SearchRegion] = updatesearch(sz, motionVector, ...
            SearchRegion, Offset, posXY);
    else
        [Offset] = updatesearch(sz, motionVector, ...
            SearchRegion, Offset, posXY);
    end

    % Translate video frame to offset the camera motion
    Stabilized = imtranslate(inputXY, Offset, 'linear');
    
    targetNum = targetNum + 1;
    if targetNum <= targetAvgNum
        targetSum = targetSum + Stabilized(round(TargetRowIndices), round(TargetColIndices));
        Target = targetSum ./ targetNum;
%         imshow(Target);
    end

    % Add black border for display
    Stabilized(:, BorderCols) = 0;
    Stabilized(BorderRows, :) = 0;
    txt = sprintf('Stabilized');
    Stabilized = insertText(Stabilized(:,:,1),[1 1],txt,'FontSize',16, ...
                    'TextColor', 'white', 'BoxOpacity', 0);

    TargetRect = [posXY.template_orig-Offset, posXY.template_size];
    SearchRegionRect = [SearchRegion, posXY.template_size + 2*posXY.search_border];

    % Draw rectangles on input to show target and search region
    inputXY = insertShape(inputXY, 'Rectangle', [TargetRect; SearchRegionRect],...
                        'Color', 'white');
    % Display the offset (displacement) values on the input image
    txt = sprintf('(%+05.1f pixels,%+05.1f pixels)', Offset);
    inputXY = insertText(inputXY(:,:,1),[1 1],txt,'FontSize',16, ...
                    'TextColor', 'white', 'BoxOpacity', 0);
    
    % Save output video to variable
    if saveOutputVideoTF
        imageStack{1,n} = [inputXY(:,:,1) Stabilized];
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
    %%%%%%%%%%%% XZ %%%%%%%%%%%%%%
    % Find location of Target in the input video frame
    if firstTimeXZ
      IdxXZ = int32(posXZ.template_center_pos);
      motionVectorXZ = [0 0];
    else
      IdxXZPrev = IdxXZ;
% IdxPrev = int32(pos.template_center_pos);

      ROI = [SearchRegionXZ, posXZ.template_size+2*posXZ.search_border];
      IdxXZ = hTMXY(inputXZ,TargetXZ,ROI);
      
      motionVectorXZ = double(IdxXZ-IdxXZPrev);
    end

    if updateSearchTF
        [OffsetXZ, SearchRegionXZ] = updatesearch(szXZ, motionVectorXZ, ...
            SearchRegionXZ, OffsetXZ, posXZ);
    else
        [OffsetXZ] = updatesearch(szXZ, motionVectorXZ, ...
            SearchRegionXZ, OffsetXZ, posXZ);
    end

    % Translate video frame to offset the camera motion
    StabilizedXZ = imtranslate(inputXZ, OffsetXZ, 'linear');
    
    targetNumXZ = targetNumXZ + 1;
    if targetNumXZ <= targetAvgNum
        targetSumXZ = targetSumXZ + StabilizedXZ(round(TargetRowIndicesXZ), round(TargetColIndicesXZ));
        TargetXZ = targetSumXZ ./ targetNumXZ;
%         imshow(TargetXZ);
    end

    % Add black border for display
    StabilizedXZ(:, BorderColsXZ) = 0;
    StabilizedXZ(BorderRowsXZ, :) = 0;
    txt = sprintf('Stabilized');
    StabilizedXZ = insertText(StabilizedXZ(:,:,1),[1 1],txt,'FontSize',10, ...
                    'TextColor', 'white', 'BoxOpacity', 0);

    TargetRect = [posXZ.template_orig-OffsetXZ, posXZ.template_size];
    SearchRegionRect = [SearchRegionXZ, posXZ.template_size + 2*posXZ.search_border];

    % Draw rectangles on input to show target and search region
    inputXZ = insertShape(inputXZ, 'Rectangle', [TargetRect; SearchRegionRect],...
                        'Color', 'white');
    % Display the offset (displacement) values on the input image
    txt = sprintf('(%+05.1f pixels,%+05.1f pixels)', OffsetXZ);
    inputXZ = insertText(inputXZ(:,:,1),[1 1],txt,'FontSize',10, ...
                    'TextColor', 'white', 'BoxOpacity', 0);
    
%     % Save output video to variable
%     if saveOutputVideoTF
%         imageStack{1,n} = [inputXZ(:,:,1) StabilizedXZ];
%         n = n + 1;
%     end
%     
    % Add pixel motion to data
    if firstTimeXZ
        targetPositionPixelXZ = [double(IdxXZ(1)), HXZ-double(IdxXZ(2))+1];
        targetPositionXZ = [0 0];
        firstTimeXZ = false;
    else
        motionVectorXZ(2) = -motionVectorXZ(2);
        targetPositionPixelXZ = [targetPositionPixelXZ;targetPositionPixelXZ(end,:)+motionVectorXZ];
    end
    %%%%%%%%%%%% YZ %%%%%%%%%%%%%%
    % Find location of Target in the input video frame
    if firstTimeYZ
      IdxYZ = int32(posYZ.template_center_pos);
      motionVectorYZ = [0 0];
    else
      IdxYZPrev = IdxYZ;
% IdxPrev = int32(pos.template_center_pos);

      ROI = [SearchRegionYZ, posYZ.template_size+2*posYZ.search_border];
      IdxYZ = hTMXY(inputYZ,TargetYZ,ROI);
      
      motionVectorYZ = double(IdxYZ-IdxYZPrev);
    end

    if updateSearchTF
        [OffsetYZ, SearchRegionYZ] = updatesearch(szYZ, motionVectorYZ, ...
            SearchRegionYZ, OffsetYZ, posYZ);
    else
        [OffsetYZ] = updatesearch(szYZ, motionVectorYZ, ...
            SearchRegionYZ, OffsetYZ, posYZ);
    end

    % Translate video frame to offset the camera motion
    StabilizedYZ = imtranslate(inputYZ, OffsetYZ, 'linear');
    
    targetNumYZ = targetNumYZ + 1;
    if targetNumYZ <= targetAvgNum
        targetSumYZ = targetSumYZ + StabilizedYZ(round(TargetRowIndicesYZ), round(TargetColIndicesYZ));
        TargetYZ = targetSumYZ ./ targetNumXZ;
%         imshow(TargetXZ);
    end

    % Add black border for display
    StabilizedYZ(:, BorderColsYZ) = 0;
    StabilizedYZ(BorderRowsYZ, :) = 0;
    txt = sprintf('Stabilized');
    StabilizedYZ = insertText(StabilizedYZ(:,:,1),[1 1],txt,'FontSize',10, ...
                    'TextColor', 'white', 'BoxOpacity', 0);

    TargetRect = [posYZ.template_orig-OffsetYZ, posYZ.template_size];
    SearchRegionRect = [SearchRegionYZ, posYZ.template_size + 2*posYZ.search_border];

    % Draw rectangles on input to show target and search region
    inputYZ = insertShape(inputYZ, 'Rectangle', [TargetRect; SearchRegionRect],...
                        'Color', 'white');
    % Display the offset (displacement) values on the input image
    txt = sprintf('(%+05.1f pixels,%+05.1f pixels)', OffsetYZ);
    inputYZ = insertText(inputYZ(:,:,1),[1 1],txt,'FontSize',10, ...
                    'TextColor', 'white', 'BoxOpacity', 0);
    
%     % Save output video to variable
%     if saveOutputVideoTF
%         imageStack{1,n} = [inputYZ(:,:,1) StabilizedYZ];
%         n = n + 1;
%     end
%     
    % Add pixel motion to data
    if firstTimeYZ
        targetPositionPixelYZ = [double(IdxYZ(1)), HYZ-double(IdxYZ(2))+1];
        targetPositionYZ = [0 0];
        firstTimeYZ = false;
    else
        motionVectorYZ(2) = -motionVectorYZ(2);
        targetPositionPixelYZ = [targetPositionPixelYZ;targetPositionPixelYZ(end,:)+motionVectorYZ];
    end
    
    %%%%%%%%%%% Display video %%%%%%%%%%%%%
    hVideoOutXY([inputXY(:,:,1) Stabilized(:,:,1)]);
    hVideoOutXZ([zeros(50,2*size(inputXZ,2)); inputXZ(:,:,1) StabilizedXZ(:,:,1); zeros(50,2*size(inputXZ,2))]);
    hVideoOutYZ([zeros(50,2*size(inputXZ,2)); inputYZ(:,:,1) StabilizedYZ(:,:,1); zeros(50,2*size(inputXZ,2))]);
    
%     [imind,cm] = rgb2ind(inputXY,256);
%     Stabilized = ind2rgb(Stabilized,cm);
%     [imind,cm] = rgb2ind([inputXY Stabilized],256);
%     if n == 1
%         imwrite(imind,cm,'Y:\Figures\PSFBeadMove3DMIPXY.gif','gif','DelayTime',.1, 'Loopcount',inf);
%     else
%         imwrite(imind,cm,'Y:\Figures\PSFBeadMove3DMIPXY.gif','gif','DelayTime',.1, 'WriteMode','append');
%     end
%     close
%     
% %     StabilizedXZ = ind2rgb(StabilizedXZ,cm);
%     [imind,cm] = rgb2ind([inputXZ StabilizedXZ],256);
%     imind = [zeros(50,2*size(inputXZ,2)); imind; zeros(50,2*size(inputXZ,2))];
%     if n == 1
%         imwrite(imind,cm,'Y:\Figures\PSFBeadMove3DMIPXZ.gif','gif','DelayTime',.1, 'Loopcount',inf);
%     else
%         imwrite(imind,cm,'Y:\Figures\PSFBeadMove3DMIPXZ.gif','gif','DelayTime',.1, 'WriteMode','append');
%     end
%     close
%     
% %     StabilizedYZ = ind2rgb(StabilizedYZ,cm);
%     [imind,cm] = rgb2ind([inputYZ StabilizedYZ],256);
%     imind = [zeros(50,2*size(inputXZ,2)); imind; zeros(50,2*size(inputXZ,2))];
%     if n == 1
%         imwrite(imind,cm,'Y:\Figures\PSFBeadMove3DMIPYZ.gif','gif','DelayTime',.1, 'Loopcount',inf);
%     else
%         imwrite(imind,cm,'Y:\Figures\PSFBeadMove3DMIPYZ.gif','gif','DelayTime',.1, 'WriteMode','append');
%     end
%     close
%     
%     n = n + 1;
end

% Perform movement calculations
movementLength = size(targetPositionPixel,1);
f = waitbar(0,'Calculating position from calibration');
[meanXCalibVec,meanYCalibVec] = getMeanCalibVecs(size(Stabilized,2),size(Stabilized,1),surfaceCalibFitX,surfaceCalibFitY);
for n = 2:movementLength
    %%%%%%%%%%%% XY %%%%%%%%%%%%%%
    % Get x movement
    micronDistTraveledX = getDistFromImCentX(targetPositionPixel(n-1,1),targetPositionPixel(n-1,2),targetPositionPixel(n,1),targetPositionPixel(n,2),midlineX,surfaceCalibFitX);
    
    % Get y movement
    micronDistTraveledY = getDistFromImCentY(targetPositionPixel(n-1,1),targetPositionPixel(n-1,2),targetPositionPixel(n,1),targetPositionPixel(n,2),midlineY,surfaceCalibFitY);
    
    moveDist(n-1,:) = [micronDistTraveledX,micronDistTraveledY];
    velocity(n-1,1) = sqrt((micronDistTraveledX)^2+(micronDistTraveledY)^2)/(layers*secondsPerFrame);
    targetPosition(n,:) = targetPosition(n-1,:)+[micronDistTraveledX,micronDistTraveledY];
    
    %%%%%%%%%%%% XZ %%%%%%%%%%%%%%
    % Get x movement
    micronDistTraveledX = getDistFromImCentX_XZ(targetPositionPixelXZ(n-1,1),targetPositionPixelXZ(n,1),midlineX,meanXCalibVec);
    % Get z movement
    micronDistTraveledZ = (targetPositionPixelXZ(n,2) - targetPositionPixelXZ(n-1,2))*micronsPerLayer;
    
    moveDistXZ(n-1,:) = [micronDistTraveledX,micronDistTraveledZ];
    velocityXZ(n-1,1) = sqrt((micronDistTraveledX)^2+(micronDistTraveledZ)^2)/(layers*secondsPerFrame);
    targetPositionXZ(n,:) = targetPositionXZ(n-1,:)+[micronDistTraveledX,micronDistTraveledZ];
    
    %%%%%%%%%%%% YZ %%%%%%%%%%%%%%
    % Get y movement
    micronDistTraveledY = -getDistFromImCentY_YZ(targetPositionPixelYZ(n-1,1),targetPositionPixelYZ(n,1),midlineY,meanYCalibVec);
    % Get z movement
    micronDistTraveledZ = (targetPositionPixelYZ(n,2) - targetPositionPixelYZ(n-1,2))*micronsPerLayer;
    
    moveDistYZ(n-1,:) = [micronDistTraveledY,micronDistTraveledZ];
    velocityYZ(n-1,1) = sqrt((micronDistTraveledY)^2+(micronDistTraveledZ)^2)/(layers*secondsPerFrame);
    targetPositionYZ(n,:) = targetPositionYZ(n-1,:)+[micronDistTraveledY,micronDistTraveledZ];
    
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
movementData.moveDistXZ = moveDistXZ;
movementData.velocityXZ = velocityXZ;
movementData.targetPositionXZ = targetPositionXZ;
movementData.moveDistYZ = moveDistYZ;
movementData.velocityYZ = velocityYZ;
movementData.targetPositionYZ = targetPositionYZ;
movementData.calibrationFileString = calibrationFileString;
movementData.frames = tifFrameBounds;
movementData.imageSize = sz;
movementData.medFiltTF = medFiltTF;
movementData.pos = posXY;
movementData.ballData = ballData;
movementData.secondsPerFrame = secondsPerFrame*layers;
movementData.objMag = objMag;
movementData.digMag = digMag;
movementData.turnabout = turnabout;
movementData.commentString = commentString;
movementData.layers = layers;
movementData.micronsPerLayer = micronsPerLayer;
movementData.hemisphere = hemisphere;

matFileName = [tifFileName(1:end-4) '3D.mat'];
save(matFileName,'movementData');

% Plot data
plotMotionTracking3D(matFileName);
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