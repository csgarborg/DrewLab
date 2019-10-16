%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: 
%________________________________________________________________________________________________________________________
%
%   Inputs: 
%
%   Outputs:        
%
%   Last Revised: June 30th, 2019    
%________________________________________________________________________________________________________________________

clear
clc

% User inputs for file information
pupilCamFileID = uigetfile('*_PupilCam.bin', 'MultiSelect', 'off');
animalID = input('Input the animal ID: ', 's'); disp(' ')
rawDataFileID = [animalID '_' pupilCamFileID(1:end - 13) '_RawData.mat'];

disp(['Loading relevant file information from ' rawDataFileID '...']); disp(' ')
try
    load(rawDataFileID)
catch
    disp([rawDataFileID ' does not appear to be in the current file path']); disp(' ')
    return
end

trialDuration = RawData.notes.trialDuration_sec;
disp([pupilCamFileID ' is ' num2str(trialDuration) ' seconds long.']); disp(' ')
startTime = input('Input the desired start time (sec): '); disp(' ')
endTime = input('Input the desired end time (sec): '); disp(' ')

if startTime >= trialDuration || startTime <= 0
    disp(['A start time of  ' num2str(startTime) ' is not a valid input']); disp(' ')
    return
elseif endTime > trialDuration || endTime <= startTime || endTime <= 0
    disp(['An end time of  ' num2str(startTime) ' is not a valid input']); disp(' ')
    return
end

imageHeight = RawData.notes.pupilCamPixelHeight;                                                                                                            
imageWidth = RawData.notes.pupilCamPixelWidth;
Fs = RawData.notes.pupilCamSamplingRate;

frameStart = floor(startTime)*Fs;
frameEnd = floor(endTime)*Fs;         
frameInds = frameStart:frameEnd;

pixelsPerFrame = imageWidth*imageHeight;
skippedPixels = pixelsPerFrame*2; % Multiply by two because there are 16 bits (2 bytes) per pixel
fid = fopen(pupilCamFileID);
fseek(fid, 0, 'eof');
fileSize = ftell(fid);
fseek(fid, 0, 'bof');
nFramesToRead = length(frameInds);
imageStack = zeros(imageWidth, imageHeight, nFramesToRead);
for a = 1:nFramesToRead
    disp(['Creating image stack: (' num2str(a) '/' num2str(nFramesToRead) ')']); disp(' ')
    fseek(fid, frameInds(a)*skippedPixels, 'bof');
    z = fread(fid, pixelsPerFrame, '*uint8', 'b');
    img = reshape(z(1:pixelsPerFrame), imageHeight, imageWidth);
    imageStack(:,:,a) = flip(imrotate(img, -90), 2);
end
fclose('all');

handle = implay(imageStack, Fs);
handle.Visual.ColorMap.UserRange = 1; 
handle.Visual.ColorMap.UserRangeMin = 0; 
handle.Visual.ColorMap.UserRangeMax = 255;