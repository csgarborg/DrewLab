%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:    calibrateFTLOpticalPowers
%
% FUNCTION:         calibrateFTLOpticalPowers(matFileName)
%
% DESCRIPTION:      Calibrates electrically focus-tunable lenses to relate
%                   diopter settings with actual measurements in microns in
%                   z for the two photon microscope (assumes diopter value
%                   = 0 is at the same level as z = 0 with positive z
%                   coming up towards the objective)
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
% WRITTEN BY:       Spencer Garborg 11/07/19
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calibrateFTLOpticalPowers(measurementStackFilename,diopterStackFilename,diopterVals,zMicronVals)

% import measurement stack
tifLength = length(imfinfo(measurementStackFilename));
% zMicronVals = [0:tifLength-1] * zStackIntMicrons;
measurementStack = [];
for i = 2:tifLength
    measurementStack(:,:,end+1) = im2double(imread(measurementStackFilename, i));
end

% import diopter stack
tifLength = length(imfinfo(diopterStackFilename));
diopterStack = [];
for i = 1:tifLength
    diopterStack(:,:,i) = im2double(imread(diopterStackFilename, i));
end

% calculate best fit for each image in diopter stack between all images in
% measurement stack
f = waitbar(0,'Calculating best match');
for n = 1:size(diopterStack,3)
    waitbar(round((n)/(size(diopterStack,3)),2),f,'Calculating best match');
    maxCVals = [];
    for i = 1:size(measurementStack,3)
        maxCVals(i) = max(max(normxcorr2(diopterStack(:,:,n),measurementStack(:,:,i))));
    end
    [~,maxI] = max(maxCVals);
    zMatchMicrons(n) = zMicronVals(maxI);
end
close(f)

% Get file name
[tokens,~] = regexpi(measurementStackFilename,'\\([^\\]*).TIF','tokens','match');
fileName = tokens{1}{1};

% Save calibration values
if exist('C:\Workspace\Code\DrewLab\calibrationValuesFTL.mat','file')
    load('C:\Workspace\Code\DrewLab\calibrationValuesFTL.mat');
end

calibrationValues.(['file_' fileName]).diopterVals = diopterVals;
calibrationValues.(['file_' fileName]).zMatchMicrons = zMatchMicrons;

save('C:\Workspace\Code\DrewLab\calibrationValuesFTL.mat','calibrationValues');

plotFTLCalibration
end