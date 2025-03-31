%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:    concatenateFTLMoviesForZCalibration
%
% FUNCTION:         concatenateFTLMoviesForZCalibration(outputFilename,inputFilenameCellArray)
%
% DESCRIPTION:      Creates single file from multiple movies at different
%                   diopter values by averaging each movie and
%                   concatenating the averaged frames in the order of the
%                   cell array input
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
function concatenateFTLMoviesForZCalibration(outputFilename,filenameStem,filenameVec)
for n = 1:length(filenameVec)
    if filenameVec(n) < 10
        inputFilenameCellArray{n} = [filenameStem '00' num2str(filenameVec(n)) '.tif'];
    elseif filenameVec(n) < 100
        inputFilenameCellArray{n} = [filenameStem '0' num2str(filenameVec(n)) '.tif'];
    else
        inputFilenameCellArray{n} = [filenameStem num2str(filenameVec(n)) '.tif'];
    end
end

for n = 1:numel(inputFilenameCellArray)
    tifLength = length(imfinfo(inputFilenameCellArray{n}));
    imageSum = zeros(512,512);
    for i = 2:tifLength
        imageSum = imageSum + im2double(imread(inputFilenameCellArray{n}, i));
    end
    imageAvg = imageSum ./ (tifLength-1);
    if n == 1
        imwrite(imageAvg,outputFilename);
    else
        imwrite(imageAvg,outputFilename,'WriteMode','append');
    end
end
end

