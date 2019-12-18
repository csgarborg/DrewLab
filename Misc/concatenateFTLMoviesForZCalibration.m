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
function concatenateFTLMoviesForZCalibration(outputFilename,inputFilenameCellArray)

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

