%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:    lowpassImageFilter2P
%
% FUNCTION:         lowpassImageFilter2P(tifFileName)
%
% DESCRIPTION:      Creates an AVI file of lowpass filtered images from a
%                   TIF file from a 2P microscope to filter out pixel
%                   flickering and allow for better motion tracking
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
% WRITTEN BY:       Spencer Garborg 1/9/19
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function aviFileName = lowpassImageFilter2P(tifFileName,redoFilterTF,framesVec)
%% Check to see if overwriting AVI file or using 

aviFileName = [tifFileName(1:end-4) '_filtered.avi'];

if exist(aviFileName,'file') && ~redoFilterTF
    return
end

%% Convert TIF image stack to matrix

if ~exist('framesVec','var')
    tifLength = length(imfinfo(tifFileName));
    framesVec = [1 tifLength];
% tifLength = 100;
else
    tifLength = framesVec(2) - framesVec(1) + 1;
end

imageStack = cell(1,tifLength);
for k = framesVec(1):framesVec(2)
    imageStack{1,k-framesVec(1)+1} = imread(tifFileName, k);
end

imageStackMat = cell2mat(imageStack);
imageStackMat = double(reshape(imageStackMat,512,512,tifLength));

%% Create filter

lpFilt = designfilt('lowpassiir','FilterOrder',8, ...
         'PassbandFrequency',5,'PassbandRipple',0.2, ...
         'SampleRate',30);
     
%% Filter images

filteredData = zeros(512,512,tifLength);
f = waitbar(0,'Filtering');
for xx = 1:512
    waitbar(round(xx/512,2),f,'Filtering');
    for yy = 1:512
        singlePixelData = squeeze(imageStackMat(xx,yy,:));
        filteredData(xx,yy,:) = filter(lpFilt,singlePixelData);
    end
end
close(f)

%% Create AVI file frame by frame

if exist(aviFileName,'file')
    delete(aviFileName)
end

aviObject = VideoWriter(aviFileName,'Uncompressed AVI');
open(aviObject);
 
f = waitbar(0,'Creating AVI file');
for k = 1:tifLength
    waitbar(round(k/tifLength,2),f,'Creating AVI file');
%     imagesc(filteredData(:,:,k));
%     filteredFrame = getframe(fig,[0.05 0.05 0.9 0.9]);
    writeVideo(aviObject,mat2gray(filteredData(:,:,k)));
end
close(f)
close(aviObject);
end