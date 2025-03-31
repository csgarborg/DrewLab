%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:    imageFilter2P
%
% FUNCTION:         imageFilter2P(tifFileName,redoFilterTF,medFiltTF,compiledTifTF,framesVec)
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

function aviFileName = imageFilter2P(tifFileName,redoFilterTF,medFiltTF,compiledTifTF,framesVec)
%% Check to see if overwriting AVI file or using 

aviFileName = [tifFileName(1:end-4) '_filtered.avi'];

if exist(aviFileName,'file') && ~redoFilterTF
    return
end

%% Convert TIF image stack to matrix
if compiledTifTF
    [imStack,tifLength] = imread_big(tifFileName);
end

if ~exist('framesVec','var')
    if compiledTifTF
        tifLength = tifLength-1;
    else
        tifLength = length(imfinfo(tifFileName)) - 1;
    end
    framesVec = [2 tifLength+1];
else
    tifLength = framesVec(2) - framesVec(1) + 1;
end

filteredData = cell(1,tifLength);
if medFiltTF
    f = waitbar(0,'Compiling and Filtering Images');
    if compiledTifTF
        for k = framesVec(1):framesVec(2)
            filteredData{1,k-framesVec(1)+1} = im2uint8(medfilt2(imStack(:,:,k)));
            waitbar(round((k-framesVec(1))/tifLength,2),f,'Compiling and Filtering Images');
        end
    else
        for k = framesVec(1):framesVec(2)
            filteredData{1,k-framesVec(1)+1} = im2uint8(medfilt2(imread(tifFileName, k)));
            waitbar(round((k-framesVec(1))/tifLength,2),f,'Compiling and Filtering Images');
        end
    end
else
    f = waitbar(0,'Compiling Unfiltered Images');
    if compiledTifTF
        for k = framesVec(1):framesVec(2)
            filteredData{1,k-framesVec(1)+1} = im2uint8(imStack(:,:,k));
            waitbar(round((k-framesVec(1))/tifLength,2),f,'Compiling Unfiltered Images');
        end
    else
        for k = framesVec(1):framesVec(2)
            filteredData{1,k-framesVec(1)+1} = im2uint8(imread(tifFileName, k));
            waitbar(round((k-framesVec(1))/tifLength,2),f,'Compiling Unfiltered Images');
        end
    end
end
close(f)

% [H,W] = size(imageStack{1,1});
% imageStackMat = cell2mat(imageStack);
% imageStackMat = double(reshape(imageStackMat,H,W,tifLength));

%% Create filter

% lpFilt = designfilt('lowpassiir','FilterOrder',3, ...
%          'PassbandFrequency',2,'PassbandRipple',0.2, ...
%          'SampleRate',30);
%      
%% Filter images

% filteredData = zeros(512,512,tifLength);
% f = waitbar(0,'Filtering');
% for xx = 1:512
%     waitbar(round(xx/512,2),f,'Filtering');
%     for yy = 1:512
%         singlePixelData = squeeze(imageStackMat(xx,yy,:));
%         filteredData(xx,yy,:) = filter(lpFilt,singlePixelData);
%     end
% end
% close(f)

% filteredData = imageStackMat;

%% Create AVI file frame by frame

if exist(aviFileName,'file')
    delete(aviFileName)
end

aviObject = VideoWriter(aviFileName,'Uncompressed AVI');
open(aviObject);
 
f = waitbar(0,'Creating input AVI file');
for k = 1:tifLength
    waitbar(round(k/tifLength,2),f,'Creating input AVI file');
%     imagesc(filteredData(:,:,k));
%     filteredFrame = getframe(fig,[0.05 0.05 0.9 0.9]);
    writeVideo(aviObject,filteredData{1,k});
end
close(f)
close(aviObject);
end