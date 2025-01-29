function imageStack = openPupilCamVid(pupilCamFileID, imageHeight, imageWidth)

pixelsPerFrame = imageWidth*imageHeight;
skippedPixels = pixelsPerFrame*2; % Multiply by two because there are 16 bits (2 bytes) per pixel
fid = fopen(pupilCamFileID);
fseek(fid, 0, 'eof');
fileSize = ftell(fid);
fseek(fid, 0, 'bof');
% nFramesToRead = floor(fileSize / (skippedPixels));
nFramesToRead = 1000;
frameInds = 1:nFramesToRead;
imageStack = zeros(imageWidth, imageHeight, nFramesToRead);
for a = 1:nFramesToRead
    disp(['Creating image stack: (' num2str(a) '/' num2str(nFramesToRead) ')']); disp(' ')
    fseek(fid, frameInds(a)*skippedPixels, 'bof');
    z = fread(fid, pixelsPerFrame, '*uint8', 'b');
    img = reshape(z(1:pixelsPerFrame), imageHeight, imageWidth);
    imageStack(:,:,a) = mat2gray(flip(imrotate(img, -90), 2),[0 255]);
end
fclose('all');
end