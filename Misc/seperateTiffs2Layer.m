tifFileName = 'C:\Workspace\21-06-07_MouseExp\210607_006.TIF';

tifLength = length(imfinfo(tifFileName));

for n = [2 3]
    finalStack = [];
    imgAvg = zeros(512,512);
    for i = n:2:tifLength
        finalStack(:,:,end+1) = im2double(imread(tifFileName, i));
    end
    finalStack(:,:,1) = [];
    if n == 2
        layer1 = finalStack;
    else
        layer2 = finalStack;
    end
end

shortStackLength = min([size(layer1,3),size(layer2,3)]);
layer1 = layer1(:,:,1:shortStackLength);
layer2 = layer2(:,:,1:shortStackLength);

saveTifStack(layer1,[tifFileName(1:end-4) '_Layer1.tif'])
saveTifStack(layer2,[tifFileName(1:end-4) '_Layer2.tif'])