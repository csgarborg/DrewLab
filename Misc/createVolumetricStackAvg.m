tifFileName = 'E:\20-01-15_MouseExp\200115_003.TIF';

tifLength = length(imfinfo(tifFileName));

finalStack = [];
for n = [2 3]
    imgAvg = zeros(512,512);
    avgNum = 0;
    for i = n:2:tifLength
        imgAvg = imgAvg + im2double(medfilt2(imread(tifFileName, i)));
        avgNum = avgNum + 1;
    end
    finalStack(:,:,end+1) = imgAvg ./ avgNum;
    imshow(finalStack(:,:,end));
    pause(1)
end
finalStack(:,:,1) = [];

saveTifStack(finalStack,'E:\20-01-15_MouseExp\200115_003_StackAverage.TIF')