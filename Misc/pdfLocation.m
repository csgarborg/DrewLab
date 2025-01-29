close all
fileList = {'F:\19-07-30_PaperExp\190730_002_processed.mat'};
% fileList = {'F:\19-07-30_PaperExp\190730_001_processed.mat','F:\19-07-30_PaperExp\190730_002_processed.mat'};

allPosDataX = [];
allPosDataY = [];
for n = 1:numel(fileList)
    load(fileList{n});
    allPosDataX = vertcat(allPosDataX,movementData.targetPosition(:,1));
    allPosDataY = vertcat(allPosDataY,movementData.targetPosition(:,2));
end

[a,b] = ksdensity(allPosDataX);
[c,d] = ksdensity(allPosDataY);

figure(1)
plot(b,a)

figure(2)
plot(d,c)