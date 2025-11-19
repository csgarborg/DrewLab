close all
load('movementDataLog.mat')
binSize = 150;

latMedX = cell2mat(moveDataMat(:,4));
latMedY = [];
for n = 1:size(moveDataMat,1)
    latMedY(n,1) = moveDataMat{n,7}(1);
end
latMedX(cell2mat(moveDataMat(:,5))<-2000) = [];
latMedY(cell2mat(moveDataMat(:,5))<-2000) = [];
edgeX = floor(min(latMedX)):binSize:(ceil(max(latMedX)/binSize)*binSize)+binSize;
binX = discretize(latMedX,edgeX);
latMedCell = cell(size(edgeX,2)-1,1);
for n = 1:size(latMedY)
    latMedCell{binX(n),1}(end+1) = latMedY(n);
end
latMedAvg = [];
for n = 1:size(latMedCell)
    if isempty(latMedCell{n,1})
        latMedAvg(n,1) = 0;
    else
        latMedAvg(n,1) = mean(latMedCell{n,1});
    end
end
% [latMedAvg,latMedMax,latMedStd] = grpstats(latMedY,binX,{@mean,@max,@std});
binCentX = edgeX(2:end)-(edgeX(2)-edgeX(1))/2;

rosCauX = cell2mat(moveDataMat(:,5));
rosCauY = [];
for n = 1:size(moveDataMat,1)
    rosCauY(n,1) = moveDataMat{n,7}(2);
end
rosCauX(cell2mat(moveDataMat(:,5))<-2000) = [];
rosCauY(cell2mat(moveDataMat(:,5))<-2000) = [];
edgeY = floor(min(rosCauX)):binSize:(ceil(max(rosCauX)/binSize)*binSize)+binSize;
binY = discretize(rosCauX,edgeY);
rosCauCell = cell(size(edgeY,2)-1,1);
for n = 1:size(rosCauY)
    rosCauCell{binY(n),1}(end+1) = rosCauY(n);
end
rosCauAvg = [];
for n = 1:size(rosCauCell)
    if isempty(rosCauCell{n,1})
        rosCauAvg(n,1) = 0;
    else
        rosCauAvg(n,1) = mean(rosCauCell{n,1});
    end
end
% [rosCauAvg,rosCauMax,rosCauStd] = grpstats(rosCauY,binY,{@mean,@max,@std});
binCentY = edgeY(2:end)-(edgeY(2)-edgeY(1))/2;

figure(1)
bar(binCentX,latMedAvg)
title('lateral/medial movement vs x location')
xlim([-3000 3000])

figure(2)
bar(binCentY,rosCauAvg)
title('rostral/caudal movement vs y location')
xlim([-1500 4000])
ylim([-1 5])


clear
load('quivDataMatSqueeze.mat')
binSize = 150;

latMedX = cell2mat(quivDataMatSqueeze(:,4));
latMedY = [];
for n = 1:size(quivDataMatSqueeze,1)
    latMedY(n,1) = quivDataMatSqueeze{n,7}(1);
end
latMedX(cell2mat(quivDataMatSqueeze(:,5))<-2000) = [];
latMedY(cell2mat(quivDataMatSqueeze(:,5))<-2000) = [];
edgeX = floor(min(latMedX)):binSize:(ceil(max(latMedX)/binSize)*binSize)+binSize;
binX = discretize(latMedX,edgeX);
latMedCell = cell(size(edgeX,2)-1,1);
for n = 1:size(latMedY)
    latMedCell{binX(n),1}(end+1) = latMedY(n);
end
latMedAvg = [];
for n = 1:size(latMedCell)
    if isempty(latMedCell{n,1})
        latMedAvg(n,1) = 0;
    else
        latMedAvg(n,1) = mean(latMedCell{n,1});
    end
end
% [latMedAvg,latMedMax,latMedStd] = grpstats(latMedY,binX,{@mean,@max,@std});
binCentX = edgeX(2:end)-(edgeX(2)-edgeX(1))/2;

rosCauX = cell2mat(quivDataMatSqueeze(:,5));
rosCauY = [];
for n = 1:size(quivDataMatSqueeze,1)
    rosCauY(n,1) = quivDataMatSqueeze{n,7}(2);
end
rosCauX(cell2mat(quivDataMatSqueeze(:,5))<-2000) = [];
rosCauY(cell2mat(quivDataMatSqueeze(:,5))<-2000) = [];
edgeY = floor(min(rosCauX)):binSize:(ceil(max(rosCauX)/binSize)*binSize)+binSize;
binY = discretize(rosCauX,edgeY);
rosCauCell = cell(size(edgeY,2)-1,1);
for n = 1:size(rosCauY)
    rosCauCell{binY(n),1}(end+1) = rosCauY(n);
end
rosCauAvg = [];
for n = 1:size(rosCauCell)
    if isempty(rosCauCell{n,1})
        rosCauAvg(n,1) = 0;
    else
        rosCauAvg(n,1) = mean(rosCauCell{n,1});
    end
end
% [rosCauAvg,rosCauMax,rosCauStd] = grpstats(rosCauY,binY,{@mean,@max,@std});
binCentY = edgeY(2:end)-(edgeY(2)-edgeY(1))/2;

figure(3)
bar(binCentX,latMedAvg)
title('lateral/medial movement vs x location')
xlim([-3000 3000])

figure(4)
bar(binCentY,rosCauAvg)
title('rostral/caudal movement vs y location')
xlim([-1500 4000])
ylim([-1 5])