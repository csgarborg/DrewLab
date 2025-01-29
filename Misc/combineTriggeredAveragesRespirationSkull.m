close all
dateCell = {'240611'};
close all
fileName = {};
for n = 1:size(dateCell,2)
    folderName = ['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/'];
    fileList = dir(folderName);
    fileNamesCell = struct2cell(fileList);
    fileNames = fileNamesCell(1,:);
    maxRun = 0;
    for i = 1:size(fileNames,2)
        if contains(fileNames{i},dateCell{n}) && str2double(fileNames{i}(8:10)) > maxRun
             maxRun = str2double(fileNames{i}(8:10));
        end
    end
    for i = 1:maxRun
        if i > 9
            runNumberStr = num2str(i);
        else
            runNumberStr = ['0' num2str(i)];
        end
        if exist(['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer2_combined.mat'],'file')
            fileName{end+1} = ['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer2_combined.mat'];
        end
    end
end
medFiltSize = 6;
timeBeforeResp = 0.25;
timeAfterResp = 1;
load(fileName{i});
[respMotionX,respMotionY,respIdxMat] = findTriggeredResps(movementData,timeBeforeResp,timeAfterResp);
timeVec = linspace(timeBeforeResp*-1,timeAfterResp,length(respMotionX(1,:)));
for i = 1:size(respIdxMat,1)
    respDataCell(i,:) = {fileName,respIdxMat(i,:),[timeVec;respMotionX(i,:)],[timeVec;respMotionY(i,:)*-1]};
end
save('respDataCellSkull_FS.mat','respDataCell')
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts(respMotionX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts(respMotionY,90);
timeVecX = linspace(timeBeforeResp*-1,timeAfterResp,length(meanX));
timeVecY = linspace(timeBeforeResp*-1,timeAfterResp,length(meanY));

h(1) = figure('Color','White');
subplot(2,1,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
plot(timeVecX,meanX,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
set(f,'facea',[.2]);
plot([0 0],[-1 5],'r')
for n = 1:size(respMotionX,1)
    plot(timeVecX,respMotionX(n,:),'Color',[0,0,1,0.1])
end
% plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
text(1,-1,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(1,2,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
title(['suppFig' 10 '\fontsize{20pt}\bf{Mean Motion During Respiration EMG Events, n = ' num2str(size(respMotionX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Brian Shift (\mum)')
ylim([-1 2])
xlim([-1*timeBeforeResp timeAfterResp])
grid on

subplot(2,1,2)
plot(timeVecY,meanY,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
set(f,'facea',[.2]);
plot([0 0],[-1 5],'r')
for n = 1:size(respMotionY,1)
    plot(timeVecY,respMotionY(n,:),'Color',[0,0,1,0.1])
end
% plot(timeVecY(brainMotionStart),meanY(brainMotionStart),'rx')
hold off
text(1,-1,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(1,2,'Rostral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
xlabel('Time (s)')
ylabel('\Delta Brian Shift (\mum)')
ylim([-1 2])
xlim([-1*timeBeforeResp timeAfterResp])
grid on
clear movementData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [respMotionsX,respMotionsY,respIdxMat] = findTriggeredResps(movementData,timeBeforeResp,timeAfterResp)
[pksResp,locs] = findpeaks(movementData.emgData(:,2),movementData.emgData(:,1),'MinPeakDistance',.8,'MinPeakHeight',1.25);
pksResp(locs<50 | locs>180) = [];
locs(locs<50 | locs>180) = [];
% close all;plot(movementData.emgData(:,1),movementData.emgData(:,2));hold on;plot(locs,pksResp,'o');hold off;xlim([50 180]);
movementTimeSeconds = (1:length(movementData.targetPosition(:,1)))*movementData.secondsPerFrame;
idxResp = [];
respIdxMat = [];
for n = 1:length(locs)
    [~,idxResp(end+1)] = min(abs(movementTimeSeconds-locs(n)));
end
RespMotionsSum = [];
for n = 1:length(idxResp)
    respTimeSec = (idxResp(n)*movementData.secondsPerFrame);
    respIdxMat(n,:) = [respTimeSec-timeBeforeResp respTimeSec respTimeSec+timeAfterResp idxResp(n)-floor(timeBeforeResp/movementData.secondsPerFrame) idxResp(n) idxResp(n)+floor(timeAfterResp/movementData.secondsPerFrame)];
    respMotionsX(n,:) = movementData.targetPosition(idxResp(n)-floor(timeBeforeResp/movementData.secondsPerFrame):idxResp(n)+floor(timeAfterResp/movementData.secondsPerFrame),1);
    if movementData.hemisphere == 2
        respMotionsX(n,:) = respMotionsX(n,:)*-1;
    end
    respMotionsX(n,:) = respMotionsX(n,:) - respMotionsX(n,1);
    
    respMotionsY(n,:) = movementData.targetPosition(idxResp(n)-floor(timeBeforeResp/movementData.secondsPerFrame):idxResp(n)+floor(timeAfterResp/movementData.secondsPerFrame),2);
    respMotionsY(n,:) = respMotionsY(n,:) - respMotionsY(n,1);
end
end