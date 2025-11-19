close all
dateCell = {'240513','240514','240604','240605','240612','240614'};
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
timeBeforeSqueeze = 2;
timeAfterSqueeze = 6;
motionEventsLocationsX = [];
motionEventsLocationsY = [];
squeezeDataCell = {};
for i = 1:numel(fileName)
    load(fileName{i});
    [squeezeMotionX,squeezeMotionY,squeezeIdxMat] = findTriggeredSqueezes(movementData,timeBeforeSqueeze,timeAfterSqueeze);
    motionEventsLocationsX(end+1,:) = squeezeMotionX;
    motionEventsLocationsY(end+1,:) = squeezeMotionY;
    timeVec = linspace(timeBeforeSqueeze*-1,timeAfterSqueeze,length(squeezeMotionX));
    squeezeDataCell(i,:) = {fileName,squeezeIdxMat,[timeVec;squeezeMotionX],[timeVec;squeezeMotionY]};
    clear movementData
end
save('squeezeDataCellSkull_FS.mat','squeezeDataCell')
load(fileName{1});
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts(motionEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts(motionEventsLocationsY,90);
timeVecX = linspace(timeBeforeSqueeze*-1,timeAfterSqueeze,length(meanX));
timeVecY = linspace(timeBeforeSqueeze*-1,timeAfterSqueeze,length(meanY));

h(1) = figure('Color','White');
subplot(2,1,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
plot(timeVecX,meanX,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
set(f,'facea',[.2]);
% plot([0 0],[-1 5],'r')
% plot([2 2],[-1 5],'r')
f = fill([0 2 2 0],[4.9 4.9 -.9 -.9],'g','Linestyle','none','FaceAlpha',0.1);
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[0,0,1,0.1])
end
% plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
text(5,-1,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(5,5,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean Motion During Locomotion Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Brian Shift (\mum)')
ylim([-1 5])
xlim([-2 6])
grid on

subplot(2,1,2)
plot(timeVecY,meanY,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
set(f,'facea',[.2]);
% plot([0 0],[-1 5],'r')
% plot([2 2],[-1 5],'r')
f = fill([0 2 2 0],[4.9 4.9 -.9 -.9],'g','Linestyle','none','FaceAlpha',0.1);
for n = 1:size(motionEventsLocationsY,1)
    plot(timeVecY,motionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
end
% plot(timeVecY(brainMotionStart),meanY(brainMotionStart),'rx')
hold off
text(5,-1,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(5,5,'Rostral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
xlabel('Time (s)')
ylabel('\Delta Brian Shift (\mum)')
ylim([-1 5])
xlim([-2 6])
grid on
clear movementData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [squeezeMotionX,squeezeMotionY,squeezeIdxMat] = findTriggeredSqueezes(movementData,timeBeforeSqueeze,timeAfterSqueeze)
[~,locs] = findpeaks(movementData.ballData(:,2),movementData.ballData(:,1),'MinPeakDistance',1.5,'MinPeakHeight',max(movementData.ballData(:,2))-.1);
movementTimeSeconds = (1:length(movementData.targetPosition(:,1)))*movementData.secondsPerFrame;
idxSqueeze = [];
squeezeIdxMat = [];
for n = 1:2:length(locs)
    [~,idxSqueeze(end+1)] = min(abs(movementTimeSeconds-locs(n)));
end
squeezeMotionsSum = [];
for n = 1:length(idxSqueeze)
    if idxSqueeze(n)+(6/movementData.secondsPerFrame) > length(movementData.targetPosition)
        disp(['skipped last squeeze event for ' movementData.fileName])
        continue
    end
    squeezeTimeSec = (idxSqueeze(n)*movementData.secondsPerFrame);
    squeezeIdxMat(n,:) = [squeezeTimeSec-timeBeforeSqueeze squeezeTimeSec squeezeTimeSec+timeAfterSqueeze idxSqueeze(n)-floor(timeBeforeSqueeze/movementData.secondsPerFrame) idxSqueeze(n) idxSqueeze(n)+floor(timeAfterSqueeze/movementData.secondsPerFrame)];
    squeezeMotionsX(n,:) = movementData.targetPosition(idxSqueeze(n)-floor(timeBeforeSqueeze/movementData.secondsPerFrame):idxSqueeze(n)+floor(timeAfterSqueeze/movementData.secondsPerFrame),1);
    if movementData.hemisphere == 2
        squeezeMotionsX(n,:) = squeezeMotionsX(n,:)*-1;
    end
    squeezeMotionsX(n,:) = squeezeMotionsX(n,:) - squeezeMotionsX(n,1);
    
    squeezeMotionsY(n,:) = movementData.targetPosition(idxSqueeze(n)-floor(timeBeforeSqueeze/movementData.secondsPerFrame):idxSqueeze(n)+floor(timeAfterSqueeze/movementData.secondsPerFrame),2);
    squeezeMotionsY(n,:) = squeezeMotionsY(n,:) - squeezeMotionsY(n,1);
end
squeezeMotionX = mean(squeezeMotionsX,1);
squeezeMotionY = mean(squeezeMotionsY,1);
end