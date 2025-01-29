close all
dateCell = {'240513','240514','240604','240605','240612','240614'};
close all
fileName = {};
fileNameMat = {};
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
        if exist(['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processe_2layerBrainInSkullDataFinal.mat'],'file')
            fileName{end+1} = ['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processe_2layerBrainInSkullDataFinal.mat'];
            fileNameMat{end+1} = [dateCell{n} '_' num2str(i)];
        end
    end
end
medFiltSize = 6;
timeBeforeSqueeze = 2;
timeAfterSqueeze = 6;
motionEventsLocationsX = [];
motionEventsLocationsY = [];
load('quivDataMatSqueeze.mat')
for i = 1:numel(fileName)
    load(fileName{i});
    [squeezeMotionX,squeezeMotionY] = findTriggeredSqueezes(movementData,timeBeforeSqueeze,timeAfterSqueeze);
    targetPositionInSkull = [squeezeMotionX' squeezeMotionY'];
    motionVec = pcaMotionAnalysis(targetPositionInSkull,false,false);
    matInd = strcmp(quivDataMatSqueeze(:,1),fileNameMat{i});
    if sum(matInd) ~= 1
        disp('error')
    end
    quivDataMatSqueeze{matInd,7} = motionVec(1:2);
    if strcmp(fileNameMat{i},'240612_20')
    h(7) = figure('Color','White');
    scatter(targetPositionInSkull(:,1),targetPositionInSkull(:,2),10)
    hold on
    drawArrow([0;0],[motionVec(1);motionVec(2)]);
    hold off
    axis equal square
    axis([-3 3 -3 3])
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
%     title(['\fontsize{20pt}\bf{Position of Brain in Skull}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('\mum')
    ylabel('\mum')
    if movementData.hemisphere == 1
        text(3,0,'Lateral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
        text(-3,0,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
        text(0,3,'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
        text(0,-3,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
    else
        text(3,0,'Medial','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
        text(-3,0,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
        text(0,3,'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
        text(0,-3,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
    end    
    end
    close all
%     motionEventsLocationsX(end+1,:) = squeezeMotionX;
%     motionEventsLocationsY(end+1,:) = squeezeMotionY;
    clear movementData
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [squeezeMotionX,squeezeMotionY] = findTriggeredSqueezes(movementData,timeBeforeSqueeze,timeAfterSqueeze)
[~,locs] = findpeaks(movementData.ballData(:,2),movementData.ballData(:,1),'MinPeakDistance',1.5,'MinPeakHeight',max(movementData.ballData(:,2))-.1);
movementTimeSeconds = (1:length(movementData.targetPosition(:,1)))*movementData.secondsPerFrame;
idxSqueeze = [];
for n = 1:2:length(locs)
    [~,idxSqueeze(end+1)] = min(abs(movementTimeSeconds-locs(n)));
end
squeezeMotionsSum = [];
for n = 1:length(idxSqueeze)
    if idxSqueeze(n)+(6/movementData.secondsPerFrame) > length(movementData.targetPosition)
        disp(['skipped last squeeze event for ' movementData.fileName])
        continue
    end
    squeezeMotionsX(n,:) = movementData.targetPosition(idxSqueeze(n)-floor(timeBeforeSqueeze/movementData.secondsPerFrame):idxSqueeze(n)+floor(timeAfterSqueeze/movementData.secondsPerFrame),1);
    if movementData.hemisphere == 2
        squeezeMotionsX(n,:) = squeezeMotionsX(n,:);
    end
    squeezeMotionsX(n,:) = squeezeMotionsX(n,:) - squeezeMotionsX(n,1);
    
    squeezeMotionsY(n,:) = movementData.targetPosition(idxSqueeze(n)-floor(timeBeforeSqueeze/movementData.secondsPerFrame):idxSqueeze(n)+floor(timeAfterSqueeze/movementData.secondsPerFrame),2);
    squeezeMotionsY(n,:) = squeezeMotionsY(n,:) - squeezeMotionsY(n,1);
end
squeezeMotionX = reshape(squeezeMotionsX,1,[]);
squeezeMotionY = reshape(squeezeMotionsY,1,[]);
end
