close all
fileOrDataCell = 1;

EMGEventsLocationsX = [];
EMGEventsLocationsY = [];
if fileOrDataCell == 1
    dateCell = {'211021','211022','211101','211102','211105','211109','211112','211116','211117','211119','211203','211216','220203','220209','220210','220211','220214','220221','220223','220303','220308','220309','220314','220318','220404','220406','220407','220429','220509','220511','220712','220714','220718','220719','220808','220809','220813','220815','220816','220822','220823','221205','221207','221208','221213'};
    close all
    fileName = {};
    for n = 12:size(dateCell,2)
        folderName = ['I:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/'];
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
            if exist(['I:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processe_2layerBrainInSkullDataFinal.mat'],'file')
                fileName{end+1} = ['I:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processe_2layerBrainInSkullDataFinal.mat'];
            end
        end
    end
    for i = 1:numel(fileName)
        load(fileName{i});
        if size(movementData.EMGEvents,1) == 0
            disp('No Ball EMG Events To Plot')
        else
            for n = 1:size(movementData.EMGEvents,1)
                EMGVectorX = movementData.targetPosition(movementData.EMGEvents(n,4):movementData.EMGEvents(n,6),1);
                EMGVectorY = movementData.targetPosition(movementData.EMGEvents(n,4):movementData.EMGEvents(n,6),2)*-1;
                if n > 1 || i > 1
                    if length(EMGVectorX) > size(EMGEventsLocationsX,2)
                        EMGVectorX = EMGVectorX(1:size(EMGEventsLocationsX,2));
                    elseif length(EMGVectorX) < size(EMGEventsLocationsX,2)
                        EMGEventsLocationsX = EMGEventsLocationsX(:,1:length(EMGVectorX));
                    end
                    if length(EMGVectorY) > size(EMGEventsLocationsY,2)
                        EMGVectorY = EMGVectorY(1:size(EMGEventsLocationsY,2));
                    elseif length(EMGVectorY) < size(EMGEventsLocationsY,2)
                        EMGEventsLocationsY = EMGEventsLocationsY(:,1:length(EMGVectorY));
                    end
                end
                EMGEventsLocationsX(end+1,:) = medfilt1(EMGVectorX-EMGVectorX(1),medFiltSize)';
                EMGEventsLocationsY(end+1,:) = medfilt1(EMGVectorY-EMGVectorY(1),medFiltSize)';
            end
        end
        clear movementData
    end
    load(fileName{1});
    [meanX,cIntFillPtsX] = getCIntMeanAndFillPts(EMGEventsLocationsX,90);
    [meanY,cIntFillPtsY] = getCIntMeanAndFillPts(EMGEventsLocationsY,90);
    meanY = -1*meanY;
    cIntFillPtsY = -1*cIntFillPtsY;
    timeVecX = linspace(round(movementData.EMGEvents(1,1)-movementData.EMGEvents(1,2)),round(movementData.EMGEvents(1,3)-movementData.EMGEvents(1,2)),length(meanX));
    timeVecY = linspace(round(movementData.EMGEvents(1,1)-movementData.EMGEvents(1,2)),round(movementData.EMGEvents(1,3)-movementData.EMGEvents(1,2)),length(meanY));
    
elseif fileOrDataCell == 2
    load('LTADataCell.mat')
    for n = 1:size(dataCellSaveNew)
        if isnan(dataCellSaveNew{n,3})
            continue
        end
        EMGVectorX = dataCellSaveNew{n,3}(2,:);
        EMGVectorY = dataCellSaveNew{n,4}(2,:);
        if n > 1
            if length(EMGVectorX) > size(EMGEventsLocationsX,2)
                EMGVectorX = EMGVectorX(1:size(EMGEventsLocationsX,2));
            elseif length(EMGVectorX) < size(EMGEventsLocationsX,2)
                EMGEventsLocationsX = EMGEventsLocationsX(:,1:length(EMGVectorX));
            end
            if length(EMGVectorY) > size(EMGEventsLocationsY,2)
                EMGVectorY = EMGVectorY(1:size(EMGEventsLocationsY,2));
            elseif length(EMGVectorY) < size(EMGEventsLocationsY,2)
                EMGEventsLocationsY = EMGEventsLocationsY(:,1:length(EMGVectorY));
            end
        end
        EMGEventsLocationsX(end+1,:) = EMGVectorX;
        EMGEventsLocationsY(end+1,:) = EMGVectorY;
    end
    load(dataCellSaveNew{1,1});
    [meanX,cIntFillPtsX] = getCIntMeanAndFillPts(EMGEventsLocationsX,90);
    [meanY,cIntFillPtsY] = getCIntMeanAndFillPts(EMGEventsLocationsY,90);
    meanY = -1*meanY;
    cIntFillPtsY = -1*cIntFillPtsY;
    timeVecX = linspace(round(movementData.EMGEvents(1,1)-movementData.EMGEvents(1,2)),round(movementData.EMGEvents(1,3)-movementData.EMGEvents(1,2)),length(meanX));
    timeVecY = linspace(round(movementData.EMGEvents(1,1)-movementData.EMGEvents(1,2)),round(movementData.EMGEvents(1,3)-movementData.EMGEvents(1,2)),length(meanY));
end

h(1) = figure('Color','White');
subplot(2,1,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
plot(timeVecX,meanX,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-maxMeanVal maxMeanVal],'r')
plot([0 0],[-1 1],'r')
for n = 1:size(EMGEventsLocationsX,1)
    plot(timeVecX,EMGEventsLocationsX(n,:),'Color',[0,1,0,0.1])
end
hold off
title(['\fontsize{20pt}\bf{Mean EMG During LocoEMG Events, n = ' num2str(size(EMGEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('X Position (\mum)')
%ylim([-maxMeanVal maxMeanVal])
ylim([-1 3])
grid on
subplot(2,1,2)
plot(timeVecY,meanY,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-maxMeanVal maxMeanVal],'r')
plot([0 0],[-1 1],'r')
for n = 1:size(EMGEventsLocationsY,1)
    plot(timeVecY,-1*EMGEventsLocationsY(n,:),'Color',[0,1,0,0.1])
end
hold off
xlabel('Time (s)')
ylabel('Y Position (\mum)')
%ylim([-maxMeanVal maxMeanVal])
ylim([-1 3])
grid on
clear movementData


stopEMGEventsLocationsX = [];
stopEMGEventsLocationsY = [];
if fileOrDataCell == 1
    dateCell = {'211021','211022','211101','211102','211105','211109','211112','211116','211117','211119','211203','211216','220203','220209','220210','220211','220214','220221','220223','220303','220308','220309','220314','220318','220404','220406','220407','220429','220509','220511','220712','220714','220718','220719','220808','220809','220813','220815','220816','220822','220823','221205','221207','221208','221213'};
    close all
    fileName = {};
    for n = 12:size(dateCell,2)
        folderName = ['I:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/'];
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
            if exist(['I:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processe_2layerBrainInSkullDataFinal.mat'],'file')
                fileName{end+1} = ['I:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processe_2layerBrainInSkullDataFinal.mat'];
            end
        end
    end
    for i = 1:numel(fileName)
        load(fileName{i});
        if size(movementData.stopEMGEvents,1) == 0
            disp('No Ball EMG Events To Plot')
        else
            for n = 1:size(movementData.stopEMGEvents,1)
                stopEMGVectorX = movementData.targetPosition(movementData.stopEMGEvents(n,4):movementData.stopEMGEvents(n,6),1);
                stopEMGVectorY = movementData.targetPosition(movementData.stopEMGEvents(n,4):movementData.stopEMGEvents(n,6),2)*-1;
                if n > 1 || i > 1
                    if length(stopEMGVectorX) > size(stopEMGEventsLocationsX,2)
                        stopEMGVectorX = stopEMGVectorX(1:size(stopEMGEventsLocationsX,2));
                    elseif length(stopEMGVectorX) < size(stopEMGEventsLocationsX,2)
                        stopEMGEventsLocationsX = stopEMGEventsLocationsX(:,1:length(stopEMGVectorX));
                    end
                    if length(stopEMGVectorY) > size(stopEMGEventsLocationsY,2)
                        stopEMGVectorY = stopEMGVectorY(1:size(stopEMGEventsLocationsY,2));
                    elseif length(stopEMGVectorY) < size(stopEMGEventsLocationsY,2)
                        stopEMGEventsLocationsY = stopEMGEventsLocationsY(:,1:length(stopEMGVectorY));
                    end
                end
                stopEMGEventsLocationsX(end+1,:) = medfilt1(stopEMGVectorX-stopEMGVectorX(1),medFiltSize)';
                stopEMGEventsLocationsY(end+1,:) = medfilt1(stopEMGVectorY-stopEMGVectorY(1),medFiltSize)';
            end
        end
        clear movementData
    end
    load(fileName{1});
    [meanX,cIntFillPtsX] = getCIntMeanAndFillPts(stopEMGEventsLocationsX,90);
    [meanY,cIntFillPtsY] = getCIntMeanAndFillPts(stopEMGEventsLocationsY,90);
    meanY = -1*meanY;
    cIntFillPtsY = -1*cIntFillPtsY;
    timeVecX = linspace(round(movementData.stopEMGEvents(1,1)-movementData.stopEMGEvents(1,2)),round(movementData.stopEMGEvents(1,3)-movementData.stopEMGEvents(1,2)),length(meanX));
    timeVecY = linspace(round(movementData.stopEMGEvents(1,1)-movementData.stopEMGEvents(1,2)),round(movementData.stopEMGEvents(1,3)-movementData.stopEMGEvents(1,2)),length(meanY));
    
elseif fileOrDataCell == 2
    load('LTADataCell.mat')
    for n = 1:size(dataCellSaveNew)
        if isnan(dataCellSaveNew{n,6})
            continue
        end
        stopEMGVectorX = dataCellSaveNew{n,6}(2,:);
        stopEMGVectorY = dataCellSaveNew{n,7}(2,:);
        if n > 1
            if length(stopEMGVectorX) > size(stopEMGEventsLocationsX,2)
                stopEMGVectorX = stopEMGVectorX(1:size(stopEMGEventsLocationsX,2));
            elseif length(stopEMGVectorX) < size(stopEMGEventsLocationsX,2)
                stopEMGEventsLocationsX = stopEMGEventsLocationsX(:,1:length(stopEMGVectorX));
            end
            if length(stopEMGVectorY) > size(stopEMGEventsLocationsY,2)
                stopEMGVectorY = stopEMGVectorY(1:size(stopEMGEventsLocationsY,2));
            elseif length(stopEMGVectorY) < size(stopEMGEventsLocationsY,2)
                stopEMGEventsLocationsY = stopEMGEventsLocationsY(:,1:length(stopEMGVectorY));
            end
        end
        stopEMGEventsLocationsX(end+1,:) = stopEMGVectorX;
        stopEMGEventsLocationsY(end+1,:) = stopEMGVectorY;
    end
    load(dataCellSaveNew{end,1});
    [meanX,cIntFillPtsX] = getCIntMeanAndFillPts(stopEMGEventsLocationsX,90);
    [meanY,cIntFillPtsY] = getCIntMeanAndFillPts(stopEMGEventsLocationsY,90);
    meanY = -1*meanY;
    cIntFillPtsY = -1*cIntFillPtsY;
    timeVecX = linspace(round(movementData.stopEMGEvents(1,1)-movementData.stopEMGEvents(1,2)),round(movementData.stopEMGEvents(1,3)-movementData.stopEMGEvents(1,2)),length(meanX));
    timeVecY = linspace(round(movementData.stopEMGEvents(1,1)-movementData.stopEMGEvents(1,2)),round(movementData.stopEMGEvents(1,3)-movementData.stopEMGEvents(1,2)),length(meanY));
end

h(2) = figure('Color','White');
subplot(2,1,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
plot(timeVecX,meanX,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-maxMeanVal maxMeanVal],'r')
plot([0 0],[-1 1],'r')
for n = 1:size(EMGEventsLocationsX,1)
    plot(timeVecX,stopEMGEventsLocationsX(n,:),'Color',[0,1,0,0.1])
end
hold off
title(['\fontsize{20pt}\bf{Mean EMG During Stopping LocoEMG Events, n = ' num2str(size(stopEMGEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('X Position (\mum)')
%ylim([-maxMeanVal maxMeanVal])
ylim([-3 1])
grid on
subplot(2,1,2)
plot(timeVecY,meanY,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-maxMeanVal maxMeanVal],'r')
plot([0 0],[-1 1],'r')
for n = 1:size(EMGEventsLocationsY,1)
    plot(timeVecY,-1*stopEMGEventsLocationsY(n,:),'Color',[0,1,0,0.1])
end
hold off
xlabel('Time (s)')
ylabel('Y Position (\mum)')
%ylim([-maxMeanVal maxMeanVal])
ylim([-3 1])
grid on
clear movementData