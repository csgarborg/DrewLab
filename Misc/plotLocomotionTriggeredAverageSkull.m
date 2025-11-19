close all
fileOrDataCell = 2;

motionEventsLocationsX = [];
motionEventsLocationsY = [];
if fileOrDataCell == 1
    dateCell = {'211021','211022','211101','211102','211105','211109','211112','211116','211117','211119','211203','211216','220203','220209','220210','220211','220214','220221','220223','220303','220308','220309','220314','220318','220404','220406','220407','220429','220509','220511','220712','220714','220718','220719','220808','220809','220813','220815','220816','220822','220823','221205','221207','221208','221213'};
    close all
    fileName = {};
    for n = 1:size(dateCell,2)
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
        if size(movementData.motionEvents,1) == 0
            disp('No Ball Motion Events To Plot')
        else
            for n = 1:size(movementData.motionEvents,1)
                motionVectorX = movementData.targetPosition(movementData.motionEvents(n,4):movementData.motionEvents(n,6),1);
                motionVectorY = movementData.targetPosition(movementData.motionEvents(n,4):movementData.motionEvents(n,6),2)*-1;
                if n > 1 || i > 1
                    if length(motionVectorX) > size(motionEventsLocationsX,2)
                        motionVectorX = motionVectorX(1:size(motionEventsLocationsX,2));
                    elseif length(motionVectorX) < size(motionEventsLocationsX,2)
                        motionEventsLocationsX = motionEventsLocationsX(:,1:length(motionVectorX));
                    end
                    if length(motionVectorY) > size(motionEventsLocationsY,2)
                        motionVectorY = motionVectorY(1:size(motionEventsLocationsY,2));
                    elseif length(motionVectorY) < size(motionEventsLocationsY,2)
                        motionEventsLocationsY = motionEventsLocationsY(:,1:length(motionVectorY));
                    end
                end
                motionEventsLocationsX(end+1,:) = medfilt1(motionVectorX-motionVectorX(1),medFiltSize)';
                motionEventsLocationsY(end+1,:) = medfilt1(motionVectorY-motionVectorY(1),medFiltSize)';
            end
        end
        clear movementData
    end
    load(fileName{1});
    [meanX,cIntFillPtsX] = getCIntMeanAndFillPts(motionEventsLocationsX,90);
    [meanY,cIntFillPtsY] = getCIntMeanAndFillPts(motionEventsLocationsY,90);
    meanY = -1*meanY;
    cIntFillPtsY = -1*cIntFillPtsY;
    timeVecX = linspace(round(movementData.motionEvents(1,1)-movementData.motionEvents(1,2)),round(movementData.motionEvents(1,3)-movementData.motionEvents(1,2)),length(meanX));
    timeVecY = linspace(round(movementData.motionEvents(1,1)-movementData.motionEvents(1,2)),round(movementData.motionEvents(1,3)-movementData.motionEvents(1,2)),length(meanY));
    
elseif fileOrDataCell == 2
    load('LTADataCellSkull_FS.mat')
    for n = 1:size(locDataCell)
        if isnan(locDataCell{n,3})
            continue
        end
        motionVectorX = locDataCell{n,3}(2,:);
        motionVectorY = locDataCell{n,4}(2,:);
        if n > 1
            if length(motionVectorX) > size(motionEventsLocationsX,2)
                motionVectorX = motionVectorX(1:size(motionEventsLocationsX,2));
            elseif length(motionVectorX) < size(motionEventsLocationsX,2)
                motionEventsLocationsX = motionEventsLocationsX(:,1:length(motionVectorX));
            end
            if length(motionVectorY) > size(motionEventsLocationsY,2)
                motionVectorY = motionVectorY(1:size(motionEventsLocationsY,2));
            elseif length(motionVectorY) < size(motionEventsLocationsY,2)
                motionEventsLocationsY = motionEventsLocationsY(:,1:length(motionVectorY));
            end
        end
        motionEventsLocationsX(end+1,:) = motionVectorX;
        motionEventsLocationsY(end+1,:) = motionVectorY;
    end
%     load(locDataCell{1,1});
    [meanX,cIntFillPtsX] = getCIntMeanAndFillPts(motionEventsLocationsX,90);
    [meanY,cIntFillPtsY] = getCIntMeanAndFillPts(motionEventsLocationsY,90);
    meanY = -1*meanY;
    cIntFillPtsY = -1*cIntFillPtsY;
    timeVecX = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(meanX));
    timeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(meanY));
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
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[0,1,0,0.1])
end
hold off
title(['\fontsize{20pt}\bf{Mean Motion During Locomotion Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('X Position (\mum)')
%ylim([-maxMeanVal maxMeanVal])
ylim([-5 5])
grid on
subplot(2,1,2)
plot(timeVecY,meanY,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-maxMeanVal maxMeanVal],'r')
plot([0 0],[-1 1],'r')
for n = 1:size(motionEventsLocationsY,1)
    plot(timeVecY,-1*motionEventsLocationsY(n,:),'Color',[0,1,0,0.1])
end
hold off
xlabel('Time (s)')
ylabel('Y Position (\mum)')
%ylim([-maxMeanVal maxMeanVal])
ylim([-5 5])
grid on
clear movementData

h(2) = figure('Color','White');
varX = var(motionEventsLocationsX);
varY = var(motionEventsLocationsY);
subplot(2,1,1)
plot(timeVecX,varX)
xlabel('Time (s)')
ylabel('X Position Variance (\mum)')
title('Variance of motion events')
%ylim([-maxMeanVal maxMeanVal])
ylim([0 2])
grid on
hold on
plot([0 0],[-1 1],'r')
hold off
subplot(2,1,2)
plot(timeVecY,varY)
xlabel('Time (s)')
ylabel('Y Position Variance (\mum)')
%ylim([-maxMeanVal maxMeanVal])
ylim([0 2])
grid on
hold on
plot([0 0],[-1 1],'r')
hold off


stopMotionEventsLocationsX = [];
stopMotionEventsLocationsY = [];
if fileOrDataCell == 1
    dateCell = {'211021','211022','211101','211102','211105','211109','211112','211116','211117','211119','211203','211216','220203','220209','220210','220211','220214','220221','220223','220303','220308','220309','220314','220318','220404','220406','220407','220429','220509','220511','220712','220714','220718','220719','220808','220809','220813','220815','220816','220822','220823','221205','221207','221208','221213'};
    close all
    fileName = {};
    for n = 1:size(dateCell,2)
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
        if size(movementData.stopMotionEvents,1) == 0
            disp('No Ball Motion Events To Plot')
        else
            for n = 1:size(movementData.stopMotionEvents,1)
                stopMotionVectorX = movementData.targetPosition(movementData.stopMotionEvents(n,4):movementData.stopMotionEvents(n,6),1);
                stopMotionVectorY = movementData.targetPosition(movementData.stopMotionEvents(n,4):movementData.stopMotionEvents(n,6),2)*-1;
                if n > 1 || i > 1
                    if length(stopMotionVectorX) > size(stopMotionEventsLocationsX,2)
                        stopMotionVectorX = stopMotionVectorX(1:size(stopMotionEventsLocationsX,2));
                    elseif length(stopMotionVectorX) < size(stopMotionEventsLocationsX,2)
                        stopMotionEventsLocationsX = stopMotionEventsLocationsX(:,1:length(stopMotionVectorX));
                    end
                    if length(stopMotionVectorY) > size(stopMotionEventsLocationsY,2)
                        stopMotionVectorY = stopMotionVectorY(1:size(stopMotionEventsLocationsY,2));
                    elseif length(stopMotionVectorY) < size(stopMotionEventsLocationsY,2)
                        stopMotionEventsLocationsY = stopMotionEventsLocationsY(:,1:length(stopMotionVectorY));
                    end
                end
                stopMotionEventsLocationsX(end+1,:) = medfilt1(stopMotionVectorX-stopMotionVectorX(1),medFiltSize)';
                stopMotionEventsLocationsY(end+1,:) = medfilt1(stopMotionVectorY-stopMotionVectorY(1),medFiltSize)';
            end
        end
        clear movementData
    end
    load(fileName{1});
    [meanX,cIntFillPtsX] = getCIntMeanAndFillPts(stopMotionEventsLocationsX,90);
    [meanY,cIntFillPtsY] = getCIntMeanAndFillPts(stopMotionEventsLocationsY,90);
    meanY = -1*meanY;
    cIntFillPtsY = -1*cIntFillPtsY;
    timeVecX = linspace(round(movementData.stopMotionEvents(1,1)-movementData.stopMotionEvents(1,2)),round(movementData.stopMotionEvents(1,3)-movementData.stopMotionEvents(1,2)),length(meanX));
    timeVecY = linspace(round(movementData.stopMotionEvents(1,1)-movementData.stopMotionEvents(1,2)),round(movementData.stopMotionEvents(1,3)-movementData.stopMotionEvents(1,2)),length(meanY));
    
elseif fileOrDataCell == 2
    load('LTADataCellSkull_FS.mat')
    for n = 1:size(locDataCell)
        if isnan(locDataCell{n,6})
            continue
        end
        stopMotionVectorX = locDataCell{n,6}(2,:);
        stopMotionVectorY = locDataCell{n,7}(2,:);
        if n > 1
            if length(stopMotionVectorX) > size(stopMotionEventsLocationsX,2)
                stopMotionVectorX = stopMotionVectorX(1:size(stopMotionEventsLocationsX,2));
            elseif length(stopMotionVectorX) < size(stopMotionEventsLocationsX,2)
                stopMotionEventsLocationsX = stopMotionEventsLocationsX(:,1:length(stopMotionVectorX));
            end
            if length(stopMotionVectorY) > size(stopMotionEventsLocationsY,2)
                stopMotionVectorY = stopMotionVectorY(1:size(stopMotionEventsLocationsY,2));
            elseif length(stopMotionVectorY) < size(stopMotionEventsLocationsY,2)
                stopMotionEventsLocationsY = stopMotionEventsLocationsY(:,1:length(stopMotionVectorY));
            end
        end
        stopMotionEventsLocationsX(end+1,:) = stopMotionVectorX;
        stopMotionEventsLocationsY(end+1,:) = stopMotionVectorY;
    end
    load(locDataCell{end,1});
    [meanX,cIntFillPtsX] = getCIntMeanAndFillPts(stopMotionEventsLocationsX,90);
    [meanY,cIntFillPtsY] = getCIntMeanAndFillPts(stopMotionEventsLocationsY,90);
    meanY = -1*meanY;
    cIntFillPtsY = -1*cIntFillPtsY;
    timeVecX = linspace(round(locDataCell{1,5}(1,1)-locDataCell{1,5}(1,2)),round(locDataCell{1,5}(1,3)-locDataCell{1,5}(1,2)),length(meanX));
    timeVecY = linspace(round(locDataCell{1,5}(1,1)-locDataCell{1,5}(1,2)),round(locDataCell{1,5}(1,3)-locDataCell{1,5}(1,2)),length(meanY));
end

h(3) = figure('Color','White');
subplot(2,1,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
plot(timeVecX,meanX,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-maxMeanVal maxMeanVal],'r')
plot([0 0],[-1 1],'r')
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,stopMotionEventsLocationsX(n,:),'Color',[0,1,0,0.1])
end
hold off
title(['\fontsize{20pt}\bf{Mean Motion During Stopping Locomotion Events, n = ' num2str(size(stopMotionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('X Position (\mum)')
%ylim([-maxMeanVal maxMeanVal])
ylim([-5 5])
grid on
subplot(2,1,2)
plot(timeVecY,meanY,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-maxMeanVal maxMeanVal],'r')
plot([0 0],[-1 1],'r')
for n = 1:size(motionEventsLocationsY,1)
    plot(timeVecY,-1*stopMotionEventsLocationsY(n,:),'Color',[0,1,0,0.1])
end
hold off
xlabel('Time (s)')
ylabel('Y Position (\mum)')
%ylim([-maxMeanVal maxMeanVal])
ylim([-5 5])
grid on
clear movementData

h(4) = figure('Color','White');
varX = var(stopMotionEventsLocationsX);
varY = var(stopMotionEventsLocationsY);
subplot(2,1,1)
plot(timeVecX,varX)
xlabel('Time (s)')
ylabel('X Position Variance (\mum)')
title('Variance of stop motion events')
%ylim([-maxMeanVal maxMeanVal])
ylim([0 2])
grid on
hold on
plot([0 0],[-1 1],'r')
hold off
subplot(2,1,2)
plot(timeVecY,varY)
xlabel('Time (s)')
ylabel('Y Position Variance (\mum)')
%ylim([-maxMeanVal maxMeanVal])
ylim([0 2])
grid on
hold on
plot([0 0],[-1 1],'r')
hold off