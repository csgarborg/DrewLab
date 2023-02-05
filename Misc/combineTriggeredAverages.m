close all
% fileNames = {'H:\21-10-21_MouseExp\211021_00','H:\21-10-22_MouseExp\211022_00','H:\21-11-01_MouseExp\211101_00','H:\21-11-09_MouseExp\211109_00','H:\21-11-12_MouseExp\211112_00','H:\21-11-16_MouseExp\211116_00','H:\21-11-17_MouseExp\211117_00','H:\21-11-19_MouseExp\211119_00','H:\21-12-03_MouseExp\211203_00'};
% vecs = {1:9,1:12,1:8,[1 2 3 5 6 7 8 9],[1:7 9:18],1:9,[1 2 3 7 10 11 13 14 15 16 17 18],[1:18 20:24],1:18};
% fileName = {};
% for j = 1:numel(fileNames)
%     for k = vecs{j}
%         if k <= 9
%             fileName{end+1} = [fileNames{j} num2str(k) '_processe_2layerBrainInSkullDataFinal.mat'];
%         else
%             fileName{end+1} = [fileNames{j}(1:end-1) num2str(k) '_processe_2layerBrainInSkullDataFinal.mat'];
%         end
%     end
% end
% fileName = {'C:\Workspace\210218_005_processed_combined.mat','C:\Workspace\210218_006_processed_combinedMean.mat','C:\Workspace\210218_007_processed_combined.mat','C:\Workspace\210218_008_processed_combined.mat'};
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
medFiltSize = 6;

motionEventsLocationsX = [];
motionEventsLocationsY = [];
for i = 1:numel(fileName)
    load(fileName{i});
    if size(movementData.motionEvents,1) == 0
        disp('No Ball Motion Events To Plot')
    else
        for n = 1:size(movementData.motionEvents,1)
            motionVectorX = movementData.targetPosition(movementData.motionEvents(n,4):movementData.motionEvents(n,6),1);
            motionVectorY = movementData.targetPosition(movementData.motionEvents(n,4):movementData.motionEvents(n,6),2);
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
timeVecX = linspace(round(movementData.motionEvents(1,1)-movementData.motionEvents(1,2)),round(movementData.motionEvents(1,3)-movementData.motionEvents(1,2)),length(meanX));
timeVecY = linspace(round(movementData.motionEvents(1,1)-movementData.motionEvents(1,2)),round(movementData.motionEvents(1,3)-movementData.motionEvents(1,2)),length(meanY));

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
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[1,0,0,0.2])
end
hold off
title(['\fontsize{20pt}\bf{Mean Motion During Locomotion Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('X Position (\mum)')
%ylim([-maxMeanVal maxMeanVal])
ylim([-1 1])
grid on
subplot(2,1,2)
plot(timeVecY,meanY,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-maxMeanVal maxMeanVal],'r')
plot([0 0],[-1 1],'r')
hold off
xlabel('Time (s)')
ylabel('Y Position (\mum)')
%ylim([-maxMeanVal maxMeanVal])
ylim([-1 1])
grid on
clear movementData

h(2) = figure('Color','White');
EMGEventsLocationsX = [];
EMGEventsLocationsY = [];
for i = 1:numel(fileName)
    load(fileName{i});
    if size(movementData.EMGEvents,1) == 0
        disp('No EMG Events To Plot')
    else
        for n = 1:size(movementData.EMGEvents,1)
            motionVectorX = movementData.targetPosition(movementData.EMGEvents(n,4):movementData.EMGEvents(n,6),1);
            motionVectorY = movementData.targetPosition(movementData.EMGEvents(n,4):movementData.EMGEvents(n,6),2);
            if n > 1 || i > 1
                if length(motionVectorX) > size(EMGEventsLocationsX,2)
                    motionVectorX = motionVectorX(1:size(EMGEventsLocationsX,2));
                elseif length(motionVectorX) < size(EMGEventsLocationsX,2)
                    EMGEventsLocationsX = EMGEventsLocationsX(:,1:length(motionVectorX));
                end
                if length(motionVectorY) > size(EMGEventsLocationsY,2)
                    motionVectorY = motionVectorY(1:size(EMGEventsLocationsY,2));
                elseif length(motionVectorY) < size(EMGEventsLocationsY,2)
                    EMGEventsLocationsY = EMGEventsLocationsY(:,1:length(motionVectorY));
                end
            end
            EMGEventsLocationsX(end+1,:) = medfilt1(motionVectorX-motionVectorX(1),medFiltSize);
            EMGEventsLocationsY(end+1,:) = medfilt1(motionVectorY-motionVectorY(1),medFiltSize);
        end
    end
    clear movementData
end
load(fileName{1});
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts(EMGEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts(EMGEventsLocationsY,90);
timeVecX = linspace(round(movementData.EMGEvents(1,1)-movementData.EMGEvents(1,2)),round(movementData.EMGEvents(1,3)-movementData.EMGEvents(1,2)),length(meanX));
timeVecY = linspace(round(movementData.EMGEvents(1,1)-movementData.EMGEvents(1,2)),round(movementData.EMGEvents(1,3)-movementData.EMGEvents(1,2)),length(meanY));
subplot(2,1,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
plot(timeVecX,meanX,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-maxMeanVal maxMeanVal],'r')
plot([0 0],[-1 1],'r')
hold off
title(['\fontsize{20pt}\bf{Mean Motion During EMG Events, n = ' num2str(size(EMGEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('X Position (\mum)')
%ylim([-maxMeanVal maxMeanVal])
ylim([-1 1])
grid on
subplot(2,1,2)
plot(timeVecY,meanY,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-maxMeanVal maxMeanVal],'r')
plot([0 0],[-1 1],'r')
hold off
xlabel('Time (s)')
ylabel('Y Position (\mum)')
%ylim([-maxMeanVal maxMeanVal])
ylim([-1 1])
grid on
clear movementData

stopMotionEventsLocationsX = [];
stopMotionEventsLocationsY = [];
for i = 1:numel(fileName)
    load(fileName{i});
    if size(movementData.stopMotionEvents,1) == 0
        disp('No Ball Motion Events To Plot')
    else
        for n = 1:size(movementData.stopMotionEvents,1)
            motionVectorX = movementData.targetPosition(movementData.stopMotionEvents(n,4):movementData.stopMotionEvents(n,6),1);
            motionVectorY = movementData.targetPosition(movementData.stopMotionEvents(n,4):movementData.stopMotionEvents(n,6),2);
            if ~isempty(stopMotionEventsLocationsX)
                if length(motionVectorX) > size(stopMotionEventsLocationsX,2)
                    motionVectorX = motionVectorX(1:size(stopMotionEventsLocationsX,2));
                elseif length(motionVectorX) < size(stopMotionEventsLocationsX,2)
                    stopMotionEventsLocationsX = stopMotionEventsLocationsX(:,1:length(motionVectorX));
                end
                if length(motionVectorY) > size(stopMotionEventsLocationsY,2)
                    motionVectorY = motionVectorY(1:size(stopMotionEventsLocationsY,2));
                elseif length(motionVectorY) < size(stopMotionEventsLocationsY,2)
                    stopMotionEventsLocationsY = stopMotionEventsLocationsY(:,1:length(motionVectorY));
                end
            end
            stopMotionEventsLocationsX(end+1,:) = medfilt1(motionVectorX-motionVectorX(1),medFiltSize)';
            stopMotionEventsLocationsY(end+1,:) = medfilt1(motionVectorY-motionVectorY(1),medFiltSize)';
        end
    end
    clear movementData
end
load(fileName{1});
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts(stopMotionEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts(stopMotionEventsLocationsY,90);
timeVecX = linspace(round(movementData.motionEvents(1,1)-movementData.motionEvents(1,2)),round(movementData.motionEvents(1,3)-movementData.motionEvents(1,2)),length(meanX));
timeVecY = linspace(round(movementData.motionEvents(1,1)-movementData.motionEvents(1,2)),round(movementData.motionEvents(1,3)-movementData.motionEvents(1,2)),length(meanY));

h(3) = figure('Color','White');
subplot(2,1,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
plot(timeVecX,meanX,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-maxMeanVal maxMeanVal],'r')
plot([0 0],[-1 1],'r')
hold off
title(['\fontsize{20pt}\bf{Mean Motion During Stopping Locomotion Events, n = ' num2str(size(stopMotionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('X Position (\mum)')
%ylim([-maxMeanVal maxMeanVal])
ylim([-1 1])
grid on
subplot(2,1,2)
plot(timeVecY,meanY,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-maxMeanVal maxMeanVal],'r')
plot([0 0],[-1 1],'r')
hold off
xlabel('Time (s)')
ylabel('Y Position (\mum)')
%ylim([-maxMeanVal maxMeanVal])
ylim([-1 1])
grid on
clear movementData

h(4) = figure('Color','White');
stopEMGEventsLocationsX = [];
stopEMGEventsLocationsY = [];
for i = 1:numel(fileName)
    load(fileName{i});
    if size(movementData.stopEMGEvents,1) == 0
        disp('No EMG Events To Plot')
    else
        for n = 1:size(movementData.stopEMGEvents,1)
            motionVectorX = movementData.targetPosition(movementData.stopEMGEvents(n,4):movementData.stopEMGEvents(n,6),1);
            motionVectorY = movementData.targetPosition(movementData.stopEMGEvents(n,4):movementData.stopEMGEvents(n,6),2);
            if ~isempty(stopEMGEventsLocationsY)
                if length(motionVectorX) > size(stopEMGEventsLocationsX,2)
                    motionVectorX = motionVectorX(1:size(stopEMGEventsLocationsX,2));
                elseif length(motionVectorX) < size(stopEMGEventsLocationsX,2)
                    stopEMGEventsLocationsX = stopEMGEventsLocationsX(:,1:length(motionVectorX));
                end
                if length(motionVectorY) > size(stopEMGEventsLocationsY,2)
                    motionVectorY = motionVectorY(1:size(stopEMGEventsLocationsY,2));
                elseif length(motionVectorY) < size(stopEMGEventsLocationsY,2)
                    stopEMGEventsLocationsY = stopEMGEventsLocationsY(:,1:length(motionVectorY));
                end
            end
            stopEMGEventsLocationsX(end+1,:) = medfilt1(motionVectorX-motionVectorX(1),medFiltSize);
            stopEMGEventsLocationsY(end+1,:) = medfilt1(motionVectorY-motionVectorY(1),medFiltSize);
        end
    end
    clear movementData
end
load(fileName{1});
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts(stopEMGEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts(stopEMGEventsLocationsY,90);
timeVecX = linspace(round(movementData.EMGEvents(1,1)-movementData.EMGEvents(1,2)),round(movementData.EMGEvents(1,3)-movementData.EMGEvents(1,2)),length(meanX));
timeVecY = linspace(round(movementData.EMGEvents(1,1)-movementData.EMGEvents(1,2)),round(movementData.EMGEvents(1,3)-movementData.EMGEvents(1,2)),length(meanY));
subplot(2,1,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
plot(timeVecX,meanX,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-maxMeanVal maxMeanVal],'r')
plot([0 0],[-1 1],'r')
hold off
title(['\fontsize{20pt}\bf{Mean Motion During Stopping EMG Events, n = ' num2str(size(stopEMGEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('X Position (\mum)')
%ylim([-maxMeanVal maxMeanVal])
ylim([-1 1])
grid on
subplot(2,1,2)
plot(timeVecY,meanY,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-maxMeanVal maxMeanVal],'r')
plot([0 0],[-1 1],'r')
hold off
xlabel('Time (s)')
ylabel('Y Position (\mum)')
%ylim([-maxMeanVal maxMeanVal])
ylim([-1 1])
grid on
clear movementData