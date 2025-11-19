load('LTADataCell_NoFiltfilt.mat')
load('procBallDataCell.mat')
dataCellOld = dataCellSaveNew;
for n = 1:size(dataCellOld,1)
    disp(num2str(n))
    if isempty(dataCellOld{n,2})
        continue
    end
    fileName = dataCellOld{n,1};
    fileName = ['D' fileName(2:end)];
    load(fileName)
%     ballData = load([fileName(1:31) '.txt']);
%     procData = smoothBallData(ballData(:,1:2),10000);
%     procData(:,2) = convertBallVoltToMPS(abs(procData(:,2)));
    procData = procBallDataCell{n,2};
    startStopInd = dataCellOld{n,2};
    [pks,loc] = findpeaks(procData(:,2),procData(:,1));
    close all
    brainMotionX = [];
    brainMotionY = [];
    newTimes = [];
    triggerShift = [];
    for i = 1:size(dataCellOld{n,2},1)
        earlyVals = pks(loc < startStopInd(i,2));
        earlyVal = earlyVals(end);
        lateVals = pks(loc > startStopInd(i,2));
        lateVal = lateVals(1);
        if lateVal >= earlyVal
            newTriggerTime = loc(length(earlyVals)+1);
        else
            newTriggerTime = loc(length(earlyVals)+1);
        end
        triggerShift(end+1) = newTriggerTime - startStopInd(i,2);
%         close all
%         startTime = startStopInd(i,1);
%         stopTime = startStopInd(i,3);
%         plot(movementData.ballData(:,1),movementData.ballData(:,2),'r')
%         hold on
%         plot(procData(:,1),procData(:,2),'b')
%         plot([startStopInd(i,2),startStopInd(i,2)],[-1 1],'k--')
%         xlim([startTime stopTime])
%         ylim([0 0.25])
%         [~,locThresh] = ginput(1);
%         newTriggerTime = procData(find(procData(:,2)>locThresh & procData(:,1)>startTime & procData(:,1)<stopTime,1),1);
%         plot([newTriggerTime,newTriggerTime],[-1 1],'g--')
%         close all
        newTimesRun = [newTriggerTime-2 newTriggerTime newTriggerTime+3];
        newTimesRun(4:6) = round(newTimesRun./movementData.secondsPerFrame);
        brainMotionXRun = movementData.targetPosition(:,1);
        brainMotionYRun = movementData.targetPosition(:,2);
        timeBrain = movementData.secondsPerFrame*[1:size(movementData.targetPosition,1)];
        [minVal,closestIdx] = min(abs(newTimesRun(2)-timeBrain));
        brainMotionX(end+1,:) = medfilt1(brainMotionXRun(closestIdx-79:closestIdx+118)-brainMotionXRun(closestIdx-79),6);
        brainMotionY(end+1,:) = medfilt1(brainMotionYRun(closestIdx-79:closestIdx+118)-brainMotionYRun(closestIdx-79),6);
        newTimes(end+1,:) = newTimesRun;
%         brainMotionX = brainMotionX(newTimes(1) < timeBrain & timeBrain < newTimes(3));
%         brainMotionY = brainMotionY(newTimes(1) < timeBrain & timeBrain < newTimes(3));

        
    end
    if size(dataCellOld{n,3},2) > 198
        dataCellSaveNew{n,3} = dataCellOld{n,3}(:,1:end-1);
        dataCellSaveNew{n,4} = dataCellOld{n,4}(:,1:end-1);
    end
    dataCellSaveNew{n,2} = newTimes;
    if movementData.hemisphere == 2
        brainMotionX = brainMotionX * -1;
    end
    dataCellSaveNew{n,3}(2,:) = mean(brainMotionX);
    dataCellSaveNew{n,4}(2,:) = mean(brainMotionY);
    
%     close all
%     plot(mean(brainMotionY),'b')
%     hold on
%     plot(dataCellOld{n,4}(2,:),'r')
%     hold off
%     startStopInd = dataCellOld{n,5};
%     [pks,loc] = findpeaks(procData(:,2),procData(:,1));
%     close all
    startStopInd = dataCellOld{n,5};
    brainMotionX = [];
    brainMotionY = [];
    newTimes = [];
    for i = 1:size(dataCellOld{n,5},1)
        earlyVals = pks(loc < startStopInd(i,2));
        earlyVal = earlyVals(end);
        lateVals = pks(loc > startStopInd(i,2));
        lateVal = lateVals(1);
        if lateVal <= earlyVal
            newTriggerTime = loc(length(earlyVals)+1);
        else
            newTriggerTime = loc(length(earlyVals));
        end
%         close all
%         startTime = startStopInd(i,1);
%         stopTime = startStopInd(i,3);
%         plot(movementData.ballData(:,1),movementData.ballData(:,2),'r')
%         hold on
%         plot(procData(:,1),procData(:,2),'b')
%         plot([startStopInd(i,2),startStopInd(i,2)],[-1 1],'k--')
%         xlim([startTime stopTime])
%         ylim([0 0.25])
%         [~,locThresh] = ginput(1);
%         newTriggerTime = procData(find(procData(:,2)>locThresh & procData(:,1)>startTime & procData(:,1)<stopTime,1,'last'),1);
%         plot([newTriggerTime,newTriggerTime],[-1 1],'g--')
%         close all
        newTimesRun = [newTriggerTime-2 newTriggerTime newTriggerTime+3];
        newTimesRun(4:6) = round(newTimesRun./movementData.secondsPerFrame);
        brainMotionXRun = movementData.targetPosition(:,1);
        brainMotionYRun = movementData.targetPosition(:,2);
        timeBrain = movementData.secondsPerFrame*[1:size(movementData.targetPosition,1)];
        [minVal,closestIdx] = min(abs(newTimesRun(2)-timeBrain));
        brainMotionX(end+1,:) = medfilt1(brainMotionXRun(closestIdx-79:closestIdx+118)-brainMotionXRun(closestIdx-79),6);
        brainMotionY(end+1,:) = medfilt1(brainMotionYRun(closestIdx-79:closestIdx+118)-brainMotionYRun(closestIdx-79),6);
        newTimes(end+1,:) = newTimesRun;
%         brainMotionX = brainMotionX(newTimes(1) < timeBrain & timeBrain < newTimes(3));
%         brainMotionY = brainMotionY(newTimes(1) < timeBrain & timeBrain < newTimes(3));

        
    end
    if size(dataCellOld{n,6},2) > 198
        dataCellSaveNew{n,6} = dataCellOld{n,6}(:,1:end-1);
        dataCellSaveNew{n,7} = dataCellOld{n,7}(:,1:end-1);
    end
    dataCellSaveNew{n,5} = newTimes;
    if movementData.hemisphere == 2
        brainMotionX = brainMotionX * -1;
    end
    dataCellSaveNew{n,6}(2,:) = mean(brainMotionX);
    dataCellSaveNew{n,7}(2,:) = mean(brainMotionY);
end
    mean(triggerShift)
% load('dataCellSaveNew.mat')
% for n = 1:size(dataCellSaveNew,1)
%     if isnan(dataCellSaveNew{n,3})
%         continue
%     end
%     dataCellSaveNew{n,3}(2,:) = dataCellSaveNew{n,3}(2,:) - dataCellSaveNew{n,3}(2,1);
%     dataCellSaveNew{n,4}(2,:) = (dataCellSaveNew{n,4}(2,:) - dataCellSaveNew{n,4}(2,1))*-1;
%     dataCellSaveNew{n,6}(2,:) = dataCellSaveNew{n,6}(2,:) - dataCellSaveNew{n,6}(2,1);
%     dataCellSaveNew{n,7}(2,:) = (dataCellSaveNew{n,7}(2,:) - dataCellSaveNew{n,7}(2,1))*-1;
% end

motionEventsLocationsX = [];
motionEventsLocationsY = [];
timeToThreshX = [];
timeToThreshY = [];
thresh = 1;
% load('LTADataCell_FS.mat')
% load('LTAFiltFIlt.mat')
locDataCell = dataCellSaveNew;
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
    
    singleTimeVecX = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(motionVectorX));
    idxToThreshXSingle = find((motionVectorX - motionVectorX(find(singleTimeVecX>0,1))) > thresh & singleTimeVecX>0,1);
%     plot(singleTimeVecX,motionVectorX - motionVectorX(find(singleTimeVecX>0,1)))
%     hold on
%     plot(singleTimeVecX(idxToThreshXSingle),0,'x')
%     plot(0,0,'x')
%     hold off
    if ~isempty(idxToThreshXSingle)
        timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
    end
    singleTimeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(motionVectorY));
    idxToThreshYSingle = find((motionVectorY - motionVectorY(find(singleTimeVecY>0,1))) > thresh & singleTimeVecY>0,1);
    if ~isempty(idxToThreshYSingle)
        timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
    end
end
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts(motionEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts(motionEventsLocationsY,90);
% meanY = -1*meanY;
% cIntFillPtsY = -1*cIntFillPtsY;
timeVecX = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(meanX));
timeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(meanY));

h(1) = figure('Color','White');
subplot(2,1,1)
histfit(timeToThreshX,20,'kernel')
title('Time for brain to displace laterally 0.5 micrometers following locomotion trigger')
xlabel('Time (s)')
xlim([-2 3])
text(-1.5,5,['thresh = ' num2str(thresh) ', n = ' num2str(length(timeToThreshX))])

subplot(2,1,2)
histfit(timeToThreshY,20,'kernel')
title('Time for brain to displace rostrally 0.5 micrometers following locomotion trigger')
xlabel('Time (s)')
xlim([-2 3])
text(-1.5,5,['thresh = ' num2str(thresh) ', n = ' num2str(length(timeToThreshY))])

h(2) = figure('Color','White');
subplot(2,2,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
plot(timeVecX,meanX,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
set(f,'facea',[.2]);
plot([0 0],[-3 3],'r')
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[0,0,1,0.1])
end
hold off
text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean Motion During Locomotion Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Brian Shift (\mum)')
ylim([-3 3])
xlim([-2 3])
grid on

subplot(2,2,3)
plot(timeVecY,meanY,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
set(f,'facea',[.2]);
plot([0 0],[-3 3],'r')
for n = 1:size(motionEventsLocationsY,1)
    plot(timeVecY,motionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
end
hold off
text(3,-3,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Rostral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
xlabel('Time (s)')
ylabel('\Delta Brian Shift (\mum)')
ylim([-3 3])
xlim([-2 3])
grid on
clear movementData

stopMotionEventsLocationsX = [];
stopMotionEventsLocationsY = [];
% load('LTADataCell_FS.mat')
% load('LTAFiltFIlt.mat')
locDataCell = dataCellSaveNew;
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
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts(stopMotionEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts(stopMotionEventsLocationsY,90);
% meanY = -1*meanY;
% cIntFillPtsY = -1*cIntFillPtsY;
timeVecX = linspace(round(locDataCell{1,5}(1,1)-locDataCell{1,5}(1,2)),round(locDataCell{1,5}(1,3)-locDataCell{1,5}(1,2)),length(meanX));
timeVecY = linspace(round(locDataCell{1,5}(1,1)-locDataCell{1,5}(1,2)),round(locDataCell{1,5}(1,3)-locDataCell{1,5}(1,2)),length(meanY));

subplot(2,2,2)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
plot(timeVecX,meanX,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
set(f,'facea',[.2]);
plot([0 0],[-3 3],'r')
for n = 1:size(stopMotionEventsLocationsX,1)
    plot(timeVecX,stopMotionEventsLocationsX(n,:),'Color',[0,0,1,0.1])
end
hold off
text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
title(['\fontsize{20pt}\bf{Mean Motion During Stopping Locomotion Events, n = ' num2str(size(stopMotionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Brian Shift (\mum)')
ylim([-3 3])
xlim([-2 3])
grid on

subplot(2,2,4)
plot(timeVecY,meanY,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
set(f,'facea',[.2]);
plot([0 0],[-3 3],'r')
for n = 1:size(stopMotionEventsLocationsY,1)
    plot(timeVecY,stopMotionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
end
hold off
text(3,-3,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Rostral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
xlabel('Time (s)')
ylabel('\Delta Brian Shift (\mum)')
ylim([-3 3])
xlim([-2 3])
grid on
clear movementData

motionEventsLocationsX = [];
motionEventsLocationsY = [];
timeToThreshX = [];
timeToThreshY = [];
% load('LTADataCell_FS.mat')
load('ETADataCell.mat')
locDataCell = EMGDataCell;
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
    
%     idxToThreshXSingle = find(motionVectorX > thresh,1);
%     if ~isempty(idxToThreshXSingle)
%         idxToThreshX(end+1) = idxToThreshXSingle;
%     end
%     idxToThreshYSingle = find(motionVectorY > thresh,1);
%     if ~isempty(idxToThreshYSingle)
%         idxToThreshY(end+1) = idxToThreshYSingle;
%     end
    
    singleTimeVecX = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(motionVectorX));
    idxToThreshXSingle = find((motionVectorX - motionVectorX(find(singleTimeVecX>0,1))) > thresh & singleTimeVecX>0,1);
    if ~isempty(idxToThreshXSingle)
        timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
    end
    singleTimeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(motionVectorY));
    idxToThreshYSingle = find((motionVectorY - motionVectorY(find(singleTimeVecY>0,1)))*-1 > thresh & singleTimeVecY>0,1);
    if ~isempty(idxToThreshYSingle)
        timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
    end
end
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts(motionEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts(motionEventsLocationsY,90);
% meanY = -1*meanY;
% cIntFillPtsY = -1*cIntFillPtsY;
timeVecX = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(meanX));
timeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(meanY));

h(3) = figure('Color','White');
subplot(2,1,1)
histfit(timeToThreshX,30,'kernel')
title('Time for brain to displace laterally 0.5 micrometers following EMG trigger')
xlabel('Time (s)')
xlim([-2 3])
text(-1.5,5,['thresh = ' num2str(thresh) ', n = ' num2str(length(timeToThreshX))])

subplot(2,1,2)
histfit(timeToThreshY,30,'kernel')
title('Time for brain to displace rostrally 0.5 micrometers following EMG trigger')
xlabel('Time (s)')
xlim([-2 3])
text(-1.5,5,['thresh = ' num2str(thresh) ', n = ' num2str(length(timeToThreshY))])