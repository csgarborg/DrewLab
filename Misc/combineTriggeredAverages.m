close all
fileName = {'C:\Workspace\210218_005_processed_combined.mat','C:\Workspace\210218_006_processed_combinedMean.mat','C:\Workspace\210218_007_processed_combined.mat','C:\Workspace\210218_008_processed_combined.mat'};
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
plot([0 0],[-maxMeanVal maxMeanVal],'r')
hold off
title(['\fontsize{20pt}\bf{Mean Motion During Locomotion Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('X Position (\mum)')
ylim([-maxMeanVal maxMeanVal])
grid on
subplot(2,1,2)
plot(timeVecY,meanY,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
set(f,'facea',[.2]);
plot([0 0],[-maxMeanVal maxMeanVal],'r')
hold off
xlabel('Time (s)')
ylabel('Y Position (\mum)')
ylim([-maxMeanVal maxMeanVal])
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
plot([0 0],[-maxMeanVal maxMeanVal],'r')
hold off
title(['\fontsize{20pt}\bf{Mean Motion During EMG Events, n = ' num2str(size(EMGEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('X Position (\mum)')
ylim([-maxMeanVal maxMeanVal])
grid on
subplot(2,1,2)
plot(timeVecY,meanY,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
set(f,'facea',[.2]);
plot([0 0],[-maxMeanVal maxMeanVal],'r')
hold off
xlabel('Time (s)')
ylabel('Y Position (\mum)')
ylim([-maxMeanVal maxMeanVal])
grid on
clear movementData