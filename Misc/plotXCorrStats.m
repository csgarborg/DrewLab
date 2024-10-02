dateCell = {'211021','211022','211101','211102','211105','211109','211112','211116','211117','211119','211203','211216','220203','220209','220210','220211','220214','220221','220223','220303','220308','220309','220314','220318','220404','220406','220407','220429','220509','220511','220712','220714','220718','220719','220808','220809','220813','220815','220816','220822','220823','221205','221207','221208','221213'};
close all
fileName = {};
for n = 1:size(dateCell,2)
    folderName = ['D:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/'];
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
        if exist(['D:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processe_2layerBrainInSkullDataFinal.mat'],'file')
            fileName{end+1} = ['D:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processe_2layerBrainInSkullDataFinal.mat'];
        end
    end
end
medFiltSize = 6;

xCorrXLoc = [];
xCorrYLoc = [];
xCorrXEMG = [];
xCorrYEMG = [];
maxlag = 500;
for i = 1:numel(fileName)
%     load([fileName{i}(1:31) '.txt'])
    load(fileName{i});
    motionVectorX = movementData.targetPosition(:,1);
    motionVectorY = movementData.targetPosition(:,2);
    if movementData.hemisphere == 2
        movementData.targetPosition(:,1) = movementData.targetPosition(:,1) * -1;
    end
%     movementData.secondsPerFrame = movementData.secondsPerFrame/2;
    movement_time=movementData.secondsPerFrame*(1:length(motionVectorX));
    locDataInterp=zeros(size(movementData.targetPosition));
    locDataInterp(:,1)=movement_time;
    locDataInterp(:,2)=interp1(movementData.ballData(:,1),abs(movementData.ballData(:,2)),movement_time,'linear');
    emgDataInterp=zeros(size(movementData.targetPosition));
    emgDataInterp(:,1)=movement_time;
    emgDataInterp(:,2)=interp1(movementData.emgData(:,1),abs(movementData.emgData(:,2)),movement_time,'linear');
    xc_1=xcorr(detrend(locDataInterp(1:(end-100),2)-mean(locDataInterp(1:(end-100),2)))', detrend(movementData.targetPosition(1:(end-100),1)-mean(movementData.targetPosition(1:(end-100),1)))',maxlag,'coeff');
    xc_2=xcorr(detrend(locDataInterp(1:(end-100),2)-mean(locDataInterp(1:(end-100),2)))', detrend(movementData.targetPosition(1:(end-100),2)-mean(movementData.targetPosition(1:(end-100),2)))',maxlag,'coeff');
    xc_3=xcorr(detrend(emgDataInterp(1:(end-100),2)-mean(emgDataInterp(1:(end-100),2)))', detrend(movementData.targetPosition(1:(end-100),1)-mean(movementData.targetPosition(1:(end-100),1)))',maxlag,'coeff');
    xc_4=xcorr(detrend(emgDataInterp(1:(end-100),2)-mean(emgDataInterp(1:(end-100),2)))', detrend(movementData.targetPosition(1:(end-100),2)-mean(movementData.targetPosition(1:(end-100),2)))',maxlag,'coeff');
    if i > 1
        if length(xc_1) > size(xCorrXLoc,2)
            xc_1 = xc_1(1:size(xCorrXLoc,2));
        elseif length(xc_1) < size(xCorrXLoc,2)
            xCorrXLoc = xCorrXLoc(:,1:length(xc_1));
        end
        if length(xc_2) > size(xCorrYLoc,2)
            xc_2 = xc_2(1:size(xCorrYLoc,2));
        elseif length(xc_2) < size(xCorrYLoc,2)
            xCorrYLoc = xCorrYLoc(:,1:length(xc_2));
        end
        if length(xc_3) > size(xCorrXEMG,2)
            xc_3 = xc_3(1:size(xCorrXEMG,2));
        elseif length(xc_3) < size(xCorrXEMG,2)
            xCorrXEMG = xCorrXEMG(:,1:length(xc_3));
        end
        if length(xc_4) > size(xCorrYEMG,2)
            xc_4 = xc_4(1:size(xCorrYEMG,2));
        elseif length(xc_4) < size(xCorrYEMG,2)
            xCorrYEMG = xCorrYEMG(:,1:length(xc_4));
        end
    end
%     if min(xc_1) > -.1 && max(xc_1) < .1
%         disp(fileName{i})
%         timeVecX = movementData.secondsPerFrame*(-maxlag:maxlag);
%         plot(timeVecX,xc_1)
%         figure
%         plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)))
%         hold on
%         plot(locDataInterp(:,1),locDataInterp(:,2))
%         figure
%         plot(movementData.emgData(:,1),abs(movementData.emgData(:,2)))
%         hold on
%         plot(emgDataInterp(:,1),emgDataInterp(:,2))
%         disp('check')
%         close
%     end
    if ~any(isnan(xc_1))
        xCorrXLoc(end+1,:) = xc_1;
    end
    if ~any(isnan(xc_2))
        xCorrYLoc(end+1,:) = xc_2;
    end
    if ~any(isnan(xc_3))
        xCorrXEMG(end+1,:) = xc_3;
    end
    if ~any(isnan(xc_4))
        xCorrYEMG(end+1,:) = xc_4;
    end
    clear movementData
end
load(fileName{1});
[meanXLoc,cIntFillPtsXLoc] = getCIntMeanAndFillPts(xCorrXLoc,90);
[meanYLoc,cIntFillPtsYLoc] = getCIntMeanAndFillPts(xCorrYLoc,90);
[meanXEMG,cIntFillPtsXEMG] = getCIntMeanAndFillPts(xCorrXEMG,90);
[meanYEMG,cIntFillPtsYEMG] = getCIntMeanAndFillPts(xCorrYEMG,90);
timeVecX = movementData.secondsPerFrame*(-maxlag:maxlag);
timeVecY = movementData.secondsPerFrame*(-maxlag:maxlag);


% targetPositionInSkull = combineBrainSkullMovement_FS(movementData,stationaryData);
% movementData.secondsPerFrame = movementData.secondsPerFrame/2;
% movementData.targetPosition = targetPositionInSkull;
% movement_time=movementData.secondsPerFrame*(1:length(movementData.targetPosition));
% movementData.locDataInterp=zeros(size(movementData.targetPosition));
% movementData.locDataInterp(:,1)=movement_time;
% movementData.locDataInterp(:,2)=interp1(movementData.ballData(:,1),abs(movementData.ballData(:,2)),movement_time,'linear');

h(1) = figure('Color','White');
subplot(2,2,1)
maxMeanVal = max(abs([meanXLoc meanYLoc cIntFillPtsXLoc cIntFillPtsYLoc]));
plot(timeVecX,meanXLoc,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsXLoc,'r','Linestyle','none');
set(f,'facea',[.2]);
for n = 1:size(xCorrXLoc,1)
    plot(timeVecX,xCorrXLoc(n,:),'Color',[0,0,1,0.1])
end
peak = max(meanXLoc);
peakTime = timeVecX(meanXLoc == peak);
plot(peakTime,peak,'kx','MarkerSize',25)
peakVals = [];
for n = 1:size(xCorrXLoc,1)
    peakVals(end+1) = xCorrXLoc(n,timeVecX == peakTime);
end
stdPeak = std(peakVals);
text(5,.8,['x=' num2str(peakTime) ', y=' num2str(peak) ', std=+-' num2str(stdPeak)])
% plot(timeVecX(brainMotionStart),meanXLoc(brainMotionStart),'rx')
hold off
title(['Figure 2d' 10 'Brain Motion and Locomotion Cross-Correlation - X'])
ylabel('Noramlized Cross-Correlation')
xlabel('Lags (s)')
ylim([0 1])


subplot(2,2,2)
maxMeanVal = max(abs([meanXLoc meanYLoc cIntFillPtsXLoc cIntFillPtsYLoc]));
plot(timeVecY,meanYLoc,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsYLoc,'r','Linestyle','none');
set(f,'facea',[.2]);
for n = 1:size(xCorrYLoc,1)
    plot(timeVecY,xCorrYLoc(n,:),'Color',[0,0,1,0.1])
end
peak = max(meanYLoc);
peakTime = timeVecY(meanYLoc == peak);
plot(peakTime,peak,'kx','MarkerSize',25)
peakVals = [];
for n = 1:size(xCorrYLoc,1)
    peakVals(end+1) = xCorrYLoc(n,timeVecY == peakTime);
end
stdPeak = std(peakVals);
text(5,.8,['x=' num2str(peakTime) ', y=' num2str(peak) ', std=+-' num2str(stdPeak)])
% plot(timeVecX(brainMotionStart),meanXLoc(brainMotionStart),'rx')
hold off
title(['Figure 2d' 10 'Brain Motion and Locomotion Cross-Correlation - Y'])
ylabel('Noramlized Cross-Correlation')
xlabel('Lags (s)')
ylim([0 1])


subplot(2,2,3)
maxMeanVal = max(abs([meanXEMG meanYEMG cIntFillPtsXEMG cIntFillPtsYEMG]));
plot(timeVecX,meanXEMG,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsXEMG,'r','Linestyle','none');
set(f,'facea',[.2]);
for n = 1:size(xCorrXEMG,1)
    plot(timeVecX,xCorrXEMG(n,:),'Color',[0,0,1,0.1])
end
peak = max(meanXEMG);
peakTime = timeVecX(meanXEMG == peak);
plot(peakTime,peak,'kx','MarkerSize',25)
peakVals = [];
for n = 1:size(xCorrXEMG,1)
    peakVals(end+1) = xCorrXEMG(n,timeVecX == peakTime);
end
stdPeak = std(peakVals);
text(5,.8,['x=' num2str(peakTime) ', y=' num2str(peak) ', std=+-' num2str(stdPeak)])
% plot(timeVecX(brainMotionStart),meanXEMG(brainMotionStart),'rx')
hold off
title(['Figure 2d' 10 'Brain Motion and EMG Cross-Correlation - X'])
ylabel('Noramlized Cross-Correlation')
xlabel('Lags (s)')
ylim([0 1])

subplot(2,2,4)
maxMeanVal = max(abs([meanYEMG meanYEMG cIntFillPtsYEMG cIntFillPtsYEMG]));
plot(timeVecY,meanYEMG,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsYEMG,'r','Linestyle','none');
set(f,'facea',[.2]);
for n = 1:size(xCorrYEMG,1)
    plot(timeVecY,xCorrYEMG(n,:),'Color',[0,0,1,0.1])
end
peak = max(meanYEMG);
peakTime = timeVecY(meanYEMG == peak);
plot(peakTime,peak,'kx','MarkerSize',25)
peakVals = [];
for n = 1:size(xCorrYEMG,1)
    peakVals(end+1) = xCorrYEMG(n,timeVecY == peakTime);
end
stdPeak = std(peakVals);
text(5,.8,['x=' num2str(peakTime) ', y=' num2str(peak) ', std=+-' num2str(stdPeak)])
% plot(timeVecY(brainMotionStart),meanYEMG(brainMotionStart),'rx')
hold off
title(['Figure 2d' 10 'Brain Motion and EMG Cross-Correlation - Y'])
ylabel('Noramlized Cross-Correlation')
xlabel('Lags (s)')
ylim([0 1])

% h(5) = figure('Color','White');
% maxlag=500;
% xc_1=xcorr(detrend(movementData.locDataInterp(1:(end-100),2))', detrend(movementData.targetPosition(1:(end-100),1))',maxlag,'coeff');
% xc_2=xcorr(detrend(movementData.locDataInterp(1:(end-100),2))', detrend(movementData.targetPosition(1:(end-100),2))',maxlag,'coeff');
% plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_1,'b')
% hold on
% plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_2,'g')
% hold off
% title(['Figure 2d' 10 'Brain Motion and Locomotion Cross-Correlation'])
% ylabel('Noramlized Cross-Correlation')
% xlabel('Lags (s)')
% legend({'L/M Brain Motion','R/C Brain Motion'})
% ylim([0 1])
% axes('Position',[.2 .7 .2 .2])
% box on
% plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_1,'b')
% hold on
% plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_2,'g')
% hold off
% axis([-2 2 .6 .7])