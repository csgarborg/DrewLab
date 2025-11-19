% load('D:\24-05-13_MouseExp\240513_005_processed_Layer1_combined.mat')
% combinedMovementDataBrain_3f = movementData;
clear
close all
% load('squeezeDataCell_FS.mat')
load('ETADataCell_FS.mat')
if exist('squeezeDataCell','var')
    locDataCell = squeezeDataCell;
elseif exist('EMGDataCell','var')
    locDataCell = EMGDataCell;
end
for i = 1:size(locDataCell,1)
load(locDataCell{i,1})
%     locIdx = find(strcmp(squeezeDataCell(:,1),'D:/21-12-16_MouseExp/211216_002_processe_2layerBrainInSkullDataFinal.mat'));
locIdx = i;
%     targetPositionInSkull = combineBrainSkullMovement_FS(movementData,stationaryData);
% movementData.secondsPerFrame = movementData.secondsPerFrame/2;
%     posL2 = [movementData.emgData(:,1), movementData.emgData(:,2)];
%     meanPosValL2 = [.5*(posL2(1:end-1,1) + posL2(2:end,1)),.5*(posL2(1:end-1,2) + posL2(2:end,2))];
%     posL2 = zipperVecs(posL2,meanPosValL2);
    emgEventsLocationsXRegTime = [];
%     motionEventsLocationsY = [];
%     deleteRows = [];
% movementData.secondsPerFrame = 0.0252*2;
locTriggerEMGIdxReg = [];
    for n = 1:size(locDataCell{locIdx,2},1)
        locTriggerTime = locDataCell{locIdx,2}(n,2);
%         locEMGOverlapCheck = find(abs(locTriggerTime - EMGDataCell{emgIdx,2}(:,2)) < timeOverlapThresh);
%         if isempty(locEMGOverlapCheck)
%             badOn = badOn+1;
%             deleteRows = [deleteRows n];
%             continue
%         elseif length(locEMGOverlapCheck) > 1
%             disp('wtf')
%         else
%             goodOn = goodOn+1;
%         end
        timeVec = (1:length(movementData.targetPosition(:,1)))*movementData.secondsPerFrame;
        locTriggerEMGIdxReg(n) = find(min(abs(locTriggerTime - timeVec)) == abs(locTriggerTime - timeVec));
%         locDataCell{i,2}(n,4:6) = [locTriggerEMGIdx-60 locTriggerEMGIdx locTriggerEMGIdx+90];
        emgVector = movementData.targetPosition(locTriggerEMGIdxReg(n)-round(2/movementData.secondsPerFrame):locTriggerEMGIdxReg(n)+round(3/movementData.secondsPerFrame),2);
%         if movementData.hemisphere == 2
%             emgVector = emgVector*-1;
%         end
%         motionVectorY = movementData.emgData(locDataCell{i,2}(n,4):locDataCell{i,2}(n,6),2);
%         if n > 1
%             if length(emgVector) > size(motionEventsLocationsX,2)
%                 emgVector = emgVector(1:size(motionEventsLocationsX,2));
%             elseif length(emgVector) < size(motionEventsLocationsX,2)
%                 motionEventsLocationsX = motionEventsLocationsX(:,1:length(emgVector));
%             end
%             if length(motionVectorY) > size(motionEventsLocationsY,2)
%                 motionVectorY = motionVectorY(1:size(motionEventsLocationsY,2));
%             elseif length(motionVectorY) < size(motionEventsLocationsY,2)
%                 motionEventsLocationsY = motionEventsLocationsY(:,1:length(motionVectorY));
%             end
%         end
        emgEventsLocationsXRegTime(end+1,:) = emgVector - emgVector(1);
%         motionEventsLocationsY(end+1,:) = motionVectorY - motionVectorY(1);
    end
    emgEventsLocationsXAltTime = [];
movementData.secondsPerFrame = 0.0252;
locTriggerEMGIdxAlt = [];
    for n = 1:size(locDataCell{locIdx,2},1)
        locTriggerTime = locDataCell{locIdx,2}(n,2);
%         locEMGOverlapCheck = find(abs(locTriggerTime - EMGDataCell{emgIdx,2}(:,2)) < timeOverlapThresh);
%         if isempty(locEMGOverlapCheck)
%             badOn = badOn+1;
%             deleteRows = [deleteRows n];
%             continue
%         elseif length(locEMGOverlapCheck) > 1
%             disp('wtf')
%         else
%             goodOn = goodOn+1;
%         end
        timeVec = (1:length(movementData.targetPosition(:,1)))*movementData.secondsPerFrame;
        locTriggerEMGIdxAlt(n) = find(min(abs(locTriggerTime - timeVec)) == abs(locTriggerTime - timeVec));
%         locDataCell{i,2}(n,4:6) = [locTriggerEMGIdx-60 locTriggerEMGIdx locTriggerEMGIdx+90];
        emgVector = movementData.targetPosition(locTriggerEMGIdxAlt(n)-round(2/movementData.secondsPerFrame):locTriggerEMGIdxAlt(n)+round(3/movementData.secondsPerFrame),2);
%         if movementData.hemisphere == 2
%             emgVector = emgVector*-1;
%         end
%         motionVectorY = movementData.emgData(locDataCell{i,2}(n,4):locDataCell{i,2}(n,6),2);
%         if n > 1
%             if length(emgVector) > size(motionEventsLocationsX,2)
%                 emgVector = emgVector(1:size(motionEventsLocationsX,2));
%             elseif length(emgVector) < size(motionEventsLocationsX,2)
%                 motionEventsLocationsX = motionEventsLocationsX(:,1:length(emgVector));
%             end
%             if length(motionVectorY) > size(motionEventsLocationsY,2)
%                 motionVectorY = motionVectorY(1:size(motionEventsLocationsY,2));
%             elseif length(motionVectorY) < size(motionEventsLocationsY,2)
%                 motionEventsLocationsY = motionEventsLocationsY(:,1:length(motionVectorY));
%             end
%         end
        emgEventsLocationsXAltTime(end+1,:) = emgVector - emgVector(1);
%         motionEventsLocationsY(end+1,:) = motionVectorY - motionVectorY(1);
    end
%     if size(emgEventsLocationsX,2) == 199
%         emgEventsLocationsX = emgEventsLocationsX(:,1:end-1);
%         motionEventsLocationsY = motionEventsLocationsY(:,1:end-1);
%     end
% if isempty(emgEventsLocationsX)
%     locDataCell{i,3} = NaN;
% else
%     locDataCell{i,3} = [linspace(-2,3,size(emgEventsLocationsX,2)); mean(emgEventsLocationsX,1)];
% end
%     locDataCell{i,2}(deleteRows,:) = [];
    
%     locDataCell{i,4}(2,:) = mean(motionEventsLocationsY);
    
%     emgEventsLocationsX = [];
%     motionEventsLocationsY = [];
%     deleteRows = [];
%     for n = 1:size(locDataCell{i,5},1)
%         locTriggerTime = locDataCell{i,5}(n,2);
%         locEMGOverlapCheck = find(abs(locTriggerTime - EMGDataCell{emgIdx,2}(:,2)) < timeOverlapThresh);
%         if isempty(locEMGOverlapCheck)
%             badOff = badOff+1;
%             deleteRows = [deleteRows n];
%             continue
%         elseif length(locEMGOverlapCheck) > 1
%             disp('wtf')
%         else
%             goodOff = goodOff+1;
%         end
%         locTriggerEMGIdx = find(min(abs(locTriggerTime - movementData.emgData(:,1))) == abs(locTriggerTime - movementData.emgData(:,1)));
%         locDataCell{i,5}(n,4:6) = [locTriggerEMGIdx-60 locTriggerEMGIdx locTriggerEMGIdx+90];
%         emgVector = movementData.emgData(locTriggerEMGIdx-60:locTriggerEMGIdx+90,2);
%         if n > 1
%             if length(emgVector) > size(emgEventsLocationsX,2)
%                 emgVector = emgVector(1:size(emgEventsLocationsX,2));
%             elseif length(emgVector) < size(emgEventsLocationsX,2)
%                 emgEventsLocationsX = emgEventsLocationsX(:,1:length(emgVector));
%             end
%             if length(motionVectorY) > size(motionEventsLocationsY,2)
%                 motionVectorY = motionVectorY(1:size(motionEventsLocationsY,2));
%             elseif length(motionVectorY) < size(motionEventsLocationsY,2)
%                 motionEventsLocationsY = motionEventsLocationsY(:,1:length(motionVectorY));
%             end
%         end
%         emgEventsLocationsX(end+1,:) = emgVector - emgVector(1) + 1;
%         motionEventsLocationsY(end+1,:) = motionVectorY - motionVectorY(1);
%     end
%     if size(emgEventsLocationsX,2) == 199
%         emgEventsLocationsX = emgEventsLocationsX(:,1:end-1);
% %         motionEventsLocationsY = motionEventsLocationsY(:,1:end-1);
%     end
% if isempty(emgEventsLocationsX)
%     locDataCell{i,6} = NaN;
% else
%     locDataCell{i,6} = [linspace(-2,3,size(emgEventsLocationsX,2)); mean(emgEventsLocationsX,1)];
% end
%     locDataCell{i,5}(deleteRows,:) = [];
%     locDataCell{i,7}(2,:) = mean(motionEventsLocationsY);
% disp('goodOn badOn goodOff badOff')
% [goodOn badOn goodOff badOff]
% [badLoc badEMG goodLocEMG]
% locDataCell(:,[4 7]) = [];
% save('LTADataCellEMG_FS.mat','locDataCell')

motionEventsLocationsXRegTime = [];
% motionEventsLocationsY = [];
timeToThreshX = [];
% timeToThreshY = [];
dispTimeThreshX = [];
% dispTimeThreshY = [];
moveThresh = 1.5;
timeThresh = -.25;
binEdgesTime = -2:.1:3;
binEdgesDisp = 0:.1:4;
% load('LTADataCellEMG_FS.mat')
for n = 1:size(emgEventsLocationsXRegTime)
%     if isnan(locDataCell{n,3})
%         continue
%     end
    motionVectorX = emgEventsLocationsXRegTime(n,:);
%     motionVectorY = locDataCell{n,4}(2,:);
%     if n > 1
%         if length(motionVectorX) > size(motionEventsLocationsX,2)
%             motionVectorX = motionVectorX(1:size(motionEventsLocationsX,2));
%         elseif length(motionVectorX) < size(motionEventsLocationsX,2)
%             motionEventsLocationsX = motionEventsLocationsX(:,1:length(motionVectorX));
%         end
%         if length(motionVectorY) > size(motionEventsLocationsY,2)
%             motionVectorY = motionVectorY(1:size(motionEventsLocationsY,2));
%         elseif length(motionVectorY) < size(motionEventsLocationsY,2)
%             motionEventsLocationsY = motionEventsLocationsY(:,1:length(motionVectorY));
%         end
%     end
    motionEventsLocationsXRegTime(end+1,:) = motionVectorX;
%     motionEventsLocationsY(end+1,:) = motionVectorY;
    
    singleTimeVecX = linspace(-2,3,length(motionVectorX));
%     dispTimeThreshX(end+1) = motionVectorX((find(singleTimeVecX>timeThresh,1))) - motionVectorX((find(singleTimeVecX>0,1)));
%     dispTimeThreshX(end+1) = motionVectorX((find(singleTimeVecX>timeThresh,1)));
%     idxToThreshXSingle = find((motionVectorX - motionVectorX(find(singleTimeVecX>0,1))) > moveThresh & singleTimeVecX>0 & singleTimeVecX<=3,1);
%     idxToThreshXSingle = find(motionVectorX > moveThresh & singleTimeVecX<=3,1);
%     if ~isempty(idxToThreshXSingle)
%         timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
%     end
%     singleTimeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(motionVectorY));
%     dispTimeThreshY(end+1) =  (motionVectorY((find(singleTimeVecY>timeThresh,1))) - motionVectorY((find(singleTimeVecY>0,1))))*-1;
%     idxToThreshYSingle = find((motionVectorY - motionVectorY(find(singleTimeVecY>0,1)))*-1 > moveThresh & singleTimeVecY>0 & singleTimeVecY<=3,1);
%     if ~isempty(idxToThreshYSingle)
%         timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
%     end
end
if isempty(motionEventsLocationsXRegTime)
    continue
end
[meanXRegTime,cIntFillPtsXRegTime] = getCIntMeanAndFillPts(motionEventsLocationsXRegTime,90);
motionEventsLocationsXAltTime = [];
for n = 1:size(emgEventsLocationsXAltTime)
%     if isnan(locDataCell{n,3})
%         continue
%     end
    motionVectorX = emgEventsLocationsXAltTime(n,:);
%     motionVectorY = locDataCell{n,4}(2,:);
%     if n > 1
%         if length(motionVectorX) > size(motionEventsLocationsX,2)
%             motionVectorX = motionVectorX(1:size(motionEventsLocationsX,2));
%         elseif length(motionVectorX) < size(motionEventsLocationsX,2)
%             motionEventsLocationsX = motionEventsLocationsX(:,1:length(motionVectorX));
%         end
%         if length(motionVectorY) > size(motionEventsLocationsY,2)
%             motionVectorY = motionVectorY(1:size(motionEventsLocationsY,2));
%         elseif length(motionVectorY) < size(motionEventsLocationsY,2)
%             motionEventsLocationsY = motionEventsLocationsY(:,1:length(motionVectorY));
%         end
%     end
    motionEventsLocationsXAltTime(end+1,:) = motionVectorX;
%     motionEventsLocationsY(end+1,:) = motionVectorY;
    
    singleTimeVecX = linspace(-2,3,length(motionVectorX));
%     dispTimeThreshX(end+1) = motionVectorX((find(singleTimeVecX>timeThresh,1))) - motionVectorX((find(singleTimeVecX>0,1)));
%     dispTimeThreshX(end+1) = motionVectorX((find(singleTimeVecX>timeThresh,1)));
%     idxToThreshXSingle = find((motionVectorX - motionVectorX(find(singleTimeVecX>0,1))) > moveThresh & singleTimeVecX>0 & singleTimeVecX<=3,1);
%     idxToThreshXSingle = find(motionVectorX > moveThresh & singleTimeVecX<=3,1);
%     if ~isempty(idxToThreshXSingle)
%         timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
%     end
%     singleTimeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(motionVectorY));
%     dispTimeThreshY(end+1) =  (motionVectorY((find(singleTimeVecY>timeThresh,1))) - motionVectorY((find(singleTimeVecY>0,1))))*-1;
%     idxToThreshYSingle = find((motionVectorY - motionVectorY(find(singleTimeVecY>0,1)))*-1 > moveThresh & singleTimeVecY>0 & singleTimeVecY<=3,1);
%     if ~isempty(idxToThreshYSingle)
%         timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
%     end
end
[meanXAltTime,cIntFillPtsXAltTime] = getCIntMeanAndFillPts(motionEventsLocationsXAltTime,90);
% [meanY,cIntFillPtsY] = getCIntMeanAndFillPts_FS(motionEventsLocationsY,90);
% meanY = -1*meanY;
% cIntFillPtsY = -1*cIntFillPtsY;
% timeVecX = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(meanX));
timeVecX = linspace(-2,3,length(meanXAltTime));
% timeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(meanY));

% stdWindowSize = 5;
% for n = stdWindowSize+1:length(meanY)
%     if std(meanY(n-stdWindowSize:n)) > .01
%         brainMotionStart = n;
%         break
%     end 
% end

idxDiff = locTriggerEMGIdxAlt - locTriggerEMGIdxReg;
h(9) = figure('Color','White');
subplot(1,2,1)
% maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
plot(timeVecX,meanXRegTime,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsXRegTime,'r','Linestyle','none');
set(f,'facea',[.2]);
plot([0 0],[0 4],'r')
subplot(1,2,2)
% maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
plot(timeVecX,meanXAltTime,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsXAltTime,'r','Linestyle','none');
set(f,'facea',[.2]);
plot([0 0],[0 4],'r')
for n = 1:size(motionEventsLocationsXRegTime,1)
    subplot(1,2,1)
    plot(timeVecX,motionEventsLocationsXRegTime(n,:),'Color',[0,0,1,0.1])
    subplot(1,2,2)
    plot(timeVecX,motionEventsLocationsXAltTime(n,:),'Color',[0,0,1,0.1])
end
subplot(1,2,1)
title('MSCAN SPF - 0.0253')
xlim([-2 3])
text(-2,4,[regexprep(num2str(idxDiff),' +',' ') '  IDX'])
text(-2,3,[regexprep(num2str(round(idxDiff*0.0252,2)),' +',' ') '  SEC'])
subplot(1,2,2)
title('FIT TO ANALOG SPF - 0.0252')
xlim([-2 3])
% plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
% text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean EMG During Locomotion Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('EKG Power (au)')
ylim([0 4])
% xlim([-2 3])
grid on

close all
% subplot(2,2,3)
% plot(timeVecY,meanY,'k')
% hold on
% f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
% set(f,'facea',[.2]);
% plot([0 0],[-3 3],'r')
% for n = 1:size(motionEventsLocationsY,1)
%     plot(timeVecY,-1*motionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
% end
% plot(timeVecY(brainMotionStart),meanY(brainMotionStart),'rx')
% hold off
% text(3,-3,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(3,3,'Rostral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% xlabel('Time (s)')
% ylabel('\Delta Brian Shift (\mum)')
% ylim([-3 3])
% xlim([-2 3])
% grid on
% clear movementData

% stopMotionEventsLocationsX = [];
% % stopMotionEventsLocationsY = [];
% for n = 1:size(locDataCell)
%     if isnan(locDataCell{n,5})
%         continue
%     end
%     stopMotionVectorX = locDataCell{n,5}(2,:);
% %     stopMotionVectorY = locDataCell{n,7}(2,:);
% %     if n > 1
% %         if length(stopMotionVectorX) > size(stopMotionEventsLocationsX,2)
% %             stopMotionVectorX = stopMotionVectorX(1:size(stopMotionEventsLocationsX,2));
% %         elseif length(stopMotionVectorX) < size(stopMotionEventsLocationsX,2)
% %             stopMotionEventsLocationsX = stopMotionEventsLocationsX(:,1:length(stopMotionVectorX));
% %         end
% %         if length(stopMotionVectorY) > size(stopMotionEventsLocationsY,2)
% %             stopMotionVectorY = stopMotionVectorY(1:size(stopMotionEventsLocationsY,2));
% %         elseif length(stopMotionVectorY) < size(stopMotionEventsLocationsY,2)
% %             stopMotionEventsLocationsY = stopMotionEventsLocationsY(:,1:length(stopMotionVectorY));
% %         end
% %     end
%     stopMotionEventsLocationsX(end+1,:) = stopMotionVectorX;
% %     stopMotionEventsLocationsY(end+1,:) = stopMotionVectorY;
% end
% [meanX,cIntFillPtsX] = getCIntMeanAndFillPts_FS(stopMotionEventsLocationsX,90);
% % [meanY,cIntFillPtsY] = getCIntMeanAndFillPts_FS(stopMotionEventsLocationsY,90);
% % meanY = -1*meanY;
% % cIntFillPtsY = -1*cIntFillPtsY;
% timeVecX = linspace(round(locDataCell{2,4}(1,1)-locDataCell{2,4}(1,2)),round(locDataCell{2,4}(1,3)-locDataCell{2,4}(1,2)),length(meanX));
% % timeVecY = linspace(round(locDataCell{1,5}(1,1)-locDataCell{1,5}(1,2)),round(locDataCell{1,5}(1,3)-locDataCell{1,5}(1,2)),length(meanY));
% 
% % subplot(2,1,2)
% % % maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
% % plot(timeVecX,meanX,'k')
% % hold on
% % f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
% % set(f,'facea',[.2]);
% % plot([0 0],[0 3],'r')
% % for n = 1:size(stopMotionEventsLocationsX,1)
% %     plot(timeVecX,stopMotionEventsLocationsX(n,:),'Color',[0,0,1,0.1])
% % end
% % hold off
% % % text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% % % text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% % title(['\fontsize{20pt}\bf{Mean EMG During Stopping Locomotion Events, n = ' num2str(size(stopMotionEventsLocationsX,1)) '}'])
% % xlabel('Time (s)')
% % ylabel('EMG Power (au)')
% % ylim([0 3])
% % xlim([-2 3])
% % grid on
% 
% % subplot(2,2,4)
% % plot(timeVecY,meanY,'k')
% % hold on
% % f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
% % set(f,'facea',[.2]);
% % plot([0 0],[-3 3],'r')
% % for n = 1:size(stopMotionEventsLocationsY,1)
% %     plot(timeVecY,-1*stopMotionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
% % end
% % hold off
% % text(3,-3,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% % text(3,3,'Rostral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% % xlabel('Time (s)')
% % ylabel('\Delta Brian Shift (\mum)')
% % ylim([-3 3])
% % xlim([-2 3])
% % grid on
% % clear movementData
% 
% h(10) = figure('Color','White');
% % subplot(2,1,1)
% histogram(timeToThreshX,binEdgesTime);
% hold on
% [pdfXVals,pdfYVals] = findKernelPDF(timeToThreshX,binEdgesTime);
% plot(pdfXVals,pdfYVals*2,'r','LineWidth',2)
% title('Time for EMG power increase one-half order of magnitude in relation to locomotion trigger')
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 20])
% mu = mean(timeToThreshX);
% sig = std(timeToThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% % subplot(2,1,2)
% % histfit(timeToThreshY,numBins,'kernel')
% % title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following locomotion trigger'])
% % xlabel('Time (s)')
% % xlim([-2 3])
% % ylim([0 25])
% % hold on
% % mu = mean(timeToThreshY);
% % sig = std(timeToThreshY);
% % plot([mu mu],[0 40],'k','LineWidth',2);
% % plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% % plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% % hold off
% % text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% h(50) = figure('Color','White');
% % subplot(2,1,1)
% histogram(dispTimeThreshX,binEdgesDisp);
% hold on
% [pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshX,binEdgesDisp);
% plot(pdfXVals,pdfYVals*2,'r','LineWidth',2)
% title('EMG power 0.25s before locomotion trigger')
% xlabel('Power (au)')
% xlim([0 4])
% ylim([0 20])
% hold on
% mu = mean(dispTimeThreshX);
% sig = std(dispTimeThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% % subplot(2,1,2)
% % histfit(dispTimeThreshY,numBins,'kernel')
% % title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following locomotion trigger'])
% % xlabel('Time (s)')
% % xlim([-3 4])
% % ylim([0 40])
% % hold on
% % mu = mean(dispTimeThreshY);
% % sig = std(dispTimeThreshY);
% % plot([mu mu],[0 40],'k','LineWidth',2);
% % plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% % plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% % hold off
% % text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% % timeToThreshX = [];
% % timeToThreshY = [];
% % dispTimeThreshX = [];
% % dispTimeThreshY = [];
% % for n = 1:size(locDataCell)
% %     if isnan(locDataCell{n,3})
% %         continue
% %     end
% %     motionVectorX = locDataCell{n,3}(2,:);
% %     motionVectorY = locDataCell{n,4}(2,:);
% %     if n > 1
% %         if length(motionVectorX) > size(motionEventsLocationsX,2)
% %             motionVectorX = motionVectorX(1:size(motionEventsLocationsX,2));
% %         elseif length(motionVectorX) < size(motionEventsLocationsX,2)
% %             motionEventsLocationsX = motionEventsLocationsX(:,1:length(motionVectorX));
% %         end
% %         if length(motionVectorY) > size(motionEventsLocationsY,2)
% %             motionVectorY = motionVectorY(1:size(motionEventsLocationsY,2));
% %         elseif length(motionVectorY) < size(motionEventsLocationsY,2)
% %             motionEventsLocationsY = motionEventsLocationsY(:,1:length(motionVectorY));
% %         end
% %     end
% %     motionEventsLocationsX(end+1,:) = motionVectorX;
% %     motionEventsLocationsY(end+1,:) = motionVectorY;
% %     
% %     singleTimeVecX = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(motionVectorX));
% %     timeThreshIdx = find(singleTimeVecX>(singleTimeVecX(brainMotionStart)+timeThresh),1);
% %     dispTimeThreshX(end+1) = motionVectorX(timeThreshIdx) - motionVectorX(brainMotionStart);
% %     idxToThreshXSingle = find((motionVectorX - motionVectorX(brainMotionStart)) > moveThresh & singleTimeVecX>singleTimeVecX(brainMotionStart) & singleTimeVecX<=3,1);
% %     if ~isempty(idxToThreshXSingle)
% %         timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
% %     end
% %     singleTimeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(motionVectorY));
% %     dispTimeThreshY(end+1) = (motionVectorY(timeThreshIdx) - motionVectorY(brainMotionStart))*-1;
% %     idxToThreshYSingle = find((motionVectorY - motionVectorY(brainMotionStart))*-1 > moveThresh & singleTimeVecY>singleTimeVecY(brainMotionStart) & singleTimeVecY<=3,1);
% %     if ~isempty(idxToThreshYSingle)
% %         timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
% %     end
% % end
% % h(13) = figure('Color','White');
% % subplot(2,1,1)
% % histfit(timeToThreshX,numBins,'kernel')
% % title(['Time for brain to displace laterally ' num2str(moveThresh) ' micrometers following brain motion start'])
% % xlabel('Time (s)')
% % xlim([-2 3])
% % ylim([0 25])
% % hold on
% % mu = mean(timeToThreshX);
% % sig = std(timeToThreshX);
% % plot([mu mu],[0 40],'k','LineWidth',2);
% % plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% % plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% % hold off
% % text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% % 
% % subplot(2,1,2)
% % histfit(timeToThreshY,numBins,'kernel')
% % title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following brain motion start'])
% % xlabel('Time (s)')
% % xlim([-2 3])
% % ylim([0 25])
% % hold on
% % mu = mean(timeToThreshY);
% % sig = std(timeToThreshY);
% % plot([mu mu],[0 40],'k','LineWidth',2);
% % plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% % plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% % hold off
% % text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% % 
% % h(50) = figure('Color','White');
% % subplot(2,1,1)
% % histfit(dispTimeThreshX,numBins,'kernel')
% % title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following brain motion start'])
% % xlabel('Displacement (\mum)')
% % xlim([-3 4])
% % ylim([0 40])
% % hold on
% % mu = mean(dispTimeThreshX);
% % sig = std(dispTimeThreshX);
% % plot([mu mu],[0 40],'k','LineWidth',2);
% % plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% % plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% % hold off
% % text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% % 
% % subplot(2,1,2)
% % histfit(dispTimeThreshY,numBins,'kernel')
% % title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following brain motion start'])
% % xlabel('Time (s)')
% % xlim([-3 4])
% % ylim([0 40])
% % hold on
% % mu = mean(dispTimeThreshY);
% % sig = std(dispTimeThreshY);
% % plot([mu mu],[0 40],'k','LineWidth',2);
% % plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% % plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% % hold off
% % text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
end