clear
close all
timeOverlapThresh = 2;
load('LTADataCell_FS.mat')
load('ETADataCell_FS.mat')
locDataCellOld = locDataCell;
goodLocEMG = 0;
badLoc = 0;
badEMG = 0;
goodOn = 0;
    badOn = 0;
    goodOff = 0;
    badOff = 0;
for i = 1:size(locDataCell,1)
    emgIdx = find(strcmp(EMGDataCell(:,1),locDataCell{i,1}));
    if isempty(locDataCell{i,2})
        locDataCell(i,2:end) = {[],NaN,NaN,[],NaN,NaN};
        badLoc = badLoc + 1;
        continue
    elseif isempty(EMGDataCell{emgIdx,2})
        locDataCell(i,2:end) = {[],NaN,NaN,[],NaN,NaN};
        badEMG = badEMG + 1;
        continue
    else
        goodLocEMG = goodLocEMG + 1;
    end
    load(locDataCell{i,1})
%     posL2 = [movementData.emgData(:,1), movementData.emgData(:,2)];
%     meanPosValL2 = [.5*(posL2(1:end-1,1) + posL2(2:end,1)),.5*(posL2(1:end-1,2) + posL2(2:end,2))];
%     posL2 = zipperVecs(posL2,meanPosValL2);
    emgEventsLocationsX = [];
%     motionEventsLocationsY = [];
    deleteRows = [];
    for n = 1:size(locDataCell{i,2},1)
        locTriggerTime = locDataCell{i,2}(n,2);
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
        locTriggerEMGIdx = find(min(abs(locTriggerTime - movementData.emgData(:,1))) == abs(locTriggerTime - movementData.emgData(:,1)));
        locDataCell{i,2}(n,4:6) = [locTriggerEMGIdx-60 locTriggerEMGIdx locTriggerEMGIdx+90];
        emgVector = movementData.emgData(locTriggerEMGIdx-60:locTriggerEMGIdx+90,2);
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
        emgEventsLocationsX(end+1,:) = emgVector - mean(emgVector(1:15)) + 1;
%         motionEventsLocationsY(end+1,:) = motionVectorY - motionVectorY(1);
    end
%     if size(emgEventsLocationsX,2) == 199
%         emgEventsLocationsX = emgEventsLocationsX(:,1:end-1);
%         motionEventsLocationsY = motionEventsLocationsY(:,1:end-1);
%     end
if isempty(emgEventsLocationsX)
    locDataCell{i,3} = NaN;
else
    locDataCell{i,3} = [linspace(-2,3,size(emgEventsLocationsX,2)); mean(emgEventsLocationsX,1)];
end
    locDataCell{i,2}(deleteRows,:) = [];
    
%     locDataCell{i,4}(2,:) = mean(motionEventsLocationsY);
    
    emgEventsLocationsX = [];
%     motionEventsLocationsY = [];
    deleteRows = [];
    for n = 1:size(locDataCell{i,5},1)
        locTriggerTime = locDataCell{i,5}(n,2);
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
        locTriggerEMGIdx = find(min(abs(locTriggerTime - movementData.emgData(:,1))) == abs(locTriggerTime - movementData.emgData(:,1)));
        locDataCell{i,5}(n,4:6) = [locTriggerEMGIdx-60 locTriggerEMGIdx locTriggerEMGIdx+90];
        emgVector = movementData.emgData(locTriggerEMGIdx-60:locTriggerEMGIdx+90,2);
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
        emgEventsLocationsX(end+1,:) = emgVector - emgVector(1) + 1;
%         motionEventsLocationsY(end+1,:) = motionVectorY - motionVectorY(1);
    end
%     if size(emgEventsLocationsX,2) == 199
%         emgEventsLocationsX = emgEventsLocationsX(:,1:end-1);
% %         motionEventsLocationsY = motionEventsLocationsY(:,1:end-1);
%     end
if isempty(emgEventsLocationsX)
    locDataCell{i,6} = NaN;
else
    locDataCell{i,6} = [linspace(-2,3,size(emgEventsLocationsX,2)); mean(emgEventsLocationsX,1)];
end
    locDataCell{i,5}(deleteRows,:) = [];
%     locDataCell{i,7}(2,:) = mean(motionEventsLocationsY);
end
disp('goodOn badOn goodOff badOff')
% [goodOn badOn goodOff badOff]
[badLoc badEMG goodLocEMG]
locDataCell(:,[4 7]) = [];
save('LTADataCellEMG_FS.mat','locDataCell')

% load('ETADataCell_FS.mat')
% for i = 1:size(EMGDataCell,1)
%     if isnan(EMGDataCell{i,3})
%         continue
%     end
%     load(['D' EMGDataCell{i,1}(2:40) 'd_Layer2_combined.mat'])
%     posL2 = [movementData.emgData(:,1), movementData.emgData(:,2)];
%     meanPosValL2 = [.5*(posL2(1:end-1,1) + posL2(2:end,1)),.5*(posL2(1:end-1,2) + posL2(2:end,2))];
%     posL2 = zipperVecs(posL2,meanPosValL2);
%     motionEventsLocationsX = [];
%     motionEventsLocationsY = [];
%     for n = 1:size(EMGDataCell{i,2},1)
%         emgVectorX = posL2(EMGDataCell{i,2}(n,4):EMGDataCell{i,2}(n,6),1);
%         if movementData.hemisphere == 2
%             emgVectorX = emgVectorX*-1;
%         end
%         motionVectorY = posL2(EMGDataCell{i,2}(n,4):EMGDataCell{i,2}(n,6),2);
%         if n > 1
%             if length(emgVectorX) > size(motionEventsLocationsX,2)
%                 emgVectorX = emgVectorX(1:size(motionEventsLocationsX,2));
%             elseif length(emgVectorX) < size(motionEventsLocationsX,2)
%                 motionEventsLocationsX = motionEventsLocationsX(:,1:length(emgVectorX));
%             end
%             if length(motionVectorY) > size(motionEventsLocationsY,2)
%                 motionVectorY = motionVectorY(1:size(motionEventsLocationsY,2));
%             elseif length(motionVectorY) < size(motionEventsLocationsY,2)
%                 motionEventsLocationsY = motionEventsLocationsY(:,1:length(motionVectorY));
%             end
%         end
%         motionEventsLocationsX(end+1,:) = emgVectorX - emgVectorX(1);
%         motionEventsLocationsY(end+1,:) = motionVectorY - motionVectorY(1);
%     end
%     if size(motionEventsLocationsX,2) == 199
%         motionEventsLocationsX = motionEventsLocationsX(:,1:end-1);
%         motionEventsLocationsY = motionEventsLocationsY(:,1:end-1);
%     end
%     EMGDataCell{i,3}(2,:) = mean(motionEventsLocationsX);
%     EMGDataCell{i,4}(2,:) = mean(motionEventsLocationsY);
%     
%     motionEventsLocationsX = [];
%     motionEventsLocationsY = [];
%     for n = 1:size(EMGDataCell{i,5},1)
%         emgVectorX = posL2(EMGDataCell{i,5}(n,4):EMGDataCell{i,5}(n,6),1);
%         motionVectorY = posL2(EMGDataCell{i,5}(n,4):EMGDataCell{i,5}(n,6),2);
%         if n > 1
%             if length(emgVectorX) > size(motionEventsLocationsX,2)
%                 emgVectorX = emgVectorX(1:size(motionEventsLocationsX,2));
%             elseif length(emgVectorX) < size(motionEventsLocationsX,2)
%                 motionEventsLocationsX = motionEventsLocationsX(:,1:length(emgVectorX));
%             end
%             if length(motionVectorY) > size(motionEventsLocationsY,2)
%                 motionVectorY = motionVectorY(1:size(motionEventsLocationsY,2));
%             elseif length(motionVectorY) < size(motionEventsLocationsY,2)
%                 motionEventsLocationsY = motionEventsLocationsY(:,1:length(motionVectorY));
%             end
%         end
%         motionEventsLocationsX(end+1,:) = emgVectorX - emgVectorX(1);
%         motionEventsLocationsY(end+1,:) = motionVectorY - motionVectorY(1);
%     end
%     if size(motionEventsLocationsX,2) == 199
%         motionEventsLocationsX = motionEventsLocationsX(:,1:end-1);
%         motionEventsLocationsY = motionEventsLocationsY(:,1:end-1);
%     end
%     EMGDataCell{i,6}(2,:) = mean(motionEventsLocationsX);
%     EMGDataCell{i,7}(2,:) = mean(motionEventsLocationsY);
% end
% save('ETADataCellSkull_FS.mat','EMGDataCell')

% function zipPos = zipperVecs(rawPos,meanPos)
% zipPos = [];
% for n = 1:size(rawPos,1)-1
%     zipPos(end+1,:) = [rawPos(n,1),rawPos(n,2)];
%     zipPos(end+1,:) = [meanPos(n,1),meanPos(n,2)];
% end
% zipPos(end+1,:) = [rawPos(end,1),rawPos(end,2)];
% end