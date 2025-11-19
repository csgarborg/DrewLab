load('LTADataCell_FS.mat')
for i = 1:size(locDataCell,1)
    if isnan(locDataCell{i,3})
        continue
    end
    load(['D' locDataCell{i,1}(2:40) 'd_Layer2_combined.mat'])
    posL2 = [movementData.targetPosition(:,1), movementData.targetPosition(:,2)];
    meanPosValL2 = [.5*(posL2(1:end-1,1) + posL2(2:end,1)),.5*(posL2(1:end-1,2) + posL2(2:end,2))];
    posL2 = zipperVecs(posL2,meanPosValL2);
    motionEventsLocationsX = [];
    motionEventsLocationsY = [];
    for n = 1:size(locDataCell{i,2},1)
        motionVectorX = posL2(locDataCell{i,2}(n,4):locDataCell{i,2}(n,6),1);
        if movementData.hemisphere == 2
            motionVectorX = motionVectorX*-1;
        end
        motionVectorY = posL2(locDataCell{i,2}(n,4):locDataCell{i,2}(n,6),2);
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
        motionEventsLocationsX(end+1,:) = motionVectorX - motionVectorX(1);
        motionEventsLocationsY(end+1,:) = motionVectorY - motionVectorY(1);
    end
    if size(motionEventsLocationsX,2) == 199
        motionEventsLocationsX = motionEventsLocationsX(:,1:end-1);
        motionEventsLocationsY = motionEventsLocationsY(:,1:end-1);
    end
    locDataCell{i,3}(2,:) = mean(motionEventsLocationsX);
    locDataCell{i,4}(2,:) = mean(motionEventsLocationsY);
    
    motionEventsLocationsX = [];
    motionEventsLocationsY = [];
    for n = 1:size(locDataCell{i,5},1)
        motionVectorX = posL2(locDataCell{i,5}(n,4):locDataCell{i,5}(n,6),1);
        motionVectorY = posL2(locDataCell{i,5}(n,4):locDataCell{i,5}(n,6),2);
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
        motionEventsLocationsX(end+1,:) = motionVectorX - motionVectorX(1);
        motionEventsLocationsY(end+1,:) = motionVectorY - motionVectorY(1);
    end
    if size(motionEventsLocationsX,2) == 199
        motionEventsLocationsX = motionEventsLocationsX(:,1:end-1);
        motionEventsLocationsY = motionEventsLocationsY(:,1:end-1);
    end
    locDataCell{i,6}(2,:) = mean(motionEventsLocationsX);
    locDataCell{i,7}(2,:) = mean(motionEventsLocationsY);
end
save('LTADataCellSkull_FS.mat','locDataCell')

load('ETADataCell_FS.mat')
for i = 1:size(EMGDataCell,1)
    if isnan(EMGDataCell{i,3})
        continue
    end
    load(['D' EMGDataCell{i,1}(2:40) 'd_Layer2_combined.mat'])
    posL2 = [movementData.targetPosition(:,1), movementData.targetPosition(:,2)];
    meanPosValL2 = [.5*(posL2(1:end-1,1) + posL2(2:end,1)),.5*(posL2(1:end-1,2) + posL2(2:end,2))];
    posL2 = zipperVecs(posL2,meanPosValL2);
    motionEventsLocationsX = [];
    motionEventsLocationsY = [];
    for n = 1:size(EMGDataCell{i,2},1)
        motionVectorX = posL2(EMGDataCell{i,2}(n,4):EMGDataCell{i,2}(n,6),1);
        if movementData.hemisphere == 2
            motionVectorX = motionVectorX*-1;
        end
        motionVectorY = posL2(EMGDataCell{i,2}(n,4):EMGDataCell{i,2}(n,6),2);
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
        motionEventsLocationsX(end+1,:) = motionVectorX - motionVectorX(1);
        motionEventsLocationsY(end+1,:) = motionVectorY - motionVectorY(1);
    end
    if size(motionEventsLocationsX,2) == 199
        motionEventsLocationsX = motionEventsLocationsX(:,1:end-1);
        motionEventsLocationsY = motionEventsLocationsY(:,1:end-1);
    end
    EMGDataCell{i,3}(2,:) = mean(motionEventsLocationsX);
    EMGDataCell{i,4}(2,:) = mean(motionEventsLocationsY);
    
    motionEventsLocationsX = [];
    motionEventsLocationsY = [];
    for n = 1:size(EMGDataCell{i,5},1)
        motionVectorX = posL2(EMGDataCell{i,5}(n,4):EMGDataCell{i,5}(n,6),1);
        motionVectorY = posL2(EMGDataCell{i,5}(n,4):EMGDataCell{i,5}(n,6),2);
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
        motionEventsLocationsX(end+1,:) = motionVectorX - motionVectorX(1);
        motionEventsLocationsY(end+1,:) = motionVectorY - motionVectorY(1);
    end
    if size(motionEventsLocationsX,2) == 199
        motionEventsLocationsX = motionEventsLocationsX(:,1:end-1);
        motionEventsLocationsY = motionEventsLocationsY(:,1:end-1);
    end
    EMGDataCell{i,6}(2,:) = mean(motionEventsLocationsX);
    EMGDataCell{i,7}(2,:) = mean(motionEventsLocationsY);
end
save('ETADataCellSkull_FS.mat','EMGDataCell')

function zipPos = zipperVecs(rawPos,meanPos)
zipPos = [];
for n = 1:size(rawPos,1)-1
    zipPos(end+1,:) = [rawPos(n,1),rawPos(n,2)];
    zipPos(end+1,:) = [meanPos(n,1),meanPos(n,2)];
end
zipPos(end+1,:) = [rawPos(end,1),rawPos(end,2)];
end