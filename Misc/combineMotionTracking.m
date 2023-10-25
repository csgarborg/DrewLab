function combineMotionTracking(fileName,versionVec,mode,brainTF)
close all
targetPositionX = [];
targetPositionY = [];
for n = 1:length(versionVec)
    load([fileName num2str(versionVec(n)) '.mat'])
    if n > 1
        if size(movementData.targetPosition,1) > size(targetPositionX,2)
            movementData.targetPosition = movementData.targetPosition(1:size(targetPositionX,2),:);
        elseif size(movementData.targetPosition,1) < size(targetPositionX,2)
            targetPositionX = targetPositionX(:,1:size(movementData.targetPosition,1));
        end
    end
    targetPositionX(n,:) = movementData.targetPosition(:,1)';
    targetPositionY(n,:) = movementData.targetPosition(:,2)';
    if n < length(versionVec)
        clear movementData
    end
end

% if brainTF
%     cutoff = 7;
% else
%     cutoff = 3;
% end
% for m = 1:size(targetPositionX,2)
%     for n = 2:size(targetPositionX,1)
%         if abs(targetPositionX(n,m)) > cutoff
%             targetPositionX(n,m) = targetPositionX(n-1,m);
%         end
%     end
% end
% for m = 1:size(targetPositionY,2)
%     for n = 2:size(targetPositionY,1)
%         if abs(targetPositionY(n,m)) > cutoff
%             targetPositionY(n,m) = targetPositionY(n-1,m);
%         end
%     end
% end

if mode == 1
    [combinedTargetPositionX,cIntFillPtsX,meanCIX,stdCIX] = getCIntMeanAndFillPts(targetPositionX,95);
    [combinedTargetPositionY,cIntFillPtsY,meanCIY,stdCIY] = getCIntMeanAndFillPts(targetPositionY,95);
elseif mode == 2
    [combinedTargetPositionX,cIntFillPtsX,meanCIX,stdCIX] = getCIntMedianAndFillPts(targetPositionX,95);
    [combinedTargetPositionY,cIntFillPtsY,meanCIY,stdCIY] = getCIntMedianAndFillPts(targetPositionY,95);
else
    disp('Specify mode for combining positional data (1 = mean, 2 = median)')
    return
end

moveDist = [diff(combinedTargetPositionX)' diff(combinedTargetPositionY)'];
for n = 1:size(moveDist,1)
    velocity(n,1) = sqrt((moveDist(n,1))^2+(moveDist(n,2))^2)/movementData.secondsPerFrame;
end
targetPosition = [combinedTargetPositionX' combinedTargetPositionY'];

% sgf = sgolayfilt(targetPosition,3,17);
% sgf = sgolayfilt(targetPosition,3,13);

movementData.moveDist = moveDist;
movementData.velocity = velocity;
movementData.targetPosition = targetPosition;
% movementData.targetPositionSGF = sgf;
movementData.cIntFillPtsX = cIntFillPtsX;
movementData.cIntFillPtsY = cIntFillPtsY;
movementData.meanCIX = meanCIX;
movementData.meanCIY = meanCIY;
movementData.stdCIX = stdCIX;
movementData.stdCIY = stdCIY;

save([fileName 'combined.mat'],'movementData');
try
    if contains(fileName,'Layer')
        plotMotionTracking2Layer([fileName 'combined.mat']);
    else
        plotMotionTracking([fileName 'combined.mat']);
    end
catch
    display('Error in producing plots, data saved')
end
end