function combineMotionTracking(fileName,versionVec,mode)
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

if mode == 1
    [combinedTargetPositionX cIntFillPtsX] = getCIntMeanAndFillPts(targetPositionX,99);
    [combinedTargetPositionY cIntFillPtsY] = getCIntMeanAndFillPts(targetPositionY,99);
elseif mode == 2
    [combinedTargetPositionX cIntFillPtsX] = getCIntMedianAndFillPts(targetPositionX,99);
    [combinedTargetPositionY cIntFillPtsY] = getCIntMedianAndFillPts(targetPositionY,99);
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
sgf = sgolayfilt(targetPosition,3,13);

movementData.moveDist = moveDist;
movementData.velocity = velocity;
movementData.targetPosition = targetPosition;
movementData.targetPositionSGF = sgf;
movementData.cIntFillPtsX = cIntFillPtsX;
movementData.cIntFillPtsY = cIntFillPtsY;

save([fileName 'combined.mat'],'movementData');
plotMotionTracking([fileName 'combined.mat']);
end