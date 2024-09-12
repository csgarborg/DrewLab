% Script to generate all figures for paper
clear;
close all;

%% figure 1e    210218_006 221207_010

combinedMovementDataBrain_1e = combineMotionTracking_FS('221207_010_Layer1_',1:3);
combinedMovementDataSkull_1e = combineMotionTracking_FS('221207_010_Layer2_',1:3);
plotMotionTrackingBrainAndSkull_FS(combinedMovementDataBrain_1e,combinedMovementDataSkull_1e);


combinedMovementData_1e = combineMotionTracking_FS('210218_006_',1:8);
plotMotionTrackingBrainOnly_FS(combinedMovementData_1e);

%% figure 2a

combinedMovementDataBrain_2a = combineMotionTracking_FS('220203_003_Layer1_',1:3);
combinedMovementDataSkull_2a = combineMotionTracking_FS('220203_003_Layer2_',1:3);
plotMotionTracking2P2LPCA_FS(combinedMovementDataBrain_2a,combinedMovementDataSkull_2a);

%% figure 2b

plotMovementQuiver_FS

%% figure 2c

combinedMovementDataBrain_2c = combineMotionTracking_FS('211216_002_Layer1_',1:5);
combinedMovementDataSkull_2c = combineMotionTracking_FS('211216_002_Layer2_',1:5);
load('211216_002_rawData.mat');
plotThermoCoherence_FS(combinedMovementDataBrain_2c,combinedMovementDataSkull_2c,rawData)

%% figure 2d

combinedMovementDataBrain_2d = combinedMovementDataBrain_2c;
combinedMovementDataSkull_2d = combinedMovementDataSkull_2c;
plotLocomotionXCorr_FS(combinedMovementDataBrain_2d,combinedMovementDataSkull_2d)

%% figure 3c

combinedMovementDataBrain_3c = combinedMovementDataBrain_2d;
combinedMovementDataSkull_3c = combinedMovementDataSkull_2d;
plotMovementBrainInSkull_FS(combinedMovementDataBrain_3c,combinedMovementDataSkull_3c,rawData)

%% figure 3d

combinedMovementDataBrain_3d = combinedMovementDataBrain_3c;
combinedMovementDataSkull_3d = combinedMovementDataSkull_3c;
plotEMGMovement_FS(combinedMovementDataBrain_3d,combinedMovementDataSkull_3d)

%% figure 3e

combinedMovementDataBrain_3e = combinedMovementDataBrain_3d;
combinedMovementDataSkull_3e = combinedMovementDataSkull_3d;
plotEMGXCorr_FS(combinedMovementDataBrain_3e,combinedMovementDataSkull_3e)

%% figure 3f

plotLocomotionTriggeredAvg_FS
plotEMGTriggeredAvg_FS

plotLocomotionTriggeredAvgSkull_FS
plotEMGTriggeredAvgSkull_FS

plotSqueezeTriggeredAvg_FS
plotSqueezeTriggeredAvgSkull_FS

plotRespTriggeredAvg_FS
plotRespTriggeredAvgSkull_FS

combinedMovementDataBrain_3f = combinedMovementDataBrain_3e;
plotLocomotionTriggeredAvgEMGSingleTrial_FS(combinedMovementDataBrain_3f)
plotLocomotionTriggeredAvgEMG_FS

plotCalibration2D_FS
plotCalibrationETL_FS

combinedMovementDataBrain_s4d = combineMotionTracking_FS('240611_001_Layer1_',1:4);
combinedMovementDataSkull_s4d = combineMotionTracking_FS('240611_001_Layer2_',1:4);
plotMotionTrackingBrainAndSkullResp_FS(combinedMovementDataBrain_s4d,combinedMovementDataSkull_s4d);

combinedMovementDataBrain_s4e = combineMotionTracking_FS('240612_002_Layer1_',1:3);
combinedMovementDataSkull_s4e = combineMotionTracking_FS('240612_002_Layer2_',1:3);
plotMotionTrackingBrainAndSkullRespComp_FS(combinedMovementDataBrain_s4e,combinedMovementDataSkull_s4e);

%% convert all figures to render for vectorization
h =  findobj('type','figure');
i = length(h);
for n = 1:i
    figure(n)
    set(gcf,'renderer','Painters')
end

%% SUBFUNCTIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMotionTrackingBrainAndSkull_FS(movementData,stationaryData)
[movementData.targetPosition,stationaryData.targetPosition] = interpBrainSkullMovement_FS(movementData,stationaryData);
if movementData.hemisphere == 2
    movementData.targetPosition(:,1) = movementData.targetPosition(:,1)*-1;
    stationaryData.targetPosition(:,1) = stationaryData.targetPosition(:,1)*-1;
end
movementData.targetPosition(:,2) = movementData.targetPosition(:,2)*-1;
stationaryData.targetPosition(:,2) = stationaryData.targetPosition(:,2)*-1;
movementData.secondsPerFrame = movementData.secondsPerFrame/2;
h(6) = figure('Color','White');
subplot(3,1,1)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'r')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'g')
hold off
title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
xlabel('Time (s)')
ylabel('X Position (\mum)')
grid on
axis([148 398 -3 6])
text(148,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(148,-3,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(3,1,2)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'r')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'g')
hold off
xlabel('Time (s)')
ylabel('Y Position (\mum)')
grid on
axis([148 398 -3 6])
text(148,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(148,-3,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(3,1,3)
plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)*2*pi*.06),'k')
xlabel('Time (s)')
ylabel('m/s')
grid on
axis([148 398 0 .3])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combinedMovementData = combineMotionTracking_FS(fileName,versionVec)
mode = 1;
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
    [combinedTargetPositionX,cIntFillPtsX,meanCIX,stdCIX] = getCIntMeanAndFillPts_FS(targetPositionX,95);
    [combinedTargetPositionY,cIntFillPtsY,meanCIY,stdCIY] = getCIntMeanAndFillPts_FS(targetPositionY,95);
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

combinedMovementData = movementData;
combinedMovementData.moveDist = moveDist;
combinedMovementData.velocity = velocity;
combinedMovementData.targetPosition = targetPosition;
combinedMovementData.cIntFillPtsX = cIntFillPtsX;
combinedMovementData.cIntFillPtsY = cIntFillPtsY;
combinedMovementData.meanCIX = meanCIX;
combinedMovementData.meanCIY = meanCIY;
combinedMovementData.stdCIX = stdCIX;
combinedMovementData.stdCIY = stdCIY;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [yMean,cIntFillPts,meanCI,stdCI] = getCIntMeanAndFillPts_FS(yMat,CIntPercentage)

N = size(yMat,1);                                      % Number of ‘Experiments’ In Data Set
if N < 2
    yMean = sgolayfilt(yMat,3,13);
else
    yMean = sgolayfilt(mean(yMat),3,13);                                    % Mean Of All Experiments At Each Value Of ‘x’
end
ySEM = std(yMat)/sqrt(N);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
cIntVal = (1 - (CIntPercentage/100)) * 0.5;
CI95 = tinv([cIntVal 1-cIntVal], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
cIntFillPts = sgolayfilt([(yMean + yCI95(1,:)) (flip(yMean) + flip(yCI95(2,:)))],3,13);
meanCI = mean(yCI95(2,:));
stdCI = std(yCI95(2,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMotionTrackingBrainOnly_FS(movementData)

subtitle = 'Figure 1e';
yLimitMotion = [-2 15];
yLimitLocomotion = [0 0.3];

h(1) = figure('Color','White');

x1 = subplot(3,1,1);
plot((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame,movementData.targetPosition(:,1),'r')
hold on;
f = fill([(1:size(movementData.targetPosition,1))*movementData.secondsPerFrame flip((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame)],movementData.cIntFillPtsX,'r','Linestyle','none');
set(f,'facea',[.2]);
hold off
title(['\fontsize{20pt}\bf{Figure 1e}'])
xlabel('Time (s)')
ylabel('X Position (\mum)')
grid on
if movementData.hemisphere == 1
    text(0,ceil(max([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)])),'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(0,floor(min([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)])),'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
else
    text(0,ceil(max([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)])),'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(0,floor(min([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)])),'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
end
ylim(yLimitMotion)

x2 = subplot(3,1,2);
plot((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame,-1*movementData.targetPosition(:,2),'b')
hold on;
f = fill([(1:size(movementData.targetPosition,1))*movementData.secondsPerFrame flip((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame)],-1*movementData.cIntFillPtsY,'r','Linestyle','none');
set(f,'facea',[.2]);
hold off
xlabel('Time (s)')
ylabel('Y Position (\mum)')
grid on
text(0,-floor(min([movementData.targetPosition(:,1);movementData.targetPosition(:,2)])),'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(0,-ceil(max([movementData.targetPosition(:,1);movementData.targetPosition(:,2)])),'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
ylim(yLimitMotion)

x3 = subplot(3,1,3);
plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)*2*pi*.06),'k')
xlabel('Time (s)')
ylabel('m/s')
grid on
ylim(yLimitLocomotion)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function targetPositionInSkull = combineBrainSkullMovement_FS(movementData,stationaryData)
posL1 = [movementData.targetPosition(:,1), movementData.targetPosition(:,2)];
meanPosValL1 = [.5*(posL1(1:end-1,1) + posL1(2:end,1)),.5*(posL1(1:end-1,2) + posL1(2:end,2))];
posL1 = zipperVecs_FS(posL1,meanPosValL1);
posL2 = [stationaryData.targetPosition(:,1), stationaryData.targetPosition(:,2)];
meanPosValL2 = [.5*(posL2(1:end-1,1) + posL2(2:end,1)),.5*(posL2(1:end-1,2) + posL2(2:end,2))];
posL2 = zipperVecs_FS(posL2,meanPosValL2);
if size(posL1,1) > size(posL2,1)
    posL1 = posL1(1:end-1,:);
    posL2 = [posL2(1,:); posL2(:,:)];
elseif size(posL1,1) < size(posL2,1)
    posL2 = posL2(1:end-1,:);
    posL1 = [posL1(1,:); posL1(:,:)];
end
targetPositionInSkull = [posL1(:,1)-posL2(:,1),-1*(posL1(:,2)-posL2(:,2))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [brainPosInterp,skullPosInterp] = interpBrainSkullMovement_FS(movementData,stationaryData)
posL1 = [movementData.targetPosition(:,1), movementData.targetPosition(:,2)];
meanPosValL1 = [.5*(posL1(1:end-1,1) + posL1(2:end,1)),.5*(posL1(1:end-1,2) + posL1(2:end,2))];
posL1 = zipperVecs_FS(posL1,meanPosValL1);
posL2 = [stationaryData.targetPosition(:,1), stationaryData.targetPosition(:,2)];
meanPosValL2 = [.5*(posL2(1:end-1,1) + posL2(2:end,1)),.5*(posL2(1:end-1,2) + posL2(2:end,2))];
posL2 = zipperVecs_FS(posL2,meanPosValL2);
if size(posL1,1) > size(posL2,1)
    posL1 = posL1(1:end-1,:);
    posL2 = [posL2(1,:); posL2(:,:)];
elseif size(posL1,1) < size(posL2,1)
    posL2 = posL2(1:end-1,:);
    posL1 = [posL1(1,:); posL1(:,:)];
end
brainPosInterp = posL1;
skullPosInterp = posL2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMotionTracking2P2LPCA_FS(movementData,stationaryData)
targetPositionInSkull = combineBrainSkullMovement_FS(movementData,stationaryData);
motionVec = pcaMotionAnalysis_FS(targetPositionInSkull);
h(2) = figure('Color','White');
s = scatter(targetPositionInSkull(:,1),targetPositionInSkull(:,2),10,'filled');
s.MarkerFaceAlpha = .1;
hold on
drawArrow_FS([0;0],[motionVec(1);motionVec(2)]);
hold off
axis equal square
axis([-6 6 -6 6])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
% title(['\fontsize{20pt}\bf{Position of Brain in Skull with PCA Vector}'])
title(['\fontsize{20pt}\bf{Figure 2a}'])
xlabel('\mum')
ylabel('\mum')
if movementData.hemisphere == 1
    text(6,0,'Lateral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(-6,0,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    text(0,6,'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(0,-6,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
else
    text(6,0,'Medial','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(-6,0,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    text(0,6,'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(0,-6,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function zipPos = zipperVecs_FS(rawPos,meanPos)
zipPos = [];
for n = 1:size(rawPos,1)-1
    zipPos(end+1,:) = [rawPos(n,1),rawPos(n,2)];
    zipPos(end+1,:) = [meanPos(n,1),meanPos(n,2)];
end
zipPos(end+1,:) = [rawPos(end,1),rawPos(end,2)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function motionVec = pcaMotionAnalysis_FS(positionData)

if ~exist('swapTF','var')
    swapTF = false;
end
if ~exist('vecReverseTF','var')
    vecReverseTF = false;
end
% Mean center data
positionDataMC(:,1) = positionData(:,1) - mean(positionData(:,1));
positionDataMC(:,2) = positionData(:,2) - mean(positionData(:,2));

% Calculate covariance matrix of data for rotation information data
C = cov(positionDataMC);

% Diagonalize covariance matrix to decorrelate new variables through
% rotation. V are eigenvectors belonging to diagonalized covariance matrix
% to express correlation between old and new data. D are eigenvalues that
% describe the variance within the calculated principal components. 

[V,D] = eig(C);

% Determine the principal component with largest variance to make the
% directional unit vec of movement
if D(1,1) > D(2,2)
    motionVec = [V(1,1) V(2,1)];
else
    motionVec = [V(1,2) V(2,2)];
end

% Calculate magnitude based on average of furthest 20% of movement data
% points from origin
numMags = ceil(size(positionData,1)*.2);
distanceVec = ((positionData(:,1).^2) + (positionData(:,2).^2)).^.5;
magnitude = mean(maxk(distanceVec,numMags));
motionVec = [motionVec * magnitude, motionVec * prctile(distanceVec,75), motionVec * prctile(distanceVec,95)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMovementQuiver_FS
load('movementDataLog_FS.mat')
xLoc = [];
yLoc = [];
UMic = [];
VMic = [];
a = {[] [] []};
b = {[] [] []};
c = {[] [] []};
d = {[] [] []};
e = {[] [] []};
[uniqueLocs,~,locIdx] = unique(cell2mat(moveDataMat(:,4:5)),'rows');
moveDataMatOrig = moveDataMat;
moveDataMat = {};
for n = 1:size(uniqueLocs)
    k = find(locIdx == n);
    if moveDataMatOrig{k(1),5} < -4000
        continue
    end
    moveDataMat(end+1,:) = moveDataMatOrig(k(1),:);
    %         xMoveVec = moveDataMatOrig{k,7}(1);
    %         yMoveVec = moveDataMatOrig{k,7}(2);
    xMoveSum = [];
    yMoveSum = [];
    for i = 1:length(k)
        xMoveSum(end+1) = moveDataMatOrig{k(i),7}(1);
        yMoveSum(end+1) = moveDataMatOrig{k(i),7}(2);
    end
    moveDataMat{end,7}(1) = mean(xMoveSum);
    moveDataMat{end,7}(2) = mean(yMoveSum);
end
skullMeanX = {[],[],[],[],[]};
skullMeanY = {[],[],[],[],[]};
skullStdX = {[],[],[],[],[]};
skullStdY = {[],[],[],[],[]};
brainMeanX = {[],[],[],[],[]};
brainMeanY = {[],[],[],[],[]};
brainStdX = {[],[],[],[],[]};
brainStdY = {[],[],[],[],[]};
h(3) = figure('Color','White');
hold on
for n = 1:size(moveDataMat,1)
    if strcmp(moveDataMat{n,3},'21a')
        arrowColor = 'r';
        a{1}(end+1) = sqrt(moveDataMat{n,7}(1)^2 + moveDataMat{n,7}(2)^2);
        catIdx = 0;
    elseif strcmp(moveDataMat{n,3},'27a')
        b{1}(end+1) = sqrt(moveDataMat{n,7}(1)^2 + moveDataMat{n,7}(2)^2);
        arrowColor = 'm';
        catIdx = 1;
    elseif strcmp(moveDataMat{n,3},'28a')
        arrowColor = 'g';
        c{1}(end+1) = sqrt(moveDataMat{n,7}(1)^2 + moveDataMat{n,7}(2)^2);
        catIdx = 2;
    elseif strcmp(moveDataMat{n,3},'39a')
        arrowColor = 'b';
        d{1}(end+1) = sqrt(moveDataMat{n,7}(1)^2 + moveDataMat{n,7}(2)^2);
        catIdx = 3;
    else
        arrowColor = 'k';
        e{1}(end+1) = sqrt(moveDataMat{n,7}(1)^2 + moveDataMat{n,7}(2)^2);
        catIdx = 4;
    end
    plot(moveDataMat{n,4},-moveDataMat{n,5},[arrowColor 'o'],'MarkerSize',5)
    quiver(moveDataMat{n,4},-moveDataMat{n,5},moveDataMat{n,7}(1),moveDataMat{n,7}(2),150,'Color',arrowColor,'LineWidth',2.5,'MaxHeadSize',40)
    skullMeanX{catIdx+1}(end+1) = moveDataMat{n,8};
    skullMeanY{catIdx+1}(end+1) = moveDataMat{n,9};
    skullStdX{catIdx+1}(end+1) = moveDataMat{n,10};
    skullStdY{catIdx+1}(end+1) = moveDataMat{n,11};
    brainMeanX{catIdx+1}(end+1) = moveDataMat{n,12};
    brainMeanY{catIdx+1}(end+1) = moveDataMat{n,13};
    brainStdX{catIdx+1}(end+1) = moveDataMat{n,14};
    brainStdY{catIdx+1}(end+1) = moveDataMat{n,15};
end
xlim([-8000 8000])
ylim([-8000 8000])
axis square
plot(0,0,'kx','MarkerSize',12)
plot(0,-2600,'kx','MarkerSize',12)
quiver(-2500,2500,3,0,150,'LineWidth',2.5,'MaxHeadSize',40)
title(['Figure 2b' 10 'Rostral, \mum'])
ylabel('Left, \mum')
xlabel('Caudal, \mum')
text(55,55,'Bregma')
text(55,-2500,'Lambda')
text(-2500,2300,'3 \mum')
rectangle('Position',[-2650 2000 800 800])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotThermoCoherence_FS(movementData,stationaryData,rawData)
% params.tapers = [10,19];   % Tapers [n, 2n - 1]
params.tapers = [15,29];
params.pad = 1;
params.Fs = 20;
params.fpass = [0.01,10];   % Pass band [0, nyquist]
params.trialave = 1;
params.err = [2,0.05];
% you want to make sure 'LH_restData' and 'RH_restData' have a mean of 0. You can use multiple trials just make sure the matrix is oriented properly

% input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
targetPositionInSkull = combineBrainSkullMovement_FS(movementData,stationaryData);
movementData.targetPosition = targetPositionInSkull;
% load('D:\21-12-16_MouseExp\211216_002_processed_Layer1_combined.mat');
Thermo = rawData(:,[1 2]); % load raw thermo data
resp = filterThermoDataHz_FS(Thermo,10000); % filter and resample thermo data from 10000 to 20 Hz
% resp = movementData.thermoData;
resp(:,2) = resp(:,2) - mean(resp(:,2));
brain = movementData.targetPosition;
brain = brain - mean(brain);
tBrain = ((1:size(brain,1)).*(1/params.Fs))-params.Fs;
if size(resp,1) > size(brain,1)
    resp = resp(1:size(brain,1),:);
elseif size(resp,1) < size(brain,1)
    brain = brain(1:size(resp,1),:);
    tBrain = tBrain(1:size(resp,1));
end
brainX = [tBrain' brain(:,1)];
brainY = [tBrain' brain(:,2)];
[Cx,~,~,~,~,fx,confCx,~,cErrx] = coherencyc_FS(resp(:,2),brainX(:,2),params);
[Cy,~,~,~,~,fy,confCy,~,cErry] = coherencyc_FS(resp(:,2),brainY(:,2),params);
[PowerX,HzX,ErrorX,PowerY,HzY,ErrorY,PowerT,HzT,ErrorT] = motionSpectrumAnalysisHz_FS(brain,resp(:,2));
% cErrx = [cErrx(1,:) flip(cErrx(2,:))];
% cErry = [cErry(1,:) flip(cErry(2,:))];

% figure(1)
% subplot(3,1,1)
% PowerX = downsample(PowerX,3);
% HzX = downsample(HzX,3);
% semilogx(HzX,PowerX,'k')
% xlabel('Frequency (Hz)')
% ylabel('Power')
% title('Brain Motion Spectrum (Lateral/Medial)')
% % hold on
% % f = fill([HzX flip(HzX)],[ErrorX(1,:) flip(ErrorX(2,:))],'r','Linestyle','none');
% % set(f,'facea',[.2]);
% % hold off
% subplot(3,1,2)
% PowerT = downsample(PowerT,3);
% HzT = downsample(HzT,3);
% semilogx(HzT,PowerT,'k')
% xlabel('Frequency (Hz)')
% ylabel('Power')
% title('Thermocouple Spectrum')
% % hold on
% % f = fill([HzT flip(HzT)],[ErrorT(1,:) flip(ErrorT(2,:))],'r','Linestyle','none');
% % set(f,'facea',[.2]);
% % hold off
% subplot(3,1,3)
% fx = downsample(fx,3);
% Cx = downsample(Cx,3);
% semilogx(fx,Cx);
% ylim([0 1])
% xlabel('Frequency (Hz)')
% ylabel('Coherence')
% hold on
% % f = fill([fx flip(fx)],[cErrx(1,:) flip(cErrx(2,:))],'r','Linestyle','none');
% % set(f,'facea',[.2]);
% yline(confCx);
% hold off
% text(.1,.8,['Conf = ' num2str(confCx)])
% title('Coherence')

h(4) = figure('Color','White');
subplot(3,1,1)
PowerY = downsample(PowerY,3);
HzY = downsample(HzY,3);
semilogx(HzY,PowerY,'k')
xlabel('Frequency (Hz)')
ylabel('Power')
title(['Figure 2c' 10 'Brain Motion Spectrum (Rostral/Caudal)'])
% hold on
% f = fill([HzY flip(HzY)],[ErrorY(1,:) flip(ErrorY(2,:))],'r','Linestyle','none');
% set(f,'facea',[.2]);
% hold off
subplot(3,1,2)
semilogx(HzT,PowerT,'k')
xlabel('Frequency (Hz)')
ylabel('Power')
title('Thermocouple Spectrum')
% hold on
% f = fill([HzT flip(HzT)],[ErrorT(1,:) flip(ErrorT(2,:))],'r','Linestyle','none');
% set(f,'facea',[.2]);
% hold off
subplot(3,1,3)
Cy = downsample(Cy,3);
fy = downsample(fy,3);
semilogx(fy,Cy);
ylim([0 1])
xlabel('Frequency (Hz)')
ylabel('Coherence')
hold on
% f = fill([fy flip(fy)],[cErry(1,:) flip(cErry(2,:))],'r','Linestyle','none');
% set(f,'facea',[.2]);
yline(confCx);
hold off
text(.1,.8,['Conf = ' num2str(confCy)])
title('Coherence')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function procData = filterThermoDataHz_FS(Thermo,sampleRate)

rawThermo = Thermo(:,2);
t = Thermo(:,1);
rawThermo=filloutliers(rawThermo, 'center','mean', 'ThresholdFactor', 4);%this takes out any large artifacts
% process EKG data
% fpass = [300,(sampleRate/2)-1];   % Hz
fpass = [2,40];   % Hz
trialDuration_sec = t(end);   % read this variable in from your data in seconds
analogSamplingRate = sampleRate;   % Hz - change if yours is different
dsFs = 20;   % Hz - downsampled frequency
analogExpectedLength = trialDuration_sec*analogSamplingRate;
trimmedThermo = rawThermo(1:min(analogExpectedLength,length(rawThermo)));
[z,p,k] = butter(5,fpass./(analogSamplingRate/2));
[sos,g] = zp2sos(z,p,k);
filtEKG = filtfilt(sos,g,trimmedThermo - mean(trimmedThermo));
% kernelWidth = 0.5;
kernelWidth = 0.005;
smoothingKernel = gausswin(kernelWidth*analogSamplingRate)/sum(gausswin(kernelWidth*analogSamplingRate));
% EKGPwr = log10(conv(filtEKG.^2,smoothingKernel,'same'));
EKGPwr = conv(filtEKG,smoothingKernel,'same');
% EKGPwr = filtEKG.^2;
resampEKG = resample(EKGPwr,dsFs,analogSamplingRate);
procEKG = resampEKG;   % save this as your final array
procT = 0:1/dsFs:t(end);
% procT = procT(1:end-1);
if length(procT) > length(procEKG)
    procT = procT(1:length(procEKG));
elseif length(procT) < length(procEKG)
    procEKG = procEKG(1:length(procT));
end
% plot(procT,procEKG)
% title('Select start and stop time (t) to use to determine mean baseline')
% tVals = ginput(2);
% close
% i = tVals(1,1) <= procT & procT <= tVals(2,1);
% baseline = mean(procEKG(i));
% procEKG = procEKG - (baseline - 1);
procData = [procT',procEKG];

% plot(procData(:,1),procData(:,2))
% title('Select upper or lower limit')
% limitVal = ginput(1);
% close
% for n = 1:size(procData,1)
%     if procData(n,2) > abs(limitVal(2))
%         procData(n,2) = abs(limitVal(2));
%     elseif procData(n,2) < -abs(limitVal(2))
%         procData(n,2) = -abs(limitVal(2));
%     end
% end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [PowerX,HzX,ErrorX,PowerY,HzY,ErrorY,PowerT,HzT,ErrorT] = motionSpectrumAnalysisHz_FS(positionData,thermo, plotTF)
dt = 1/20;
tapers = [2 3];
params.Fs = 1/dt;
params.tapers = tapers;
params.fpass = [0.025 7.5];
params.err = [2 0.05]; %use jack-knife resampling confidence intervals p = 0.05
% (could also try multitaper frequency domain bootstrapping (MFDB) for less noisy CI - but jackknife used by so many people in NVC)

[PowerX, HzX, ErrorX] = mtspectrumc_FS(positionData(:,1),params);
[PowerY, HzY, ErrorY] = mtspectrumc_FS(positionData(:,2),params);
[PowerT, HzT, ErrorT] = mtspectrumc_FS(thermo,params);

if ~exist('plotTF','var')
    plotTF = false;
end

if plotTF
    h(1) = figure('Color','White');
    semilogy(HzX,PowerX,'k')
    hold on
    f = fill([HzX flip(HzX)],[ErrorX(1,:) flip(ErrorX(2,:))],'r','Linestyle','none');
    set(f,'facea',[.2]);
    hold off
    title(['\fontsize{20pt}\bf{X Position Frequency Domain}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Frequency (Hz)')
    ylabel('Power')
    
    h(2) = figure('Color','White');
    semilogy(HzY,PowerY,'k')
    hold on
    f = fill([HzY flip(HzY)],[ErrorY(1,:) flip(ErrorY(2,:))],'r','Linestyle','none');
    set(f,'facea',[.2]);
    hold off
    title(['\fontsize{20pt}\bf{Y Position Frequency Domain}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Frequency (Hz)')
    ylabel('Power')
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotLocomotionXCorr_FS(movementData,stationaryData)
targetPositionInSkull = combineBrainSkullMovement_FS(movementData,stationaryData);
movementData.secondsPerFrame = movementData.secondsPerFrame/2;
movementData.targetPosition = targetPositionInSkull;
movement_time=movementData.secondsPerFrame*(1:length(movementData.targetPosition));
movementData.locDataInterp=zeros(size(movementData.targetPosition));
movementData.locDataInterp(:,1)=movement_time;
movementData.locDataInterp(:,2)=interp1(movementData.ballData(:,1),abs(movementData.ballData(:,2)),movement_time,'linear');

h(5) = figure('Color','White');
maxlag=500;
xc_1=xcorr(detrend(movementData.locDataInterp(1:(end-100),2))', detrend(movementData.targetPosition(1:(end-100),1))',maxlag,'coeff');
xc_2=xcorr(detrend(movementData.locDataInterp(1:(end-100),2))', detrend(movementData.targetPosition(1:(end-100),2))',maxlag,'coeff');
plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_1,'b')
hold on
plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_2,'g')
hold off
title(['Figure 2d' 10 'Brain Motion and Locomotion Cross-Correlation'])
ylabel('Noramlized Cross-Correlation')
xlabel('Lags (s)')
legend({'L/M Brain Motion','R/C Brain Motion'})
ylim([0 1])
axes('Position',[.2 .7 .2 .2])
box on
plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_1,'b')
hold on
plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_2,'g')
hold off
axis([-2 2 .6 .7])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMovementBrainInSkull_FS(movementData,stationaryData,rawData)
targetPositionInSkull = combineBrainSkullMovement_FS(movementData,stationaryData);
movementData.secondsPerFrame = movementData.secondsPerFrame/2;
h(6) = figure('Color','White');
subplot(4,1,1)
plot([1:size(targetPositionInSkull,1)]*movementData.secondsPerFrame,targetPositionInSkull(:,1),'r')
title(['Figure 3c' 10 '\fontsize{20pt}\bf{Position of Brain in Skull}'])
xlabel('Time (s)')
ylabel('X Position (\mum)')
grid on
axis([0 600 -2 6])
if movementData.hemisphere == 1
    text(0,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(0,-2,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
else
    text(0,6,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(0,-2,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
end
subplot(4,1,2)
plot([1:size(targetPositionInSkull,1)]*movementData.secondsPerFrame,targetPositionInSkull(:,2),'g')
xlabel('Time (s)')
ylabel('Y Position (\mum)')
grid on
axis([0 600 -2 6])
text(0,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(0,-2,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(4,1,3)
if all(movementData.emgData(:,2) == 0)
    title('\fontsize{20pt}\bf{No EMG Data}')
else
    semilogy(movementData.emgData(:,1),movementData.emgData(:,2),'b')
    xlabel('Time (s)')
    ylabel('EMG Power (A.U.)')
    grid on
    xlim([0 600])
end
subplot(4,1,4)
plot(rawData(:,1),rawData(:,3),'k')
axis([0 600 -1 1])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotEMGMovement_FS(movementData,stationaryData)
movementData.targetPosition = combineBrainSkullMovement_FS(movementData,stationaryData);
movementData.secondsPerFrame = movementData.secondsPerFrame/2;
best_lag=round(1/movementData.secondsPerFrame);
movement_time=movementData.secondsPerFrame*(1:length(movementData.targetPosition));
movementData.emgDataInterp=zeros(size(movementData.targetPosition));
movementData.emgDataInterp(:,1)=movement_time;
movementData.emgDataInterp(:,2)=interp1(movementData.emgData(:,1),movementData.emgData(:,2),movement_time,'linear');
emg_bins=.5:.1:3.5;
pixel_bins=-3:.5:7;
h(7) = figure('Color','White');
subplot(2,2,2)
scatter(movementData.emgDataInterp(1:(end-best_lag),2),movementData.targetPosition(best_lag:end-1,2))
title('Figure 3d')
xlabel('Log EMG Power')
ylabel('Brain Shift, \mum')
text(3.5,6,'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
text(3.5,-2,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
ylim([-2 6])
xlim([0.5 3.5])

subplot(2,2,4)
hh2=histogram2(movementData.emgDataInterp(1:(end-best_lag),2),movementData.targetPosition(best_lag:end-1,2), emg_bins, pixel_bins,...
    'DisplayStyle','tile','ShowEmptyBins','on');
imagesc(emg_bins, pixel_bins,log(hh2.Values'+1))
axis xy
colormap hot
xlabel('Log EMG Power')
ylabel('Brain Shift, \mum')
ylim([-2 6])
xlim([0.5 3.5])

subplot(2,2,1)
scatter(movementData.emgDataInterp(1:(end-best_lag),2),movementData.targetPosition(best_lag:end-1,1))
title('Figure 3d')
xlabel('Log EMG Power')
ylabel('Brain Shift, \mum')
text(3.5,6,'Lateral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
text(3.5,-2,'Medial','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
ylim([-2 6])
xlim([0.5 3.5])

subplot(2,2,3)
hh2=histogram2(movementData.emgDataInterp(1:(end-best_lag),2),movementData.targetPosition(best_lag:end-1,1), emg_bins, pixel_bins,...
    'DisplayStyle','tile','ShowEmptyBins','on');
imagesc(emg_bins, pixel_bins,log(hh2.Values'+1))
axis xy
colormap hot
xlabel('Log EMG Power')
ylabel('Brain Shift, \mum')
ylim([-2 6])
xlim([0.5 3.5])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotEMGXCorr_FS(movementData,stationaryData)
targetPositionInSkull = combineBrainSkullMovement_FS(movementData,stationaryData);
movementData.secondsPerFrame = movementData.secondsPerFrame/2;
movementData.targetPosition = targetPositionInSkull;
movement_time=movementData.secondsPerFrame*(1:length(movementData.targetPosition));
movementData.emgDataInterp=zeros(size(movementData.targetPosition));
movementData.emgDataInterp(:,1)=movement_time;
movementData.emgDataInterp(:,2)=interp1(movementData.emgData(:,1),abs(movementData.emgData(:,2)),movement_time,'linear');

h(8) = figure('Color','White');
maxlag=500;
xc_1=xcorr(detrend(movementData.emgDataInterp(1:(end-100),2))', detrend(movementData.targetPosition(1:(end-100),1))',maxlag,'coeff');
xc_2=xcorr(detrend(movementData.emgDataInterp(1:(end-100),2))', detrend(movementData.targetPosition(1:(end-100),2))',maxlag,'coeff');
plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_1,'b')
hold on
plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_2,'g')
hold off
title(['Figure 3e' 10 'Brain Motion and EMG Cross-Correlation'])
ylabel('Noramlized Cross-Correlation')
xlabel('Lags (s)')
legend({'L/M Brain Motion','R/C Brain Motion'})
ylim([0 1])
axes('Position',[.2 .7 .2 .2])
box on
plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_1,'b')
hold on
plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_2,'g')
hold off
axis([-2 2 .7 .8])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotLocomotionTriggeredAvg_FS
motionEventsLocationsX = [];
motionEventsLocationsY = [];
timeToThreshX = [];
timeToThreshY = [];
dispTimeThreshX = [];
dispTimeThreshY = [];
moveThresh = .75;
timeThresh = 1.5;
numBins = 20;
load('LTADataCell_FS.mat')
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
    dispTimeThreshX(end+1) = motionVectorX((find(singleTimeVecX>timeThresh,1))) - motionVectorX((find(singleTimeVecX>0,1)));
    idxToThreshXSingle = find((motionVectorX - motionVectorX(find(singleTimeVecX>0,1))) > moveThresh & singleTimeVecX>0 & singleTimeVecX<=3,1);
    if ~isempty(idxToThreshXSingle)
        timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
    end
    singleTimeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(motionVectorY));
    dispTimeThreshY(end+1) =  (motionVectorY((find(singleTimeVecY>timeThresh,1))) - motionVectorY((find(singleTimeVecY>0,1))))*-1;
    idxToThreshYSingle = find((motionVectorY - motionVectorY(find(singleTimeVecY>0,1)))*-1 > moveThresh & singleTimeVecY>0 & singleTimeVecY<=3,1);
    if ~isempty(idxToThreshYSingle)
        timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
    end
end
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts_FS(motionEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts_FS(motionEventsLocationsY,90);
meanY = -1*meanY;
cIntFillPtsY = -1*cIntFillPtsY;
timeVecX = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(meanX));
timeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(meanY));

stdWindowSize = 5;
for n = stdWindowSize+1:length(meanY)
    if std(meanY(n-stdWindowSize:n)) > .01
        brainMotionStart = n;
        break
    end 
end

h(9) = figure('Color','White');
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
plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
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
    plot(timeVecY,-1*motionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
end
plot(timeVecY(brainMotionStart),meanY(brainMotionStart),'rx')
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
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts_FS(stopMotionEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts_FS(stopMotionEventsLocationsY,90);
meanY = -1*meanY;
cIntFillPtsY = -1*cIntFillPtsY;
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
    plot(timeVecY,-1*stopMotionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
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

h(10) = figure('Color','White');
subplot(2,1,1)
histfit(timeToThreshX,numBins,'kernel')
title(['Time for brain to displace laterally ' num2str(moveThresh) ' micrometers following locomotion trigger'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 25])
hold on
mu = mean(timeToThreshX);
sig = std(timeToThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(timeToThreshY,numBins,'kernel')
title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following locomotion trigger'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 25])
hold on
mu = mean(timeToThreshY);
sig = std(timeToThreshY);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

h(50) = figure('Color','White');
subplot(2,1,1)
histfit(dispTimeThreshX,numBins,'kernel')
title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following locomotion trigger'])
xlabel('Displacement (\mum)')
xlim([-3 4])
ylim([0 40])
hold on
mu = mean(dispTimeThreshX);
sig = std(dispTimeThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(dispTimeThreshY,numBins,'kernel')
title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following locomotion trigger'])
xlabel('Time (s)')
xlim([-3 4])
ylim([0 40])
hold on
mu = mean(dispTimeThreshY);
sig = std(dispTimeThreshY);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

timeToThreshX = [];
timeToThreshY = [];
dispTimeThreshX = [];
dispTimeThreshY = [];
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
    timeThreshIdx = find(singleTimeVecX>(singleTimeVecX(brainMotionStart)+timeThresh),1);
    dispTimeThreshX(end+1) = motionVectorX(timeThreshIdx) - motionVectorX(brainMotionStart);
    idxToThreshXSingle = find((motionVectorX - motionVectorX(brainMotionStart)) > moveThresh & singleTimeVecX>singleTimeVecX(brainMotionStart) & singleTimeVecX<=3,1);
    if ~isempty(idxToThreshXSingle)
        timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
    end
    singleTimeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(motionVectorY));
    dispTimeThreshY(end+1) = (motionVectorY(timeThreshIdx) - motionVectorY(brainMotionStart))*-1;
    idxToThreshYSingle = find((motionVectorY - motionVectorY(brainMotionStart))*-1 > moveThresh & singleTimeVecY>singleTimeVecY(brainMotionStart) & singleTimeVecY<=3,1);
    if ~isempty(idxToThreshYSingle)
        timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
    end
end
h(13) = figure('Color','White');
subplot(2,1,1)
histfit(timeToThreshX,numBins,'kernel')
title(['Time for brain to displace laterally ' num2str(moveThresh) ' micrometers following brain motion start'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 25])
hold on
mu = mean(timeToThreshX);
sig = std(timeToThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(timeToThreshY,numBins,'kernel')
title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following brain motion start'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 25])
hold on
mu = mean(timeToThreshY);
sig = std(timeToThreshY);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

h(50) = figure('Color','White');
subplot(2,1,1)
histfit(dispTimeThreshX,numBins,'kernel')
title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following brain motion start'])
xlabel('Displacement (\mum)')
xlim([-3 4])
ylim([0 40])
hold on
mu = mean(dispTimeThreshX);
sig = std(dispTimeThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(dispTimeThreshY,numBins,'kernel')
title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following brain motion start'])
xlabel('Time (s)')
xlim([-3 4])
ylim([0 40])
hold on
mu = mean(dispTimeThreshY);
sig = std(dispTimeThreshY);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotEMGTriggeredAvg_FS
motionEventsLocationsX = [];
motionEventsLocationsY = [];
timeToThreshX = [];
timeToThreshY = [];
dispTimeThreshX = [];
dispTimeThreshY = [];
moveThresh = .75;
timeThresh = 1.5;
numBins = 20;
load('ETADataCell_FS.mat')
for n = 1:size(EMGDataCell)
    if isnan(EMGDataCell{n,3})
        continue
    end
    motionVectorX = EMGDataCell{n,3}(2,:);
    motionVectorY = EMGDataCell{n,4}(2,:);
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
    
    singleTimeVecX = linspace(round(EMGDataCell{1,2}(1,1)-EMGDataCell{1,2}(1,2)),round(EMGDataCell{1,2}(1,3)-EMGDataCell{1,2}(1,2)),length(motionVectorX));
    dispTimeThreshX(end+1) = motionVectorX((find(singleTimeVecX>timeThresh,1))) - motionVectorX((find(singleTimeVecX>0,1)));
    idxToThreshXSingle = find((motionVectorX - motionVectorX(find(singleTimeVecX>0,1))) > moveThresh & singleTimeVecX>0 & singleTimeVecX<=3,1);
    if ~isempty(idxToThreshXSingle)
        timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
    end
    singleTimeVecY = linspace(round(EMGDataCell{1,2}(1,1)-EMGDataCell{1,2}(1,2)),round(EMGDataCell{1,2}(1,3)-EMGDataCell{1,2}(1,2)),length(motionVectorY));
    dispTimeThreshY(end+1) =  (motionVectorY((find(singleTimeVecY>timeThresh,1))) - motionVectorY((find(singleTimeVecY>0,1))))*-1;
    idxToThreshYSingle = find((motionVectorY - motionVectorY(find(singleTimeVecY>0,1)))*-1 > moveThresh & singleTimeVecY>0 & singleTimeVecY<=3,1);
    if ~isempty(idxToThreshYSingle)
        timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
    end
end
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts_FS(motionEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts_FS(motionEventsLocationsY,90);
meanY = -1*meanY;
cIntFillPtsY = -1*cIntFillPtsY;
timeVecX = linspace(round(EMGDataCell{1,2}(1,1)-EMGDataCell{1,2}(1,2)),round(EMGDataCell{1,2}(1,3)-EMGDataCell{1,2}(1,2)),length(meanX));
timeVecY = linspace(round(EMGDataCell{1,2}(1,1)-EMGDataCell{1,2}(1,2)),round(EMGDataCell{1,2}(1,3)-EMGDataCell{1,2}(1,2)),length(meanY));

stdWindowSize = 5;
for n = stdWindowSize+1:length(meanY)
    if std(meanY(n-stdWindowSize:n)) > .01
        brainMotionStart = n;
        break
    end 
end

h(9) = figure('Color','White');
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
plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean Motion During EMG Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
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
    plot(timeVecY,-1*motionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
end
plot(timeVecY(brainMotionStart),meanY(brainMotionStart),'rx')
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
for n = 1:size(EMGDataCell)
    if isnan(EMGDataCell{n,6})
        continue
    end
    stopMotionVectorX = EMGDataCell{n,6}(2,:);
    stopMotionVectorY = EMGDataCell{n,7}(2,:);
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
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts_FS(stopMotionEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts_FS(stopMotionEventsLocationsY,90);
meanY = -1*meanY;
cIntFillPtsY = -1*cIntFillPtsY;
timeVecX = linspace(round(EMGDataCell{1,5}(1,1)-EMGDataCell{1,5}(1,2)),round(EMGDataCell{1,5}(1,3)-EMGDataCell{1,5}(1,2)),length(meanX));
timeVecY = linspace(round(EMGDataCell{1,5}(1,1)-EMGDataCell{1,5}(1,2)),round(EMGDataCell{1,5}(1,3)-EMGDataCell{1,5}(1,2)),length(meanY));

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
title(['\fontsize{20pt}\bf{Mean Motion During Stopping EMG Events, n = ' num2str(size(stopMotionEventsLocationsX,1)) '}'])
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
    plot(timeVecY,-1*stopMotionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
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

h(10) = figure('Color','White');
subplot(2,1,1)
histfit(timeToThreshX,numBins,'kernel')
title(['Time for brain to displace laterally ' num2str(moveThresh) ' micrometers following EMG trigger'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 25])
hold on
mu = mean(timeToThreshX);
sig = std(timeToThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(timeToThreshY,numBins,'kernel')
title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following EMG trigger'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 25])
hold on
mu = mean(timeToThreshY);
sig = std(timeToThreshY);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

h(50) = figure('Color','White');
subplot(2,1,1)
histfit(dispTimeThreshX,numBins,'kernel')
title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following EMG trigger'])
xlabel('Displacement (\mum)')
xlim([-3 4])
ylim([0 40])
hold on
mu = mean(dispTimeThreshX);
sig = std(dispTimeThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(dispTimeThreshY,numBins,'kernel')
title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following EMG trigger'])
xlabel('Time (s)')
xlim([-3 4])
ylim([0 40])
hold on
mu = mean(dispTimeThreshY);
sig = std(dispTimeThreshY);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

timeToThreshX = [];
timeToThreshY = [];
dispTimeThreshX = [];
dispTimeThreshY = [];
for n = 1:size(EMGDataCell)
    if isnan(EMGDataCell{n,3})
        continue
    end
    motionVectorX = EMGDataCell{n,3}(2,:);
    motionVectorY = EMGDataCell{n,4}(2,:);
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
    
    singleTimeVecX = linspace(round(EMGDataCell{1,2}(1,1)-EMGDataCell{1,2}(1,2)),round(EMGDataCell{1,2}(1,3)-EMGDataCell{1,2}(1,2)),length(motionVectorX));
    timeThreshIdx = find(singleTimeVecX>(singleTimeVecX(brainMotionStart)+timeThresh),1);
    dispTimeThreshX(end+1) = motionVectorX(timeThreshIdx) - motionVectorX(brainMotionStart);
    idxToThreshXSingle = find((motionVectorX - motionVectorX(brainMotionStart)) > moveThresh & singleTimeVecX>singleTimeVecX(brainMotionStart) & singleTimeVecX<=3,1);
    if ~isempty(idxToThreshXSingle)
        timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
    end
    singleTimeVecY = linspace(round(EMGDataCell{1,2}(1,1)-EMGDataCell{1,2}(1,2)),round(EMGDataCell{1,2}(1,3)-EMGDataCell{1,2}(1,2)),length(motionVectorY));
    dispTimeThreshY(end+1) = (motionVectorY(timeThreshIdx) - motionVectorY(brainMotionStart))*-1;
    idxToThreshYSingle = find((motionVectorY - motionVectorY(brainMotionStart))*-1 > moveThresh & singleTimeVecY>singleTimeVecY(brainMotionStart) & singleTimeVecY<=3,1);
    if ~isempty(idxToThreshYSingle)
        timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
    end
end
h(13) = figure('Color','White');
subplot(2,1,1)
histfit(timeToThreshX,numBins,'kernel')
title(['Time for brain to displace laterally ' num2str(moveThresh) ' micrometers following brain motion start'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 25])
hold on
mu = mean(timeToThreshX);
sig = std(timeToThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(timeToThreshY,numBins,'kernel')
title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following brain motion start'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 25])
hold on
mu = mean(timeToThreshY);
sig = std(timeToThreshY);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

h(50) = figure('Color','White');
subplot(2,1,1)
histfit(dispTimeThreshX,numBins,'kernel')
title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following brain motion start'])
xlabel('Displacement (\mum)')
xlim([-3 4])
ylim([0 40])
hold on
mu = mean(dispTimeThreshX);
sig = std(dispTimeThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(dispTimeThreshY,numBins,'kernel')
title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following brain motion start'])
xlabel('Time (s)')
xlim([-3 4])
ylim([0 40])
hold on
mu = mean(dispTimeThreshY);
sig = std(dispTimeThreshY);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotLocomotionTriggeredAvgSkull_FS
motionEventsLocationsX = [];
motionEventsLocationsY = [];
timeToThreshX = [];
timeToThreshY = [];
dispTimeThreshX = [];
dispTimeThreshY = [];
moveThresh = .1;
timeThresh = 1.5;
numBins = 20;
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
    
    singleTimeVecX = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(motionVectorX));
    dispTimeThreshX(end+1) = motionVectorX((find(singleTimeVecX>timeThresh,1))) - motionVectorX((find(singleTimeVecX>0,1)));
    idxToThreshXSingle = find(abs((motionVectorX - motionVectorX(find(singleTimeVecX>0,1)))) > moveThresh & singleTimeVecX>0 & singleTimeVecX<=3,1);
    if ~isempty(idxToThreshXSingle)
        timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
    end
    singleTimeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(motionVectorY));
    dispTimeThreshY(end+1) =  (motionVectorY((find(singleTimeVecY>timeThresh,1))) - motionVectorY((find(singleTimeVecY>0,1))))*-1;
    idxToThreshYSingle = find(abs((motionVectorY - motionVectorY(find(singleTimeVecY>0,1)))*-1) > moveThresh & singleTimeVecY>0 & singleTimeVecY<=3,1);
    if ~isempty(idxToThreshYSingle)
        timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
    end
end
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts_FS(motionEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts_FS(motionEventsLocationsY,90);
meanY = -1*meanY;
cIntFillPtsY = -1*cIntFillPtsY;
timeVecX = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(meanX));
timeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(meanY));

% stdWindowSize = 5;
% for n = stdWindowSize+1:length(meanY)
%     if std(meanY(n-stdWindowSize:n)) > .01
%         brainMotionStart = n;
%         break
%     end 
% end

h(9) = figure('Color','White');
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
% plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean Skull Motion During Locomotion Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Skull Shift (\mum)')
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
    plot(timeVecY,-1*motionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
end
% plot(timeVecY(brainMotionStart),meanY(brainMotionStart),'rx')
hold off
text(3,-3,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Rostral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
xlabel('Time (s)')
ylabel('\Delta Skull Shift (\mum)')
ylim([-3 3])
xlim([-2 3])
grid on
clear movementData

stopMotionEventsLocationsX = [];
stopMotionEventsLocationsY = [];
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
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts_FS(stopMotionEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts_FS(stopMotionEventsLocationsY,90);
meanY = -1*meanY;
cIntFillPtsY = -1*cIntFillPtsY;
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
title(['\fontsize{20pt}\bf{Mean Skull Motion During Stopping Locomotion Events, n = ' num2str(size(stopMotionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Skull Shift (\mum)')
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
    plot(timeVecY,-1*stopMotionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
end
hold off
text(3,-3,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Rostral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
xlabel('Time (s)')
ylabel('\Delta Skull Shift (\mum)')
ylim([-3 3])
xlim([-2 3])
grid on
clear movementData

h(10) = figure('Color','White');
subplot(2,1,1)
histfit(timeToThreshX,numBins,'kernel')
title(['Time for skull to displace laterally ' num2str(moveThresh) ' micrometers following locomotion trigger'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 30])
hold on
mu = mean(timeToThreshX);
sig = std(timeToThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(timeToThreshY,numBins,'kernel')
title(['Time for skull to displace rostrally ' num2str(moveThresh) ' micrometers following locomotion trigger'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 30])
hold on
mu = mean(timeToThreshY);
sig = std(timeToThreshY);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

h(50) = figure('Color','White');
subplot(2,1,1)
histfit(dispTimeThreshX,numBins,'kernel')
title(['Lateral displacement of skull after ' num2str(timeThresh) ' s following locomotion trigger'])
xlabel('Displacement (\mum)')
xlim([-2 3])
ylim([0 60])
hold on
mu = mean(dispTimeThreshX);
sig = std(dispTimeThreshX);
plot([mu mu],[0 60],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(dispTimeThreshY,numBins,'kernel')
title(['Rostral displacement of skull after ' num2str(timeThresh) ' s following locomotion trigger'])
xlabel('Displacement (\mum)')
xlim([-2 3])
ylim([0 60])
hold on
mu = mean(dispTimeThreshY);
sig = std(dispTimeThreshY);
plot([mu mu],[0 60],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

% timeToThreshX = [];
% timeToThreshY = [];
% dispTimeThreshX = [];
% dispTimeThreshY = [];
% for n = 1:size(locDataCell)
%     if isnan(locDataCell{n,3})
%         continue
%     end
%     motionVectorX = locDataCell{n,3}(2,:);
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
%     motionEventsLocationsX(end+1,:) = motionVectorX;
%     motionEventsLocationsY(end+1,:) = motionVectorY;
%     
%     singleTimeVecX = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(motionVectorX));
%     timeThreshIdx = find(singleTimeVecX>(singleTimeVecX(brainMotionStart)+timeThresh),1);
%     dispTimeThreshX(end+1) = motionVectorX(timeThreshIdx) - motionVectorX(brainMotionStart);
%     idxToThreshXSingle = find((motionVectorX - motionVectorX(brainMotionStart)) > moveThresh & singleTimeVecX>singleTimeVecX(brainMotionStart) & singleTimeVecX<=3,1);
%     if ~isempty(idxToThreshXSingle)
%         timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
%     end
%     singleTimeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(motionVectorY));
%     dispTimeThreshY(end+1) = (motionVectorY(timeThreshIdx) - motionVectorY(brainMotionStart))*-1;
%     idxToThreshYSingle = find((motionVectorY - motionVectorY(brainMotionStart))*-1 > moveThresh & singleTimeVecY>singleTimeVecY(brainMotionStart) & singleTimeVecY<=3,1);
%     if ~isempty(idxToThreshYSingle)
%         timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
%     end
% end
% h(13) = figure('Color','White');
% subplot(2,1,1)
% histfit(timeToThreshX,numBins,'kernel')
% title(['Time for brain to displace laterally ' num2str(moveThresh) ' micrometers following brain motion start'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 25])
% hold on
% mu = mean(timeToThreshX);
% sig = std(timeToThreshX);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% subplot(2,1,2)
% histfit(timeToThreshY,numBins,'kernel')
% title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following brain motion start'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 25])
% hold on
% mu = mean(timeToThreshY);
% sig = std(timeToThreshY);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% h(50) = figure('Color','White');
% subplot(2,1,1)
% histfit(dispTimeThreshX,numBins,'kernel')
% title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following brain motion start'])
% xlabel('Displacement (\mum)')
% xlim([-3 4])
% ylim([0 40])
% hold on
% mu = mean(dispTimeThreshX);
% sig = std(dispTimeThreshX);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% subplot(2,1,2)
% histfit(dispTimeThreshY,numBins,'kernel')
% title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following brain motion start'])
% xlabel('Time (s)')
% xlim([-3 4])
% ylim([0 40])
% hold on
% mu = mean(dispTimeThreshY);
% sig = std(dispTimeThreshY);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotEMGTriggeredAvgSkull_FS
motionEventsLocationsX = [];
motionEventsLocationsY = [];
timeToThreshX = [];
timeToThreshY = [];
dispTimeThreshX = [];
dispTimeThreshY = [];
moveThresh = .1;
timeThresh = 1.5;
numBins = 20;
load('ETADataCellSkull_FS.mat')
for n = 1:size(EMGDataCell)
    if isnan(EMGDataCell{n,3})
        continue
    end
    motionVectorX = EMGDataCell{n,3}(2,:);
    motionVectorY = EMGDataCell{n,4}(2,:);
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
    
    singleTimeVecX = linspace(round(EMGDataCell{1,2}(1,1)-EMGDataCell{1,2}(1,2)),round(EMGDataCell{1,2}(1,3)-EMGDataCell{1,2}(1,2)),length(motionVectorX));
    dispTimeThreshX(end+1) = motionVectorX((find(singleTimeVecX>timeThresh,1))) - motionVectorX((find(singleTimeVecX>0,1)));
    idxToThreshXSingle = find(abs((motionVectorX - motionVectorX(find(singleTimeVecX>0,1)))) > moveThresh & singleTimeVecX>0 & singleTimeVecX<=3,1);
    if ~isempty(idxToThreshXSingle)
        timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
    end
    singleTimeVecY = linspace(round(EMGDataCell{1,2}(1,1)-EMGDataCell{1,2}(1,2)),round(EMGDataCell{1,2}(1,3)-EMGDataCell{1,2}(1,2)),length(motionVectorY));
    dispTimeThreshY(end+1) =  (motionVectorY((find(singleTimeVecY>timeThresh,1))) - motionVectorY((find(singleTimeVecY>0,1))))*-1;
    idxToThreshYSingle = find(abs((motionVectorY - motionVectorY(find(singleTimeVecY>0,1)))*-1) > moveThresh & singleTimeVecY>0 & singleTimeVecY<=3,1);
    if ~isempty(idxToThreshYSingle)
        timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
    end
end
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts_FS(motionEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts_FS(motionEventsLocationsY,90);
meanY = -1*meanY;
cIntFillPtsY = -1*cIntFillPtsY;
timeVecX = linspace(round(EMGDataCell{1,2}(1,1)-EMGDataCell{1,2}(1,2)),round(EMGDataCell{1,2}(1,3)-EMGDataCell{1,2}(1,2)),length(meanX));
timeVecY = linspace(round(EMGDataCell{1,2}(1,1)-EMGDataCell{1,2}(1,2)),round(EMGDataCell{1,2}(1,3)-EMGDataCell{1,2}(1,2)),length(meanY));

% stdWindowSize = 5;
% for n = stdWindowSize+1:length(meanY)
%     if std(meanY(n-stdWindowSize:n)) > .01
%         brainMotionStart = n;
%         break
%     end 
% end

h(9) = figure('Color','White');
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
% plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean Skull Motion During EMG Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Skull Shift (\mum)')
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
    plot(timeVecY,-1*motionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
end
% plot(timeVecY(brainMotionStart),meanY(brainMotionStart),'rx')
hold off
text(3,-3,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Rostral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
xlabel('Time (s)')
ylabel('\Delta Skull Shift (\mum)')
ylim([-3 3])
xlim([-2 3])
grid on
clear movementData

stopMotionEventsLocationsX = [];
stopMotionEventsLocationsY = [];
for n = 1:size(EMGDataCell)
    if isnan(EMGDataCell{n,6})
        continue
    end
    stopMotionVectorX = EMGDataCell{n,6}(2,:);
    stopMotionVectorY = EMGDataCell{n,7}(2,:);
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
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts_FS(stopMotionEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts_FS(stopMotionEventsLocationsY,90);
meanY = -1*meanY;
cIntFillPtsY = -1*cIntFillPtsY;
timeVecX = linspace(round(EMGDataCell{1,5}(1,1)-EMGDataCell{1,5}(1,2)),round(EMGDataCell{1,5}(1,3)-EMGDataCell{1,5}(1,2)),length(meanX));
timeVecY = linspace(round(EMGDataCell{1,5}(1,1)-EMGDataCell{1,5}(1,2)),round(EMGDataCell{1,5}(1,3)-EMGDataCell{1,5}(1,2)),length(meanY));

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
title(['\fontsize{20pt}\bf{Mean Skull Motion During Stopping EMG Events, n = ' num2str(size(stopMotionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Skull Shift (\mum)')
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
    plot(timeVecY,-1*stopMotionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
end
hold off
text(3,-3,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Rostral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
xlabel('Time (s)')
ylabel('\Delta Skull Shift (\mum)')
ylim([-3 3])
xlim([-2 3])
grid on
clear movementData

h(10) = figure('Color','White');
subplot(2,1,1)
histfit(timeToThreshX,numBins,'kernel')
title(['Time for skull to displace laterally ' num2str(moveThresh) ' micrometers following EMG trigger'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 30])
hold on
mu = mean(timeToThreshX);
sig = std(timeToThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(timeToThreshY,numBins,'kernel')
title(['Time for skull to displace rostrally ' num2str(moveThresh) ' micrometers following EMG trigger'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 30])
hold on
mu = mean(timeToThreshY);
sig = std(timeToThreshY);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

h(50) = figure('Color','White');
subplot(2,1,1)
histfit(dispTimeThreshX,numBins,'kernel')
title(['Lateral displacement of skull after ' num2str(timeThresh) ' s following EMG trigger'])
xlabel('Displacement (\mum)')
xlim([-3 4])
ylim([0 60])
hold on
mu = mean(dispTimeThreshX);
sig = std(dispTimeThreshX);
plot([mu mu],[0 60],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(dispTimeThreshY,numBins,'kernel')
title(['Rostral displacement of skull after ' num2str(timeThresh) ' s following EMG trigger'])
xlabel('Displacement (\mum)')
xlim([-3 4])
ylim([0 60])
hold on
mu = mean(dispTimeThreshY);
sig = std(dispTimeThreshY);
plot([mu mu],[0 60],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

% timeToThreshX = [];
% timeToThreshY = [];
% dispTimeThreshX = [];
% dispTimeThreshY = [];
% for n = 1:size(EMGDataCell)
%     if isnan(EMGDataCell{n,3})
%         continue
%     end
%     motionVectorX = EMGDataCell{n,3}(2,:);
%     motionVectorY = EMGDataCell{n,4}(2,:);
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
%     motionEventsLocationsX(end+1,:) = motionVectorX;
%     motionEventsLocationsY(end+1,:) = motionVectorY;
%     
%     singleTimeVecX = linspace(round(EMGDataCell{1,2}(1,1)-EMGDataCell{1,2}(1,2)),round(EMGDataCell{1,2}(1,3)-EMGDataCell{1,2}(1,2)),length(motionVectorX));
%     timeThreshIdx = find(singleTimeVecX>(singleTimeVecX(brainMotionStart)+timeThresh),1);
%     dispTimeThreshX(end+1) = motionVectorX(timeThreshIdx) - motionVectorX(brainMotionStart);
%     idxToThreshXSingle = find((motionVectorX - motionVectorX(brainMotionStart)) > moveThresh & singleTimeVecX>singleTimeVecX(brainMotionStart) & singleTimeVecX<=3,1);
%     if ~isempty(idxToThreshXSingle)
%         timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
%     end
%     singleTimeVecY = linspace(round(EMGDataCell{1,2}(1,1)-EMGDataCell{1,2}(1,2)),round(EMGDataCell{1,2}(1,3)-EMGDataCell{1,2}(1,2)),length(motionVectorY));
%     dispTimeThreshY(end+1) = (motionVectorY(timeThreshIdx) - motionVectorY(brainMotionStart))*-1;
%     idxToThreshYSingle = find((motionVectorY - motionVectorY(brainMotionStart))*-1 > moveThresh & singleTimeVecY>singleTimeVecY(brainMotionStart) & singleTimeVecY<=3,1);
%     if ~isempty(idxToThreshYSingle)
%         timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
%     end
% end
% h(13) = figure('Color','White');
% subplot(2,1,1)
% histfit(timeToThreshX,numBins,'kernel')
% title(['Time for brain to displace laterally ' num2str(moveThresh) ' micrometers following brain motion start'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 25])
% hold on
% mu = mean(timeToThreshX);
% sig = std(timeToThreshX);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% subplot(2,1,2)
% histfit(timeToThreshY,numBins,'kernel')
% title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following brain motion start'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 25])
% hold on
% mu = mean(timeToThreshY);
% sig = std(timeToThreshY);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% h(50) = figure('Color','White');
% subplot(2,1,1)
% histfit(dispTimeThreshX,numBins,'kernel')
% title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following brain motion start'])
% xlabel('Displacement (\mum)')
% xlim([-3 4])
% ylim([0 40])
% hold on
% mu = mean(dispTimeThreshX);
% sig = std(dispTimeThreshX);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% subplot(2,1,2)
% histfit(dispTimeThreshY,numBins,'kernel')
% title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following brain motion start'])
% xlabel('Time (s)')
% xlim([-3 4])
% ylim([0 40])
% hold on
% mu = mean(dispTimeThreshY);
% sig = std(dispTimeThreshY);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotSqueezeTriggeredAvg_FS
motionEventsLocationsX = [];
motionEventsLocationsY = [];
timeToThreshX = [];
timeToThreshY = [];
dispTimeThreshX = [];
dispTimeThreshY = [];
moveThresh = .75;
timeThresh = 1.5;
numBins = 20;
load('squeezeDataCell_FS.mat')
for n = 1:size(squeezeDataCell)
    if isnan(squeezeDataCell{n,3})
        continue
    end
    motionVectorX = squeezeDataCell{n,3}(2,:);
    motionVectorY = squeezeDataCell{n,4}(2,:);
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
    
    singleTimeVecX = linspace(round(squeezeDataCell{1,2}(1,1)-squeezeDataCell{1,2}(1,2)),round(squeezeDataCell{1,2}(1,3)-squeezeDataCell{1,2}(1,2)),length(motionVectorX));
    dispTimeThreshX(end+1) = motionVectorX((find(singleTimeVecX>timeThresh,1))) - motionVectorX((find(singleTimeVecX>0,1)));
    idxToThreshXSingle = find((motionVectorX - motionVectorX(find(singleTimeVecX>0,1))) > moveThresh & singleTimeVecX>0 & singleTimeVecX<=3,1);
    if ~isempty(idxToThreshXSingle)
        timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
    end
    singleTimeVecY = linspace(round(squeezeDataCell{1,2}(1,1)-squeezeDataCell{1,2}(1,2)),round(squeezeDataCell{1,2}(1,3)-squeezeDataCell{1,2}(1,2)),length(motionVectorY));
    dispTimeThreshY(end+1) =  (motionVectorY((find(singleTimeVecY>timeThresh,1))) - motionVectorY((find(singleTimeVecY>0,1))))*-1;
    idxToThreshYSingle = find((motionVectorY - motionVectorY(find(singleTimeVecY>0,1)))*-1 > moveThresh & singleTimeVecY>0 & singleTimeVecY<=3,1);
    if ~isempty(idxToThreshYSingle)
        timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
    end
end
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts_FS(motionEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts_FS(motionEventsLocationsY,90);
meanY = -1*meanY;
cIntFillPtsY = -1*cIntFillPtsY;
timeVecX = linspace(round(squeezeDataCell{1,2}(1,1)-squeezeDataCell{1,2}(1,2)),round(squeezeDataCell{1,2}(1,3)-squeezeDataCell{1,2}(1,2)),length(meanX));
timeVecY = linspace(round(squeezeDataCell{1,2}(1,1)-squeezeDataCell{1,2}(1,2)),round(squeezeDataCell{1,2}(1,3)-squeezeDataCell{1,2}(1,2)),length(meanY));

stdWindowSize = 5;
for n = stdWindowSize+1:length(meanY)
    if std(meanY(n-stdWindowSize:n)) > .01
        brainMotionStart = n;
        break
    end 
end

h(9) = figure('Color','White');
subplot(2,1,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
plot(timeVecX,meanX,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
set(f,'facea',[.2]);
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[0,0,1,0.1])
end
plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
f = fill([0 2 2 0],[4.9 4.9 -.9 -.9],'g','Linestyle','none','FaceAlpha',0.1);
hold off
text(5,-1,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(5,5,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean Motion During Squeeze Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Brian Shift (\mum)')
ylim([-1 5])
xlim([-2 6])
grid on

subplot(2,1,2)
plot(timeVecY,meanY,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
set(f,'facea',[.2]);
for n = 1:size(motionEventsLocationsY,1)
    plot(timeVecY,-1*motionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
end
plot(timeVecY(brainMotionStart),meanY(brainMotionStart),'rx')
f = fill([0 2 2 0],[4.9 4.9 -.9 -.9],'g','Linestyle','none','FaceAlpha',0.1);
hold off
text(5,-1,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(5,5,'Rostral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
xlabel('Time (s)')
ylabel('\Delta Brian Shift (\mum)')
ylim([-1 5])
xlim([-2 6])
grid on
clear movementData

h(10) = figure('Color','White');
subplot(2,1,1)
histfit(timeToThreshX,numBins,'kernel')
title(['Time for brain to displace laterally ' num2str(moveThresh) ' micrometers following squeeze'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 25])
hold on
mu = mean(timeToThreshX);
sig = std(timeToThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(timeToThreshY,numBins,'kernel')
title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following squeeze'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 25])
hold on
mu = mean(timeToThreshY);
sig = std(timeToThreshY);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

h(50) = figure('Color','White');
subplot(2,1,1)
histfit(dispTimeThreshX,numBins,'kernel')
title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following squeeze'])
xlabel('Displacement (\mum)')
xlim([-3 4])
ylim([0 40])
hold on
mu = mean(dispTimeThreshX);
sig = std(dispTimeThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(dispTimeThreshY,numBins,'kernel')
title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following squeeze'])
xlabel('Time (s)')
xlim([-3 4])
ylim([0 40])
hold on
mu = mean(dispTimeThreshY);
sig = std(dispTimeThreshY);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

timeToThreshX = [];
timeToThreshY = [];
dispTimeThreshX = [];
dispTimeThreshY = [];
for n = 1:size(squeezeDataCell)
    if isnan(squeezeDataCell{n,3})
        continue
    end
    motionVectorX = squeezeDataCell{n,3}(2,:);
    motionVectorY = squeezeDataCell{n,4}(2,:);
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
    
    singleTimeVecX = linspace(round(squeezeDataCell{1,2}(1,1)-squeezeDataCell{1,2}(1,2)),round(squeezeDataCell{1,2}(1,3)-squeezeDataCell{1,2}(1,2)),length(motionVectorX));
    timeThreshIdx = find(singleTimeVecX>(singleTimeVecX(brainMotionStart)+timeThresh),1);
    dispTimeThreshX(end+1) = motionVectorX(timeThreshIdx) - motionVectorX(brainMotionStart);
    idxToThreshXSingle = find((motionVectorX - motionVectorX(brainMotionStart)) > moveThresh & singleTimeVecX>singleTimeVecX(brainMotionStart) & singleTimeVecX<=3,1);
    if ~isempty(idxToThreshXSingle)
        timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
    end
    singleTimeVecY = linspace(round(squeezeDataCell{1,2}(1,1)-squeezeDataCell{1,2}(1,2)),round(squeezeDataCell{1,2}(1,3)-squeezeDataCell{1,2}(1,2)),length(motionVectorY));
    dispTimeThreshY(end+1) = (motionVectorY(timeThreshIdx) - motionVectorY(brainMotionStart))*-1;
    idxToThreshYSingle = find((motionVectorY - motionVectorY(brainMotionStart))*-1 > moveThresh & singleTimeVecY>singleTimeVecY(brainMotionStart) & singleTimeVecY<=3,1);
    if ~isempty(idxToThreshYSingle)
        timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
    end
end
h(13) = figure('Color','White');
subplot(2,1,1)
histfit(timeToThreshX,numBins,'kernel')
title(['Time for brain to displace laterally ' num2str(moveThresh) ' micrometers following brain motion start'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 25])
hold on
mu = mean(timeToThreshX);
sig = std(timeToThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(timeToThreshY,numBins,'kernel')
title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following brain motion start'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 25])
hold on
mu = mean(timeToThreshY);
sig = std(timeToThreshY);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

h(50) = figure('Color','White');
subplot(2,1,1)
histfit(dispTimeThreshX,numBins,'kernel')
title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following brain motion start'])
xlabel('Displacement (\mum)')
xlim([-3 4])
ylim([0 40])
hold on
mu = mean(dispTimeThreshX);
sig = std(dispTimeThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(dispTimeThreshY,numBins,'kernel')
title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following brain motion start'])
xlabel('Time (s)')
xlim([-3 4])
ylim([0 40])
hold on
mu = mean(dispTimeThreshY);
sig = std(dispTimeThreshY);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotSqueezeTriggeredAvgSkull_FS
motionEventsLocationsX = [];
motionEventsLocationsY = [];
timeToThreshX = [];
timeToThreshY = [];
dispTimeThreshX = [];
dispTimeThreshY = [];
moveThresh = .1;
timeThresh = 1.5;
numBins = 20;
load('squeezeDataCellSkull_FS.mat')
for n = 1:size(squeezeDataCell)
    if isnan(squeezeDataCell{n,3})
        continue
    end
    motionVectorX = squeezeDataCell{n,3}(2,:);
    motionVectorY = squeezeDataCell{n,4}(2,:);
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
    
    singleTimeVecX = linspace(round(squeezeDataCell{1,2}(1,1)-squeezeDataCell{1,2}(1,2)),round(squeezeDataCell{1,2}(1,3)-squeezeDataCell{1,2}(1,2)),length(motionVectorX));
    dispTimeThreshX(end+1) = motionVectorX((find(singleTimeVecX>timeThresh,1))) - motionVectorX((find(singleTimeVecX>0,1)));
    idxToThreshXSingle = find(abs((motionVectorX - motionVectorX(find(singleTimeVecX>0,1)))) > moveThresh & singleTimeVecX>0 & singleTimeVecX<=3,1);
    if ~isempty(idxToThreshXSingle)
        timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
    end
    singleTimeVecY = linspace(round(squeezeDataCell{1,2}(1,1)-squeezeDataCell{1,2}(1,2)),round(squeezeDataCell{1,2}(1,3)-squeezeDataCell{1,2}(1,2)),length(motionVectorY));
    dispTimeThreshY(end+1) =  (motionVectorY((find(singleTimeVecY>timeThresh,1))) - motionVectorY((find(singleTimeVecY>0,1))))*-1;
    idxToThreshYSingle = find(abs((motionVectorY - motionVectorY(find(singleTimeVecY>0,1)))*-1) > moveThresh & singleTimeVecY>0 & singleTimeVecY<=3,1);
    if ~isempty(idxToThreshYSingle)
        timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
    end
end
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts_FS(motionEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts_FS(motionEventsLocationsY,90);
meanY = -1*meanY;
cIntFillPtsY = -1*cIntFillPtsY;
timeVecX = linspace(round(squeezeDataCell{1,2}(1,1)-squeezeDataCell{1,2}(1,2)),round(squeezeDataCell{1,2}(1,3)-squeezeDataCell{1,2}(1,2)),length(meanX));
timeVecY = linspace(round(squeezeDataCell{1,2}(1,1)-squeezeDataCell{1,2}(1,2)),round(squeezeDataCell{1,2}(1,3)-squeezeDataCell{1,2}(1,2)),length(meanY));

% stdWindowSize = 5;
% for n = stdWindowSize+1:length(meanY)
%     if std(meanY(n-stdWindowSize:n)) > .01
%         brainMotionStart = n;
%         break
%     end 
% end

h(9) = figure('Color','White');
subplot(2,1,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
plot(timeVecX,meanX,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
set(f,'facea',[.2]);
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[0,0,1,0.1])
end
% plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
f = fill([0 2 2 0],[4.9 4.9 -.9 -.9],'g','Linestyle','none','FaceAlpha',0.1);
hold off
text(5,-1,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(5,5,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean Skull Motion During Squeeze Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Skull Shift (\mum)')
ylim([-1 5])
xlim([-2 6])
grid on

subplot(2,1,2)
plot(timeVecY,meanY,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
set(f,'facea',[.2]);
for n = 1:size(motionEventsLocationsY,1)
    plot(timeVecY,-1*motionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
end
% plot(timeVecY(brainMotionStart),meanY(brainMotionStart),'rx')
f = fill([0 2 2 0],[4.9 4.9 -.9 -.9],'g','Linestyle','none','FaceAlpha',0.1);
hold off
text(5,-1,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(5,5,'Rostral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
xlabel('Time (s)')
ylabel('\Delta Skull Shift (\mum)')
ylim([-1 5])
xlim([-2 6])
grid on
clear movementData

h(10) = figure('Color','White');
subplot(2,1,1)
histfit(timeToThreshX,numBins,'kernel')
title(['Time for skull to displace laterally ' num2str(moveThresh) ' micrometers following squeeze'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 30])
hold on
mu = mean(timeToThreshX);
sig = std(timeToThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(timeToThreshY,numBins,'kernel')
title(['Time for skull to displace rostrally ' num2str(moveThresh) ' micrometers following squeeze'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 30])
hold on
mu = mean(timeToThreshY);
sig = std(timeToThreshY);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

h(50) = figure('Color','White');
subplot(2,1,1)
histfit(dispTimeThreshX,numBins,'kernel')
title(['Lateral displacement of skull after ' num2str(timeThresh) ' s following squeeze'])
xlabel('Displacement (\mum)')
xlim([-2 3])
ylim([0 60])
hold on
mu = mean(dispTimeThreshX);
sig = std(dispTimeThreshX);
plot([mu mu],[0 60],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(dispTimeThreshY,numBins,'kernel')
title(['Rostral displacement of skull after ' num2str(timeThresh) ' s following squeeze'])
xlabel('Displacement (\mum)')
xlim([-2 3])
ylim([0 60])
hold on
mu = mean(dispTimeThreshY);
sig = std(dispTimeThreshY);
plot([mu mu],[0 60],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotRespTriggeredAvg_FS
motionEventsLocationsX = [];
motionEventsLocationsY = [];
timeToThreshX = [];
timeToThreshY = [];
dispTimeThreshX = [];
dispTimeThreshY = [];
moveThresh = .5;
timeThresh = .5;
numBins = 20;
load('respDataCell_FS.mat')
for n = 1:size(respDataCell)
    if isnan(respDataCell{n,3})
        continue
    end
    motionVectorX = respDataCell{n,3}(2,:);
    motionVectorY = respDataCell{n,4}(2,:);
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
    
    singleTimeVecX = linspace(round(respDataCell{1,2}(1,1)-respDataCell{1,2}(1,2)),round(respDataCell{1,2}(1,3)-respDataCell{1,2}(1,2)),length(motionVectorX));
    dispTimeThreshX(end+1) = motionVectorX((find(singleTimeVecX>timeThresh,1))) - motionVectorX((find(singleTimeVecX>0,1)));
    idxToThreshXSingle = find((motionVectorX - motionVectorX(find(singleTimeVecX>0,1))) > moveThresh & singleTimeVecX>0 & singleTimeVecX<=3,1);
    if ~isempty(idxToThreshXSingle)
        timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
    end
    singleTimeVecY = linspace(round(respDataCell{1,2}(1,1)-respDataCell{1,2}(1,2)),round(respDataCell{1,2}(1,3)-respDataCell{1,2}(1,2)),length(motionVectorY));
    dispTimeThreshY(end+1) =  (motionVectorY((find(singleTimeVecY>timeThresh,1))) - motionVectorY((find(singleTimeVecY>0,1))))*-1;
    idxToThreshYSingle = find((motionVectorY - motionVectorY(find(singleTimeVecY>0,1)))*-1 > moveThresh & singleTimeVecY>0 & singleTimeVecY<=3,1);
    if ~isempty(idxToThreshYSingle)
        timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
    end
end
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts_FS(motionEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts_FS(motionEventsLocationsY,90);
meanY = -1*meanY;
cIntFillPtsY = -1*cIntFillPtsY;
timeVecX = linspace(round((respDataCell{1,2}(1,1)-respDataCell{1,2}(1,2))*4)/4,round(respDataCell{1,2}(1,3)-respDataCell{1,2}(1,2)),length(meanX));
timeVecY = linspace(round((respDataCell{1,2}(1,1)-respDataCell{1,2}(1,2))*4)/4,round(respDataCell{1,2}(1,3)-respDataCell{1,2}(1,2)),length(meanY));

stdWindowSize = 5;
for n = stdWindowSize+1:length(meanY)
    if std(meanY(n-stdWindowSize:n)) > .01
        brainMotionStart = n;
        break
    end 
end

h(9) = figure('Color','White');
subplot(2,1,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
plot(timeVecX,meanX,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
set(f,'facea',[.2]);
plot([0 0],[-1 3],'r')
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[0,0,1,0.1])
end
plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
text(5,-1,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(5,5,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean Motion During Respiration Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Brian Shift (\mum)')
ylim([-1 3])
xlim([-.25 1])
grid on

subplot(2,1,2)
plot(timeVecY,meanY,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
set(f,'facea',[.2]);
plot([0 0],[-1 3],'r')
for n = 1:size(motionEventsLocationsY,1)
    plot(timeVecY,-1*motionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
end
plot(timeVecY(brainMotionStart),meanY(brainMotionStart),'rx')
hold off
text(5,-1,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(5,5,'Rostral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
xlabel('Time (s)')
ylabel('\Delta Brian Shift (\mum)')
ylim([-1 3])
xlim([-.25 1])
grid on
clear movementData

h(10) = figure('Color','White');
subplot(2,1,1)
histfit(timeToThreshX,numBins,'kernel')
title(['Time for brain to displace laterally ' num2str(moveThresh) ' micrometers following respiration'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 25])
hold on
mu = mean(timeToThreshX);
sig = std(timeToThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(timeToThreshY,numBins,'kernel')
title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following respiration'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 25])
hold on
mu = mean(timeToThreshY);
sig = std(timeToThreshY);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

h(50) = figure('Color','White');
subplot(2,1,1)
histfit(dispTimeThreshX,numBins,'kernel')
title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following respiration'])
xlabel('Displacement (\mum)')
xlim([-3 4])
ylim([0 40])
hold on
mu = mean(dispTimeThreshX);
sig = std(dispTimeThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(dispTimeThreshY,numBins,'kernel')
title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following respiration'])
xlabel('Time (s)')
xlim([-3 4])
ylim([0 40])
hold on
mu = mean(dispTimeThreshY);
sig = std(dispTimeThreshY);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotRespTriggeredAvgSkull_FS
motionEventsLocationsX = [];
motionEventsLocationsY = [];
timeToThreshX = [];
timeToThreshY = [];
dispTimeThreshX = [];
dispTimeThreshY = [];
moveThresh = .1;
timeThresh = .5;
numBins = 20;
load('respDataCellSkull_FS.mat')
for n = 1:size(respDataCell)
    if isnan(respDataCell{n,3})
        continue
    end
    motionVectorX = respDataCell{n,3}(2,:);
    motionVectorY = respDataCell{n,4}(2,:);
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
    
    singleTimeVecX = linspace(round(respDataCell{1,2}(1,1)-respDataCell{1,2}(1,2)),round(respDataCell{1,2}(1,3)-respDataCell{1,2}(1,2)),length(motionVectorX));
    dispTimeThreshX(end+1) = motionVectorX((find(singleTimeVecX>timeThresh,1))) - motionVectorX((find(singleTimeVecX>0,1)));
    idxToThreshXSingle = find(abs((motionVectorX - motionVectorX(find(singleTimeVecX>0,1)))) > moveThresh & singleTimeVecX>0 & singleTimeVecX<=3,1);
    if ~isempty(idxToThreshXSingle)
        timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
    end
    singleTimeVecY = linspace(round(respDataCell{1,2}(1,1)-respDataCell{1,2}(1,2)),round(respDataCell{1,2}(1,3)-respDataCell{1,2}(1,2)),length(motionVectorY));
    dispTimeThreshY(end+1) =  (motionVectorY((find(singleTimeVecY>timeThresh,1))) - motionVectorY((find(singleTimeVecY>0,1))))*-1;
    idxToThreshYSingle = find(abs((motionVectorY - motionVectorY(find(singleTimeVecY>0,1)))*-1) > moveThresh & singleTimeVecY>0 & singleTimeVecY<=3,1);
    if ~isempty(idxToThreshYSingle)
        timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
    end
end
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts_FS(motionEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts_FS(motionEventsLocationsY,90);
meanY = -1*meanY;
cIntFillPtsY = -1*cIntFillPtsY;
timeVecX = linspace(round((respDataCell{1,2}(1,1)-respDataCell{1,2}(1,2))*4)/4,round(respDataCell{1,2}(1,3)-respDataCell{1,2}(1,2)),length(meanX));
timeVecY = linspace(round((respDataCell{1,2}(1,1)-respDataCell{1,2}(1,2))*4)/4,round(respDataCell{1,2}(1,3)-respDataCell{1,2}(1,2)),length(meanY));

% stdWindowSize = 5;
% for n = stdWindowSize+1:length(meanY)
%     if std(meanY(n-stdWindowSize:n)) > .01
%         brainMotionStart = n;
%         break
%     end 
% end

h(9) = figure('Color','White');
subplot(2,1,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
plot(timeVecX,meanX,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
plot([0 0],[-1 3],'r')
set(f,'facea',[.2]);
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[0,0,1,0.1])
end
% plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
text(5,-1,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(5,5,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean Skull Motion During Respiration Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Skull Shift (\mum)')
ylim([-1 3])
xlim([-.25 1])
grid on

subplot(2,1,2)
plot(timeVecY,meanY,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
set(f,'facea',[.2]);
plot([0 0],[-1 3],'r')
for n = 1:size(motionEventsLocationsY,1)
    plot(timeVecY,-1*motionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
end
% plot(timeVecY(brainMotionStart),meanY(brainMotionStart),'rx')
hold off
text(5,-1,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(5,5,'Rostral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
xlabel('Time (s)')
ylabel('\Delta Skull Shift (\mum)')
ylim([-1 3])
xlim([-.25 1])
grid on
clear movementData

h(10) = figure('Color','White');
subplot(2,1,1)
histfit(timeToThreshX,numBins,'kernel')
title(['Time for skull to displace laterally ' num2str(moveThresh) ' micrometers following respiration'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 30])
hold on
mu = mean(timeToThreshX);
sig = std(timeToThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(timeToThreshY,numBins,'kernel')
title(['Time for skull to displace rostrally ' num2str(moveThresh) ' micrometers following respiration'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 30])
hold on
mu = mean(timeToThreshY);
sig = std(timeToThreshY);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

h(50) = figure('Color','White');
subplot(2,1,1)
histfit(dispTimeThreshX,numBins,'kernel')
title(['Lateral displacement of skull after ' num2str(timeThresh) ' s following respiration'])
xlabel('Displacement (\mum)')
xlim([-2 3])
ylim([0 60])
hold on
mu = mean(dispTimeThreshX);
sig = std(dispTimeThreshX);
plot([mu mu],[0 60],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(2,1,2)
histfit(dispTimeThreshY,numBins,'kernel')
title(['Rostral displacement of skull after ' num2str(timeThresh) ' s following respiration'])
xlabel('Displacement (\mum)')
xlim([-2 3])
ylim([0 60])
hold on
mu = mean(dispTimeThreshY);
sig = std(dispTimeThreshY);
plot([mu mu],[0 60],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotLocomotionTriggeredAvgEMG_FS
motionEventsLocationsX = [];
% motionEventsLocationsY = [];
timeToThreshX = [];
% timeToThreshY = [];
dispTimeThreshX = [];
% dispTimeThreshY = [];
moveThresh = .75;
timeThresh = 1.5;
numBins = 20;
load('LTADataCellEMG_FS.mat')
for n = 1:size(locDataCell)
    if isnan(locDataCell{n,3})
        continue
    end
    motionVectorX = locDataCell{n,3}(2,:);
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
    motionEventsLocationsX(end+1,:) = motionVectorX;
%     motionEventsLocationsY(end+1,:) = motionVectorY;
    
    singleTimeVecX = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(motionVectorX));
    dispTimeThreshX(end+1) = motionVectorX((find(singleTimeVecX>timeThresh,1))) - motionVectorX((find(singleTimeVecX>0,1)));
    idxToThreshXSingle = find((motionVectorX - motionVectorX(find(singleTimeVecX>0,1))) > moveThresh & singleTimeVecX>0 & singleTimeVecX<=3,1);
    if ~isempty(idxToThreshXSingle)
        timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
    end
%     singleTimeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(motionVectorY));
%     dispTimeThreshY(end+1) =  (motionVectorY((find(singleTimeVecY>timeThresh,1))) - motionVectorY((find(singleTimeVecY>0,1))))*-1;
%     idxToThreshYSingle = find((motionVectorY - motionVectorY(find(singleTimeVecY>0,1)))*-1 > moveThresh & singleTimeVecY>0 & singleTimeVecY<=3,1);
%     if ~isempty(idxToThreshYSingle)
%         timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
%     end
end
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts_FS(motionEventsLocationsX,90);
% [meanY,cIntFillPtsY] = getCIntMeanAndFillPts_FS(motionEventsLocationsY,90);
% meanY = -1*meanY;
% cIntFillPtsY = -1*cIntFillPtsY;
timeVecX = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(meanX));
% timeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(meanY));

% stdWindowSize = 5;
% for n = stdWindowSize+1:length(meanY)
%     if std(meanY(n-stdWindowSize:n)) > .01
%         brainMotionStart = n;
%         break
%     end 
% end

h(9) = figure('Color','White');
% subplot(2,1,1)
% maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
plot(timeVecX,meanX,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
set(f,'facea',[.2]);
plot([0 0],[0 4],'r')
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[0,0,1,0.1])
end
% plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
% text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean EMG During Locomotion Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('EKG Power (au)')
ylim([0 4])
xlim([-2 3])
grid on

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

stopMotionEventsLocationsX = [];
% stopMotionEventsLocationsY = [];
for n = 1:size(locDataCell)
    if isnan(locDataCell{n,5})
        continue
    end
    stopMotionVectorX = locDataCell{n,5}(2,:);
%     stopMotionVectorY = locDataCell{n,7}(2,:);
%     if n > 1
%         if length(stopMotionVectorX) > size(stopMotionEventsLocationsX,2)
%             stopMotionVectorX = stopMotionVectorX(1:size(stopMotionEventsLocationsX,2));
%         elseif length(stopMotionVectorX) < size(stopMotionEventsLocationsX,2)
%             stopMotionEventsLocationsX = stopMotionEventsLocationsX(:,1:length(stopMotionVectorX));
%         end
%         if length(stopMotionVectorY) > size(stopMotionEventsLocationsY,2)
%             stopMotionVectorY = stopMotionVectorY(1:size(stopMotionEventsLocationsY,2));
%         elseif length(stopMotionVectorY) < size(stopMotionEventsLocationsY,2)
%             stopMotionEventsLocationsY = stopMotionEventsLocationsY(:,1:length(stopMotionVectorY));
%         end
%     end
    stopMotionEventsLocationsX(end+1,:) = stopMotionVectorX;
%     stopMotionEventsLocationsY(end+1,:) = stopMotionVectorY;
end
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts_FS(stopMotionEventsLocationsX,90);
% [meanY,cIntFillPtsY] = getCIntMeanAndFillPts_FS(stopMotionEventsLocationsY,90);
% meanY = -1*meanY;
% cIntFillPtsY = -1*cIntFillPtsY;
timeVecX = linspace(round(locDataCell{2,4}(1,1)-locDataCell{2,4}(1,2)),round(locDataCell{2,4}(1,3)-locDataCell{2,4}(1,2)),length(meanX));
% timeVecY = linspace(round(locDataCell{1,5}(1,1)-locDataCell{1,5}(1,2)),round(locDataCell{1,5}(1,3)-locDataCell{1,5}(1,2)),length(meanY));

% subplot(2,1,2)
% % maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
% plot(timeVecX,meanX,'k')
% hold on
% f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
% set(f,'facea',[.2]);
% plot([0 0],[0 3],'r')
% for n = 1:size(stopMotionEventsLocationsX,1)
%     plot(timeVecX,stopMotionEventsLocationsX(n,:),'Color',[0,0,1,0.1])
% end
% hold off
% % text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% % text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% title(['\fontsize{20pt}\bf{Mean EMG During Stopping Locomotion Events, n = ' num2str(size(stopMotionEventsLocationsX,1)) '}'])
% xlabel('Time (s)')
% ylabel('EMG Power (au)')
% ylim([0 3])
% xlim([-2 3])
% grid on

% subplot(2,2,4)
% plot(timeVecY,meanY,'k')
% hold on
% f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
% set(f,'facea',[.2]);
% plot([0 0],[-3 3],'r')
% for n = 1:size(stopMotionEventsLocationsY,1)
%     plot(timeVecY,-1*stopMotionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
% end
% hold off
% text(3,-3,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(3,3,'Rostral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% xlabel('Time (s)')
% ylabel('\Delta Brian Shift (\mum)')
% ylim([-3 3])
% xlim([-2 3])
% grid on
% clear movementData

h(10) = figure('Color','White');
% subplot(2,1,1)
histfit(timeToThreshX,numBins,'kernel')
title(['Time for brain to displace laterally ' num2str(moveThresh) ' micrometers following locomotion trigger'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 25])
hold on
mu = mean(timeToThreshX);
sig = std(timeToThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

% subplot(2,1,2)
% histfit(timeToThreshY,numBins,'kernel')
% title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following locomotion trigger'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 25])
% hold on
% mu = mean(timeToThreshY);
% sig = std(timeToThreshY);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

h(50) = figure('Color','White');
% subplot(2,1,1)
histfit(dispTimeThreshX,numBins,'kernel')
title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following locomotion trigger'])
xlabel('Displacement (\mum)')
xlim([-3 4])
ylim([0 40])
hold on
mu = mean(dispTimeThreshX);
sig = std(dispTimeThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

% subplot(2,1,2)
% histfit(dispTimeThreshY,numBins,'kernel')
% title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following locomotion trigger'])
% xlabel('Time (s)')
% xlim([-3 4])
% ylim([0 40])
% hold on
% mu = mean(dispTimeThreshY);
% sig = std(dispTimeThreshY);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

% timeToThreshX = [];
% timeToThreshY = [];
% dispTimeThreshX = [];
% dispTimeThreshY = [];
% for n = 1:size(locDataCell)
%     if isnan(locDataCell{n,3})
%         continue
%     end
%     motionVectorX = locDataCell{n,3}(2,:);
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
%     motionEventsLocationsX(end+1,:) = motionVectorX;
%     motionEventsLocationsY(end+1,:) = motionVectorY;
%     
%     singleTimeVecX = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(motionVectorX));
%     timeThreshIdx = find(singleTimeVecX>(singleTimeVecX(brainMotionStart)+timeThresh),1);
%     dispTimeThreshX(end+1) = motionVectorX(timeThreshIdx) - motionVectorX(brainMotionStart);
%     idxToThreshXSingle = find((motionVectorX - motionVectorX(brainMotionStart)) > moveThresh & singleTimeVecX>singleTimeVecX(brainMotionStart) & singleTimeVecX<=3,1);
%     if ~isempty(idxToThreshXSingle)
%         timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
%     end
%     singleTimeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(motionVectorY));
%     dispTimeThreshY(end+1) = (motionVectorY(timeThreshIdx) - motionVectorY(brainMotionStart))*-1;
%     idxToThreshYSingle = find((motionVectorY - motionVectorY(brainMotionStart))*-1 > moveThresh & singleTimeVecY>singleTimeVecY(brainMotionStart) & singleTimeVecY<=3,1);
%     if ~isempty(idxToThreshYSingle)
%         timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
%     end
% end
% h(13) = figure('Color','White');
% subplot(2,1,1)
% histfit(timeToThreshX,numBins,'kernel')
% title(['Time for brain to displace laterally ' num2str(moveThresh) ' micrometers following brain motion start'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 25])
% hold on
% mu = mean(timeToThreshX);
% sig = std(timeToThreshX);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% subplot(2,1,2)
% histfit(timeToThreshY,numBins,'kernel')
% title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following brain motion start'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 25])
% hold on
% mu = mean(timeToThreshY);
% sig = std(timeToThreshY);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% h(50) = figure('Color','White');
% subplot(2,1,1)
% histfit(dispTimeThreshX,numBins,'kernel')
% title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following brain motion start'])
% xlabel('Displacement (\mum)')
% xlim([-3 4])
% ylim([0 40])
% hold on
% mu = mean(dispTimeThreshX);
% sig = std(dispTimeThreshX);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% subplot(2,1,2)
% histfit(dispTimeThreshY,numBins,'kernel')
% title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following brain motion start'])
% xlabel('Time (s)')
% xlim([-3 4])
% ylim([0 40])
% hold on
% mu = mean(dispTimeThreshY);
% sig = std(dispTimeThreshY);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotLocomotionTriggeredAvgEMGSingleTrial_FS(movementData)
load('LTADataCell_FS.mat')
    locIdx = find(strcmp(locDataCell(:,1),'D:/21-12-16_MouseExp/211216_002_processe_2layerBrainInSkullDataFinal.mat'));
%     targetPositionInSkull = combineBrainSkullMovement_FS(movementData,stationaryData);
% movementData.secondsPerFrame = movementData.secondsPerFrame/2;
%     posL2 = [movementData.emgData(:,1), movementData.emgData(:,2)];
%     meanPosValL2 = [.5*(posL2(1:end-1,1) + posL2(2:end,1)),.5*(posL2(1:end-1,2) + posL2(2:end,2))];
%     posL2 = zipperVecs(posL2,meanPosValL2);
    emgEventsLocationsX = [];
%     motionEventsLocationsY = [];
%     deleteRows = [];
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
        locTriggerEMGIdx = find(min(abs(locTriggerTime - movementData.emgData(:,1))) == abs(locTriggerTime - movementData.emgData(:,1)));
%         locDataCell{i,2}(n,4:6) = [locTriggerEMGIdx-60 locTriggerEMGIdx locTriggerEMGIdx+90];
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

motionEventsLocationsX = [];
% motionEventsLocationsY = [];
timeToThreshX = [];
% timeToThreshY = [];
dispTimeThreshX = [];
% dispTimeThreshY = [];
moveThresh = .75;
timeThresh = 1.5;
numBins = 20;
load('LTADataCellEMG_FS.mat')
for n = 1:size(emgEventsLocationsX)
%     if isnan(locDataCell{n,3})
%         continue
%     end
    motionVectorX = emgEventsLocationsX(n,:);
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
    motionEventsLocationsX(end+1,:) = motionVectorX;
%     motionEventsLocationsY(end+1,:) = motionVectorY;
    
    singleTimeVecX = linspace(-2,3,length(motionVectorX));
    dispTimeThreshX(end+1) = motionVectorX((find(singleTimeVecX>timeThresh,1))) - motionVectorX((find(singleTimeVecX>0,1)));
    idxToThreshXSingle = find((motionVectorX - motionVectorX(find(singleTimeVecX>0,1))) > moveThresh & singleTimeVecX>0 & singleTimeVecX<=3,1);
    if ~isempty(idxToThreshXSingle)
        timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
    end
%     singleTimeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(motionVectorY));
%     dispTimeThreshY(end+1) =  (motionVectorY((find(singleTimeVecY>timeThresh,1))) - motionVectorY((find(singleTimeVecY>0,1))))*-1;
%     idxToThreshYSingle = find((motionVectorY - motionVectorY(find(singleTimeVecY>0,1)))*-1 > moveThresh & singleTimeVecY>0 & singleTimeVecY<=3,1);
%     if ~isempty(idxToThreshYSingle)
%         timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
%     end
end
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts_FS(motionEventsLocationsX,90);
% [meanY,cIntFillPtsY] = getCIntMeanAndFillPts_FS(motionEventsLocationsY,90);
% meanY = -1*meanY;
% cIntFillPtsY = -1*cIntFillPtsY;
timeVecX = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(meanX));
% timeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(meanY));

% stdWindowSize = 5;
% for n = stdWindowSize+1:length(meanY)
%     if std(meanY(n-stdWindowSize:n)) > .01
%         brainMotionStart = n;
%         break
%     end 
% end

h(9) = figure('Color','White');
% subplot(2,1,1)
% maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
plot(timeVecX,meanX,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
set(f,'facea',[.2]);
plot([0 0],[0 4],'r')
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[0,0,1,0.1])
end
% plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
% text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean EMG During Locomotion Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('EKG Power (au)')
ylim([0 4])
xlim([-2 3])
grid on

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

stopMotionEventsLocationsX = [];
% stopMotionEventsLocationsY = [];
for n = 1:size(locDataCell)
    if isnan(locDataCell{n,5})
        continue
    end
    stopMotionVectorX = locDataCell{n,5}(2,:);
%     stopMotionVectorY = locDataCell{n,7}(2,:);
%     if n > 1
%         if length(stopMotionVectorX) > size(stopMotionEventsLocationsX,2)
%             stopMotionVectorX = stopMotionVectorX(1:size(stopMotionEventsLocationsX,2));
%         elseif length(stopMotionVectorX) < size(stopMotionEventsLocationsX,2)
%             stopMotionEventsLocationsX = stopMotionEventsLocationsX(:,1:length(stopMotionVectorX));
%         end
%         if length(stopMotionVectorY) > size(stopMotionEventsLocationsY,2)
%             stopMotionVectorY = stopMotionVectorY(1:size(stopMotionEventsLocationsY,2));
%         elseif length(stopMotionVectorY) < size(stopMotionEventsLocationsY,2)
%             stopMotionEventsLocationsY = stopMotionEventsLocationsY(:,1:length(stopMotionVectorY));
%         end
%     end
    stopMotionEventsLocationsX(end+1,:) = stopMotionVectorX;
%     stopMotionEventsLocationsY(end+1,:) = stopMotionVectorY;
end
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts_FS(stopMotionEventsLocationsX,90);
% [meanY,cIntFillPtsY] = getCIntMeanAndFillPts_FS(stopMotionEventsLocationsY,90);
% meanY = -1*meanY;
% cIntFillPtsY = -1*cIntFillPtsY;
timeVecX = linspace(round(locDataCell{2,4}(1,1)-locDataCell{2,4}(1,2)),round(locDataCell{2,4}(1,3)-locDataCell{2,4}(1,2)),length(meanX));
% timeVecY = linspace(round(locDataCell{1,5}(1,1)-locDataCell{1,5}(1,2)),round(locDataCell{1,5}(1,3)-locDataCell{1,5}(1,2)),length(meanY));

% subplot(2,1,2)
% % maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
% plot(timeVecX,meanX,'k')
% hold on
% f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
% set(f,'facea',[.2]);
% plot([0 0],[0 3],'r')
% for n = 1:size(stopMotionEventsLocationsX,1)
%     plot(timeVecX,stopMotionEventsLocationsX(n,:),'Color',[0,0,1,0.1])
% end
% hold off
% % text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% % text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% title(['\fontsize{20pt}\bf{Mean EMG During Stopping Locomotion Events, n = ' num2str(size(stopMotionEventsLocationsX,1)) '}'])
% xlabel('Time (s)')
% ylabel('EMG Power (au)')
% ylim([0 3])
% xlim([-2 3])
% grid on

% subplot(2,2,4)
% plot(timeVecY,meanY,'k')
% hold on
% f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
% set(f,'facea',[.2]);
% plot([0 0],[-3 3],'r')
% for n = 1:size(stopMotionEventsLocationsY,1)
%     plot(timeVecY,-1*stopMotionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
% end
% hold off
% text(3,-3,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(3,3,'Rostral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% xlabel('Time (s)')
% ylabel('\Delta Brian Shift (\mum)')
% ylim([-3 3])
% xlim([-2 3])
% grid on
% clear movementData

h(10) = figure('Color','White');
% subplot(2,1,1)
histfit(timeToThreshX,numBins,'kernel')
title(['Time for brain to displace laterally ' num2str(moveThresh) ' micrometers following locomotion trigger'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 25])
hold on
mu = mean(timeToThreshX);
sig = std(timeToThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

% subplot(2,1,2)
% histfit(timeToThreshY,numBins,'kernel')
% title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following locomotion trigger'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 25])
% hold on
% mu = mean(timeToThreshY);
% sig = std(timeToThreshY);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

h(50) = figure('Color','White');
% subplot(2,1,1)
histfit(dispTimeThreshX,numBins,'kernel')
title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following locomotion trigger'])
xlabel('Displacement (\mum)')
xlim([-3 4])
ylim([0 40])
hold on
mu = mean(dispTimeThreshX);
sig = std(dispTimeThreshX);
plot([mu mu],[0 40],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

% subplot(2,1,2)
% histfit(dispTimeThreshY,numBins,'kernel')
% title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following locomotion trigger'])
% xlabel('Time (s)')
% xlim([-3 4])
% ylim([0 40])
% hold on
% mu = mean(dispTimeThreshY);
% sig = std(dispTimeThreshY);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

% timeToThreshX = [];
% timeToThreshY = [];
% dispTimeThreshX = [];
% dispTimeThreshY = [];
% for n = 1:size(locDataCell)
%     if isnan(locDataCell{n,3})
%         continue
%     end
%     motionVectorX = locDataCell{n,3}(2,:);
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
%     motionEventsLocationsX(end+1,:) = motionVectorX;
%     motionEventsLocationsY(end+1,:) = motionVectorY;
%     
%     singleTimeVecX = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(motionVectorX));
%     timeThreshIdx = find(singleTimeVecX>(singleTimeVecX(brainMotionStart)+timeThresh),1);
%     dispTimeThreshX(end+1) = motionVectorX(timeThreshIdx) - motionVectorX(brainMotionStart);
%     idxToThreshXSingle = find((motionVectorX - motionVectorX(brainMotionStart)) > moveThresh & singleTimeVecX>singleTimeVecX(brainMotionStart) & singleTimeVecX<=3,1);
%     if ~isempty(idxToThreshXSingle)
%         timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
%     end
%     singleTimeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(motionVectorY));
%     dispTimeThreshY(end+1) = (motionVectorY(timeThreshIdx) - motionVectorY(brainMotionStart))*-1;
%     idxToThreshYSingle = find((motionVectorY - motionVectorY(brainMotionStart))*-1 > moveThresh & singleTimeVecY>singleTimeVecY(brainMotionStart) & singleTimeVecY<=3,1);
%     if ~isempty(idxToThreshYSingle)
%         timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
%     end
% end
% h(13) = figure('Color','White');
% subplot(2,1,1)
% histfit(timeToThreshX,numBins,'kernel')
% title(['Time for brain to displace laterally ' num2str(moveThresh) ' micrometers following brain motion start'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 25])
% hold on
% mu = mean(timeToThreshX);
% sig = std(timeToThreshX);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% subplot(2,1,2)
% histfit(timeToThreshY,numBins,'kernel')
% title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following brain motion start'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 25])
% hold on
% mu = mean(timeToThreshY);
% sig = std(timeToThreshY);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% h(50) = figure('Color','White');
% subplot(2,1,1)
% histfit(dispTimeThreshX,numBins,'kernel')
% title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following brain motion start'])
% xlabel('Displacement (\mum)')
% xlim([-3 4])
% ylim([0 40])
% hold on
% mu = mean(dispTimeThreshX);
% sig = std(dispTimeThreshX);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% subplot(2,1,2)
% histfit(dispTimeThreshY,numBins,'kernel')
% title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following brain motion start'])
% xlabel('Time (s)')
% xlim([-3 4])
% ylim([0 40])
% hold on
% mu = mean(dispTimeThreshY);
% sig = std(dispTimeThreshY);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotCalibration2D_FS
load('calibrationValues_FS.mat')
x = calibrationValues.file_210408_001;
widthData = [x.pixelMidpointXWidth x.pixelMidpointYWidth x.width x.leftBoundWidth x.rightBoundWidth x.polarityWidth];
heightData = [x.pixelMidpointXHeight x.pixelMidpointYHeight x.height x.topBoundHeight x.bottomBoundHeight x.polarityHeight];
barSize = 6;
holeSize = 19;
micronLengthValWidth = ones(size(widthData,1),1);
micronLengthValWidth(widthData(:,6) == 1,1) = barSize;
micronLengthValWidth(widthData(:,6) == -1,1) = holeSize;
micronLengthValHeight = ones(size(heightData,1),1);
micronLengthValHeight(heightData(:,6) == 1,1) = barSize;
micronLengthValHeight(heightData(:,6) == -1,1) = holeSize;
W = x.imageWidthPixels;
H = x.imageHeightPixels;

h(1) = figure('Color','White');
scatter(widthData(:,1),widthData(:,2),50,micronLengthValWidth./widthData(:,3),'filled');
hold on
rectangle('Position',[125 60 265 385])
hold off
colorbar;
k = colorbar;
ylabel(k,'\mum/pixel','FontSize',15)
hold off
axis equal
axis([1 512 1 512])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title(['\fontsize{20pt}\bf{Position of X Calibration Measurements}'])
xlabel('Pixels')
ylabel('Pixels')

h(2) = figure('Color','White');
scatter(heightData(:,1),heightData(:,2),50,micronLengthValHeight./heightData(:,3),'filled');
hold on
rectangle('Position',[125 60 265 385])
hold off
colorbar;
k = colorbar;
ylabel(k,'\mum/pixel','FontSize',15)
hold off
axis equal
axis([1 512 1 512])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title(['\fontsize{20pt}\bf{Position of Y Calibration Measurements}'])
xlabel('Pixels')
ylabel('Pixels')

% plot 3d surface plot of x calibration values across image
h(3) = figure('Color','White');
% micronLengthValWidth = ones(size(widthData,1),1);
% micronLengthValWidth(widthData(:,6) == 1,1) = barSize;
% micronLengthValWidth(widthData(:,6) == -1,1) = holeSize;
% surfaceCalibFitX = fit(widthData(:,1:2),micronLengthValWidth./widthData(:,3),'poly55');
% plot(x.surfaceCalibFitX,widthData(:,1:2),micronLengthValWidth./widthData(:,3));
plot(x.surfaceCalibFitX);
axis([1 W 1 H 0.5 ceil(max(micronLengthValWidth./widthData(:,3)))])
v = [-2 -5 0.8];
view(v);
title('\fontsize{20pt}\bf{Surface Calibration Fit (X)}')
xlabel('Pixels (X)')
ylabel('Pixels (Y)')
zlabel('\mum/Pixel')

% plot 3d surface plot of y calibration values across image
h(4) = figure('Color','White');
% micronLengthValHeight = ones(size(heightData,1),1);
% micronLengthValHeight(heightData(:,6) == 1,1) = barSize;
% micronLengthValHeight(heightData(:,6) == -1,1) = holeSize;
% surfaceCalibFitY = fit(heightData(:,1:2),micronLengthValHeight./heightData(:,3),'poly55');
% plot(x.surfaceCalibFitY,heightData(:,1:2),micronLengthValHeight./heightData(:,3));
plot(x.surfaceCalibFitY);
axis([1 W 1 H 0.5 ceil(max(micronLengthValHeight./heightData(:,3)))])
v = [-2 -5 0.8];
view(v);
title('\fontsize{20pt}\bf{Surface Calibration Fit (Y)}')
xlabel('Pixels (X)')
ylabel('Pixels (Y)')
zlabel('\mum/Pixel')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotCalibrationETL_FS
load('calibrationValuesFTL_FS.mat');
% inputZVec = [0.23 0.7];
fn = fieldnames(calibrationValues);
allDataX = [];
allDataY = [];
for n = 1:numel(fn)
    allDataX = [allDataX calibrationValues.(fn{n}).diopterVals];
    allDataY = [allDataY calibrationValues.(fn{n}).zMatchMicrons];
end
coefficients = polyfit(allDataX, allDataY, 3);
xFit = linspace(min(allDataX), max(allDataX), 1000);
yFit = polyval(coefficients , xFit);
yZeroDiff = 0;
yFit = yFit - yZeroDiff;
uniqueX = unique(allDataX);
for n = 1:length(uniqueX)
    i = allDataX == uniqueX(n);
    matchingY = allDataY(i);
    stdX(n) = uniqueX(n);
    stdYPlus(n) =  polyval(coefficients , uniqueX(n)) - yZeroDiff + std(matchingY);
    stdYMinus(n) = polyval(coefficients , uniqueX(n)) - yZeroDiff - std(matchingY);
end

h(1) = figure('Color','White');
subplot(2,1,1)
hold on
for n = 1:numel(fn)
    plot(calibrationValues.(fn{n}).diopterVals,calibrationValues.(fn{n}).zMatchMicrons - yZeroDiff,'-*');
end
f = fill([stdX, fliplr(stdX)], [stdYPlus, fliplr(stdYMinus)], 'r','Linestyle','none');
set(f,'facea',[.2]);
hold off
xlabel('Diopter Input (meters^{-1})')
ylabel('\Delta Z (microns)')
title('Controller Diopter Values vs. Focal Plane Position in Z')
grid on

subplot(2,1,2)
plot(0.23,0,'kx','MarkerSize',20)
hold on
plot(xFit, yFit, 'LineWidth', 2);
if exist('inputZVec','var')
    yFitInput = polyval(coefficients , inputZVec);
    for n = 1:length(yFitInput)
        plot(inputZVec(n),yFitInput(n),'k*','MarkerSize',15);
    end
    hold off
    text(.25,ceil(max(yFitInput))-.2,['Diopter Values Output: ' num2str(round(yFitInput,2))])
end
xlabel('Diopter Input (meters^{-1})')
ylabel('\Delta Z (microns)')
title('Controller Diopter Values vs. Focal Plane Position in Z')
grid on

load('calibrationValues_FS.mat');
diopterValsPower = -1.27:.3:1.73;
powerVals = [48 48 49 49 49 49 49 49 49 50 50]; % 16x, 30% power, mW
diopterVals = [.23 .53 .83 1.13 1.43 1.73 -.07 -0.37 -0.67 -0.97 -1.27];
widthVals = [];
heightVals = [];
borderDistance = 220;
lowBoundPixel = borderDistance;
highBoundPixel = 512-borderDistance;
for n = 1:11
    if n < 10
        fieldName = ['file_210408_00' num2str(n)];
    else
        fieldName = ['file_210408_0' num2str(n)];
    end
    widthLocX = calibrationValues.(fieldName).pixelMidpointXWidth;
    widthLocY = calibrationValues.(fieldName).pixelMidpointYWidth;
    widthInd = lowBoundPixel <= widthLocX & widthLocX <= highBoundPixel & lowBoundPixel <= widthLocY & widthLocY <= highBoundPixel;
    meanWidth = calibrationValues.(fieldName).negativePolarityLengthInMicrons/median(calibrationValues.(fieldName).width(widthInd));
    widthVals(end+1) = meanWidth;
    
    heightLocX = calibrationValues.(fieldName).pixelMidpointXHeight;
    heightLocY = calibrationValues.(fieldName).pixelMidpointYHeight;
    heightInd = lowBoundPixel <= heightLocX & heightLocX <= highBoundPixel & lowBoundPixel <= heightLocY & heightLocY <= highBoundPixel;
    meanHeight = calibrationValues.(fieldName).negativePolarityLengthInMicrons/median(calibrationValues.(fieldName).height(heightInd));
    heightVals(end+1) = meanHeight;
end
[diopterVals,i] = sort(diopterVals);
widthVals = widthVals(i);
heightVals = heightVals(i);

h(2) = figure('Color','White');
subplot(3,1,1)
plot(diopterVals,widthVals,'mo')
xlabel('Diopter Input (1/m)')
ylabel('\mum/Pixel in X')
xlim([-1.5 2])
ylim([.7 .9])
lsline
subplot(3,1,2)
plot(diopterVals,heightVals,'go')
xlabel('Diopter Input (1/m)')
ylabel('\mum/Pixel in Y')
xlim([-1.5 2])
ylim([.5 .7])
lsline
subplot(3,1,3)
plot(diopterValsPower,powerVals,'ko')
xlabel('Diopter Input (1/m)')
ylabel('Laser Power Through Objective (mW)')
xlim([-1.5 2])
ylim([46 52])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMotionTrackingBrainAndSkullResp_FS(movementData,stationaryData)
[movementData.targetPosition,stationaryData.targetPosition] = interpBrainSkullMovement_FS(movementData,stationaryData);
if movementData.hemisphere == 2
    movementData.targetPosition(:,1) = movementData.targetPosition(:,1)*-1;
    stationaryData.targetPosition(:,1) = stationaryData.targetPosition(:,1)*-1;
end
movementData.targetPosition(:,2) = movementData.targetPosition(:,2)*-1;
stationaryData.targetPosition(:,2) = stationaryData.targetPosition(:,2)*-1;
movementData.targetPosition(:,1) = movementData.targetPosition(:,1)-mean(movementData.targetPosition(1:50,1));
movementData.targetPosition(:,2) = movementData.targetPosition(:,2)-mean(movementData.targetPosition(1:50,2));
stationaryData.targetPosition(:,1) = stationaryData.targetPosition(:,1)-mean(stationaryData.targetPosition(1:50,1));
stationaryData.targetPosition(:,2) = stationaryData.targetPosition(:,2)-mean(stationaryData.targetPosition(1:50,2));
movementData.secondsPerFrame = movementData.secondsPerFrame/2;
h(6) = figure('Color','White');
subplot(4,1,1)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
hold off
title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
xlabel('Time (s)')
ylabel('X Position (\mum)')
grid on
axis([15 225 -3 3])
text(148,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(148,-3,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(4,1,2)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
hold off
xlabel('Time (s)')
ylabel('Y Position (\mum)')
grid on
axis([15 225 -3 3])
text(148,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(148,-3,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(4,1,3)
plot(movementData.emgData(:,1),movementData.emgData(:,2),'k')
xlabel('Time (s)')
ylabel('Abdominal EMG (au)')
grid on
axis([15 225 0.5 2])
subplot(4,1,4)
plot(movementData.videoRespiration(:,1)-1.7,movementData.videoRespiration(:,2),'b')
xlabel('Time (s)')
ylabel('Respiration (mean pixel intensity)')
grid on
axis([15 225 50 100])

h(6) = figure('Color','White');
subplot(4,1,1)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
hold off
title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
xlabel('Time (s)')
ylabel('X Position (\mum)')
grid on
axis([15 25 -3 3])
text(148,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(148,-3,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(4,1,2)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
hold off
xlabel('Time (s)')
ylabel('Y Position (\mum)')
grid on
axis([15 25 -3 3])
text(148,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(148,-3,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(4,1,3)
plot(movementData.emgData(:,1),movementData.emgData(:,2),'k')
xlabel('Time (s)')
ylabel('Abdominal EMG (au)')
grid on
axis([15 25 0.5 2])
subplot(4,1,4)
plot(movementData.videoRespiration(:,1)-1.7,movementData.videoRespiration(:,2),'b')
xlabel('Time (s)')
ylabel('Respiration (mean pixel intensity)')
grid on
axis([15 25 50 100])

h(6) = figure('Color','White');
subplot(4,1,1)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
hold off
title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
xlabel('Time (s)')
ylabel('X Position (\mum)')
grid on
axis([135 145 -3 3])
text(148,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(148,-3,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(4,1,2)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
hold off
xlabel('Time (s)')
ylabel('Y Position (\mum)')
grid on
axis([135 145 -3 3])
text(148,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(148,-3,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(4,1,3)
plot(movementData.emgData(:,1),movementData.emgData(:,2),'k')
xlabel('Time (s)')
ylabel('Abdominal EMG (au)')
grid on
axis([135 145 0.5 2])
subplot(4,1,4)
plot(movementData.videoRespiration(:,1)-1.7,movementData.videoRespiration(:,2),'b')
xlabel('Time (s)')
ylabel('Respiration (mean pixel intensity)')
grid on
axis([135 145 50 100])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMotionTrackingBrainAndSkullRespComp_FS(movementData,stationaryData)
[movementData.targetPosition,stationaryData.targetPosition] = interpBrainSkullMovement_FS(movementData,stationaryData);
if movementData.hemisphere == 2
    movementData.targetPosition(:,1) = movementData.targetPosition(:,1)*-1;
    stationaryData.targetPosition(:,1) = stationaryData.targetPosition(:,1)*-1;
end
movementData.targetPosition(:,2) = movementData.targetPosition(:,2)*-1;
stationaryData.targetPosition(:,2) = stationaryData.targetPosition(:,2)*-1;
movementData.targetPosition(:,1) = movementData.targetPosition(:,1)-mean(movementData.targetPosition(2480:2530,1));
movementData.targetPosition(:,2) = movementData.targetPosition(:,2)-mean(movementData.targetPosition(2480:2530,2));
stationaryData.targetPosition(:,1) = stationaryData.targetPosition(:,1)-mean(stationaryData.targetPosition(2480:2530,1));
stationaryData.targetPosition(:,2) = stationaryData.targetPosition(:,2)-mean(stationaryData.targetPosition(2480:2530,2));
movementData.secondsPerFrame = movementData.secondsPerFrame/2;
h(6) = figure('Color','White');
subplot(4,1,1)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
hold off
title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
xlabel('Time (s)')
ylabel('X Position (\mum)')
grid on
axis([62 272 -6 6])
text(62,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(62,-6,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(4,1,2)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
hold off
xlabel('Time (s)')
ylabel('Y Position (\mum)')
grid on
axis([62 272 -6 6])
text(62,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(62,-6,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(4,1,3)
plot(movementData.emgData(:,1),movementData.emgData(:,2),'k')
xlabel('Time (s)')
ylabel('Abdominal EMG (au)')
grid on
axis([62 272 0.5 3])
subplot(4,1,4)
plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)*2*pi*.06),'b')
xlabel('Time (s)')
ylabel('Locomotion (m/s)')
grid on
axis([62 272 0 0.2])

% h(6) = figure('Color','White');
% subplot(4,1,1)
% plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
% hold on
% plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
% hold off
% title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
% xlabel('Time (s)')
% ylabel('X Position (\mum)')
% grid on
% axis([15 25 -3 3])
% text(148,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(148,-3,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(4,1,2)
% plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
% hold on
% plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
% hold off
% xlabel('Time (s)')
% ylabel('Y Position (\mum)')
% grid on
% axis([15 25 -3 3])
% text(148,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(148,-3,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(4,1,3)
% plot(movementData.emgData(:,1),movementData.emgData(:,2),'k')
% xlabel('Time (s)')
% ylabel('Abdominal EMG (au)')
% grid on
% axis([15 25 0.5 2])
% subplot(4,1,4)
% plot(movementData.videoRespiration(:,1)-1.7,movementData.videoRespiration(:,2),'b')
% xlabel('Time (s)')
% ylabel('Respiration (mean pixel intensity)')
% grid on
% axis([15 25 50 100])
% 
% h(6) = figure('Color','White');
% subplot(4,1,1)
% plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
% hold on
% plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
% hold off
% title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
% xlabel('Time (s)')
% ylabel('X Position (\mum)')
% grid on
% axis([135 145 -3 3])
% text(148,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(148,-3,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(4,1,2)
% plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
% hold on
% plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
% hold off
% xlabel('Time (s)')
% ylabel('Y Position (\mum)')
% grid on
% axis([135 145 -3 3])
% text(148,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(148,-3,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(4,1,3)
% plot(movementData.emgData(:,1),movementData.emgData(:,2),'k')
% xlabel('Time (s)')
% ylabel('Abdominal EMG (au)')
% grid on
% axis([135 145 0.5 2])
% subplot(4,1,4)
% plot(movementData.videoRespiration(:,1)-1.7,movementData.videoRespiration(:,2),'b')
% xlabel('Time (s)')
% ylabel('Respiration (mean pixel intensity)')
% grid on
% axis([135 145 50 100])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hArrow = drawArrow_FS(p0,p1,color)
% drawArrow(p0,p1)
%
% Draws a simple arrow in 2D, from p0 to p1.
%
% INPUTS:
%   p0 = [x0; y0] = position of the tail
%   p1 = [x1; y1] = position of the tip
%   color = arrow color. Optional: default is black 
%       --> can be 'r','g','b','c','m','y','w', 'k' or a 1x3 color vector
%
% OUTPUTS:
%   hArrow = handle to the patch object representing the arrow
%

% Defaults:
if nargin == 2
   color = 'k'; 
end

% Parameters:
W1 = 0.08;   % half width of the arrow head, normalized by length of arrow
W2 = 0.014;  % half width of the arrow shaft
L1 = 0.18;   % Length of the arrow head, normalized by length of arrow
L2 = 0.13;  % Length of the arrow inset

% Unpack the tail and tip of the arrow
x0 = p0(1);
y0 = p0(2);
x1 = p1(1);
y1 = p1(2);

% Start by drawing an arrow from 0 to 1 on the x-axis
P = [...
    0, (1-L2), (1-L1), 1, (1-L1), (1-L2), 0;
    W2,    W2,     W1, 0,    -W1,    -W2, -W2];

% Scale,rotate, shift and plot:
dx = x1-x0;
dy = y1-y0;
Length = sqrt(dx*dx + dy*dy);
Angle = atan2(-dy,dx);
P = Length*P;   %Scale
P = [cos(Angle), sin(Angle); -sin(Angle), cos(Angle)]*P;  %Rotate
P = p0(:)*ones(1,7) + P;  %Shift

% Plot!
hArrow = patch(P(1,:), P(2,:),color);  axis equal;
hArrow.EdgeColor = color;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C,phi,S12,S1,S2,f,confC,phistd,Cerr] = coherencyc_FS(data1,data2,params)
%________________________________________________________________________________________________________________________
% Utilized in analysis by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Code unchanged with the exception of this title block for record keeping
%
%   Last Opened: February 23rd, 2019
%________________________________________________________________________________________________________________________
%
% Multi-taper coherency,cross-spectrum and individual spectra - continuous process
%
% Usage:
% [C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(data1,data2,params)
% Input: 
% Note units have to be consistent. See chronux.m for more information.
%       data1 (in form samples x trials) -- required
%       data2 (in form samples x trials) -- required
%       params: structure with fields tapers, pad, Fs, fpass, err, trialave
%       - optional
%           tapers : precalculated tapers from dpss or in the one of the following
%                    forms: 
%                    (1) A numeric vector [TW K] where TW is the
%                        time-bandwidth product and K is the number of
%                        tapers to be used (less than or equal to
%                        2TW-1). 
%                    (2) A numeric vector [W T p] where W is the
%                        bandwidth, T is the duration of the data and p 
%                        is an integer such that 2TW-p tapers are used. In
%                        this form there is no default i.e. to specify
%                        the bandwidth, you have to specify T and p as
%                        well. Note that the units of W and T have to be
%                        consistent: if W is in Hz, T must be in seconds
%                        and vice versa. Note that these units must also
%                        be consistent with the units of params.Fs: W can
%                        be in Hz if and only if params.Fs is in Hz.
%                        The default is to use form 1 with TW=3 and K=5
%
%	        pad		    (padding factor for the FFT) - optional (can take values -1,0,1,2...). 
%                    -1 corresponds to no padding, 0 corresponds to padding
%                    to the next highest power of 2 etc.
%			      	 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
%			      	 to 512 points, if pad=1, we pad to 1024 points etc.
%			      	 Defaults to 0.
%           Fs   (sampling frequency) - optional. Default 1.
%           fpass    (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%           err  (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
%                                   [0 p] or 0 - no error bars) - optional. Default 0.
%           trialave (average over trials when 1, don't average when 0) - optional. Default 0
% Output:
%       C (magnitude of coherency - frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       phi (phase of coherency - frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       S12 (cross spectrum -  frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       S1 (spectrum 1 - frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       S2 (spectrum 2 - frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       f (frequencies)
%       confC (confidence level for C at 1-p %) - only for err(1)>=1
%       phistd - theoretical/jackknife (depending on err(1)=1/err(1)=2) standard deviation for phi. 
%                Note that phi + 2 phistd and phi - 2 phistd will give 95% confidence
%                bands for phi - only for err(1)>=1 
%       Cerr  (Jackknife error bars for C - use only for Jackknife - err(1)=2)

if nargin < 2; error('Need data1 and data2'); end;
data1=change_row_to_column_FS(data1);
data2=change_row_to_column_FS(data2);
if nargin < 3; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave]=getparams_FS(params);
if nargout > 8 && err(1)~=2; 
    error('Cerr computed only for Jackknife. Correct inputs and run again');
end;
if nargout > 6 && err(1)==0;
%   Errors computed only if err(1) is nonzero. Need to change params and run again.
    error('When errors are desired, err(1) has to be non-zero.');
end;
N=check_consistency_FS(data1,data2);
nfft=max(2^(nextpow2(N)+pad),N);
[f,findx]=getfgrid_FS(Fs,nfft,fpass); 
tapers=dpsschk_FS(tapers,N,Fs); % check tapers
J1=mtfftc_FS(data1,tapers,nfft,Fs);
J2=mtfftc_FS(data2,tapers,nfft,Fs);
J1=J1(findx,:,:); J2=J2(findx,:,:);
S12=squeeze(mean(conj(J1).*J2,2));
S1=squeeze(mean(conj(J1).*J1,2));
S2=squeeze(mean(conj(J2).*J2,2));
if trialave; S12=squeeze(mean(S12,2)); S1=squeeze(mean(S1,2)); S2=squeeze(mean(S2,2)); end;
C12=S12./sqrt(S1.*S2);
C=abs(C12); 
phi=angle(C12);
if nargout>=9; 
     [confC,phistd,Cerr]=coherr_FS(C,J1,J2,err,trialave);
elseif nargout==8;
     [confC,phistd]=coherr_FS(C,J1,J2,err,trialave);
end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = change_row_to_column_FS(data)
%________________________________________________________________________________________________________________________
% Utilized in analysis by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Code unchanged with the exception of this title block for record keeping
%
%   Last Opened: February 23rd, 2019
%________________________________________________________________________________________________________________________
%
% Helper routine to transform 1d arrays into column vectors that are needed
% by other routines in Chronux
%
% Usage: data=change_row_to_column(data)
% 
% Inputs:
% data -- required. If data is a matrix, it is assumed that it is of the
% form samples x channels/trials and it is returned without change. If it
% is a vector, it is transformed to a column vector. If it is a struct
% array of dimension 1, it is again returned as a column vector. If it is a
% struct array with multiple dimensions, it is returned without change
% Note that the routine only looks at the first field of a struct array.
% 
% Ouputs:
% data (in the form samples x channels/trials)
%

dtmp=[];
if isstruct(data);
   C=length(data);
   if C==1;
      fnames=fieldnames(data);
      eval(['dtmp=data.' fnames{1} ';'])
      data=dtmp(:);
   end
else
  [N,C]=size(data);
  if N==1 || C==1;
    data=data(:);
  end;
end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tapers,pad,Fs,fpass,err,trialave,params] = getparams_FS(params)
%________________________________________________________________________________________________________________________
% Utilized in analysis by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Code unchanged with the exception of this title block for record keeping
%
%   Last Opened: February 23rd, 2019
%________________________________________________________________________________________________________________________
%
% Helper function to convert structure params to variables used by the
% various routines - also performs checks to ensure that parameters are
% defined; returns default values if they are not defined.
%
% Usage: [tapers,pad,Fs,fpass,err,trialave,params]=getparams(params)
%
% Inputs:
%       params: structure with fields tapers, pad, Fs, fpass, err, trialave
%           - optional
%             tapers : precalculated tapers from dpss or in the one of the following
%                       forms:  
%                       (1) A numeric vector [TW K] where TW is the
%                           time-bandwidth product and K is the number of
%                           tapers to be used (less than or equal to
%                           2TW-1). 
%                       (2) A numeric vector [W T p] where W is the
%                           bandwidth, T is the duration of the data and p 
%                           is an integer such that 2TW-p tapers are used. In
%                           this form there is no default i.e. to specify
%                           the bandwidth, you have to specify T and p as
%                           well. Note that the units of W and T have to be
%			                consistent: if W is in Hz, T must be in seconds
% 			                and vice versa. Note that these units must also
%			                be consistent with the units of params.Fs: W can
%		    	            be in Hz if and only if params.Fs is in Hz.
%                           The default is to use form 1 with TW=3 and K=5
%
%	        pad		    (padding factor for the FFT) - optional (can take values -1,0,1,2...). 
%                    -1 corresponds to no padding, 0 corresponds to padding
%                    to the next highest power of 2 etc.
%			      	 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
%			      	 to 512 points, if pad=1, we pad to 1024 points etc.
%			      	 Defaults to 0.
%           Fs   (sampling frequency) - optional. Default 1.
%           fpass    (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%           err  (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
%                                   [0 p] or 0 - no error bars) - optional. Default 0.
%           trialave (average over trials when 1, don't average when 0) - optional. Default 0
% Outputs: 
% The fields listed above as well as the struct params. The fields are used
% by some routines and the struct is used by others. Though returning both
% involves overhead, it is a safer, simpler thing to do.

if ~isfield(params,'tapers') || isempty(params.tapers);  %If the tapers don't exist
     display('tapers unspecified, defaulting to params.tapers=[3 5]');
     params.tapers=[3 5];
end;
if ~isempty(params) && length(params.tapers)==3 
    % Compute timebandwidth product
    TW = params.tapers(2)*params.tapers(1);
    % Compute number of tapers
    K  = floor(2*TW - params.tapers(3));
    params.tapers = [TW  K];
end

if ~isfield(params,'pad') || isempty(params.pad);
    params.pad=0;
end;
if ~isfield(params,'Fs') || isempty(params.Fs);
    params.Fs=1;
end;
if ~isfield(params,'fpass') || isempty(params.fpass);
    params.fpass=[0 params.Fs/2];
end;
if ~isfield(params,'err') || isempty(params.err);
    params.err=0;
end;
if ~isfield(params,'trialave') || isempty(params.trialave);
    params.trialave=0;
end;

tapers=params.tapers;
pad=params.pad;
Fs=params.Fs;
fpass=params.fpass;
err=params.err;
trialave=params.trialave;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N,C] = check_consistency_FS(data1,data2,sp)
%________________________________________________________________________________________________________________________
% Utilized in analysis by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Code unchanged with the exception of this title block for record keeping
%
%   Last Opened: February 23rd, 2019
%________________________________________________________________________________________________________________________
%
% Helper routine to check consistency of data dimensions
% Usage: [N,C]=check_consistency(data1,data2,sp)
% Inputs:
% data1 - first dataset
% data2 - second dataset
% sp - optional argument to be input as 1 when one of the two data sets is
% spikes times stored as a 1d array.
% Outputs:
% Dimensions of the datasets - data1 or data2 (note that 
%    routine stops with an error message if dimensions don't match - [N,C]
%    N left empty for structure arrays
N1=[]; N2=[];
if nargin < 3 || isempty(sp); sp=0; end;
if isstruct(data1);
    C1=length(data1);
else
    [N1,C1]=size(data1);
end;
if isstruct(data2);
    C2=length(data2);
else
    [N2,C2]=size(data2);
end;
if C1~=C2; error('inconsistent dimensions'); end;
if sp==0;
   if ~isstruct(data1) && ~isstruct(data2);
      if N1~=N2; error('inconsistent dimensions'); end;
   end;
end;
N=N1; C=C1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f,findx] = getfgrid_FS(Fs,nfft,fpass)
%________________________________________________________________________________________________________________________
% Utilized in analysis by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Code unchanged with the exception of this title block for record keeping
%
%   Last Opened: February 23rd, 2019
%________________________________________________________________________________________________________________________
%
% Helper function that gets the frequency grid associated with a given fft based computation
% Called by spectral estimation routines to generate the frequency axes 
% Usage: [f,findx]=getfgrid(Fs,nfft,fpass)
% Inputs:
% Fs        (sampling frequency associated with the data)-required
% nfft      (number of points in fft)-required
% fpass     (band of frequencies at which the fft is being calculated [fmin fmax] in Hz)-required
% Outputs:
% f         (frequencies)
% findx     (index of the frequencies in the full frequency grid). e.g.: If
% Fs=1000, and nfft=1048, an fft calculation generates 512 frequencies
% between 0 and 500 (i.e. Fs/2) Hz. Now if fpass=[0 100], findx will
% contain the indices in the frequency grid corresponding to frequencies <
% 100 Hz. In the case fpass=[0 500], findx=[1 512].
if nargin < 3; error('Need all arguments'); end;
df=Fs/nfft;
f=0:df:Fs; % all possible frequencies
f=f(1:nfft);
if length(fpass)~=1;
   findx=find(f>=fpass(1) & f<=fpass(end));
else
   [fmin,findx]=min(abs(f-fpass));
   clear fmin
end;
f=f(findx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tapers,eigs] = dpsschk_FS(tapers,N,Fs)
%________________________________________________________________________________________________________________________
% Utilized in analysis by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Code unchanged with the exception of this title block for record keeping
%
%   Last Opened: February 23rd, 2019
%________________________________________________________________________________________________________________________
%
% Helper function to calculate tapers and, if precalculated tapers are supplied, 
% to check that they (the precalculated tapers) the same length in time as
% the time series being studied. The length of the time series is specified
% as the second input argument N. Thus if precalculated tapers have
% dimensions [N1 K], we require that N1=N.
% Usage: tapers=dpsschk(tapers,N,Fs)
% Inputs:
% tapers        (tapers in the form of: 
%                                   (i) precalculated tapers or,
%                                   (ii) [NW K] - time-bandwidth product, number of tapers) 
%
% N             (number of samples)
% Fs            (sampling frequency - this is required for nomalization of
%                                     tapers: we need tapers to be such
%                                     that integral of the square of each taper equals 1
%                                     dpss computes tapers such that the
%                                     SUM of squares equals 1 - so we need
%                                     to multiply the dpss computed tapers
%                                     by sqrt(Fs) to get the right
%                                     normalization)
% Outputs: 
% tapers        (calculated or precalculated tapers)
% eigs          (eigenvalues) 
if nargin < 3; error('Need all arguments'); end
sz=size(tapers);
if sz(1)==1 && sz(2)==2;
    [tapers,eigs]=dpss(N,tapers(1),tapers(2));
    tapers = tapers*sqrt(Fs);
elseif N~=sz(1);
    error('seems to be an error in your dpss calculation; the number of time points is different from the length of the tapers');
end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J = mtfftc_FS(data,tapers,nfft,Fs)
%________________________________________________________________________________________________________________________
% Utilized in analysis by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Code unchanged with the exception of this title block for record keeping
%
%   Last Opened: February 23rd, 2019
%________________________________________________________________________________________________________________________
%
% Multi-taper fourier transform - continuous data
%
% Usage:
% J=mtfftc(data,tapers,nfft,Fs) - all arguments required
% Input: 
%       data (in form samples x channels/trials or a single vector) 
%       tapers (precalculated tapers from dpss) 
%       nfft (length of padded data)
%       Fs   (sampling frequency)
%                                   
% Output:
%       J (fft in form frequency index x taper index x channels/trials)
if nargin < 4; error('Need all input arguments'); end;
data=change_row_to_column_FS(data);
[NC,C]=size(data); % size of data
[NK K]=size(tapers); % size of tapers
if NK~=NC; error('length of tapers is incompatible with length of data'); end;
tapers=tapers(:,:,ones(1,C)); % add channel indices to tapers
data=data(:,:,ones(1,K)); % add taper indices to data
data=permute(data,[1 3 2]); % reshape data to get dimensions to match those of tapers
data_proj=data.*tapers; % product of data with tapers
J=fft(data_proj,nfft)/Fs;   % fft of projected data
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [confC,phistd,Cerr] = coherr_FS(C,J1,J2,err,trialave,numsp1,numsp2)
%________________________________________________________________________________________________________________________
% Utilized in analysis by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Code unchanged with the exception of this title block for record keeping
%
%   Last Opened: February 23rd, 2019
%________________________________________________________________________________________________________________________
%
% Function to compute lower and upper confidence intervals on the coherency 
% given the tapered fourier transforms, errchk, trialave.
%
% Usage: [confC,phistd,Cerr]=coherr(C,J1,J2,err,trialave,numsp1,numsp2)
% Inputs:
% C     - coherence
% J1,J2 - tapered fourier transforms 
% err - [errtype p] (errtype=1 - asymptotic estimates; errchk=2 - Jackknife estimates; 
%                   p - p value for error estimates)
% trialave - 0: no averaging over trials/channels
%            1 : perform trial averaging
% numsp1    - number of spikes for data1. supply only if finite size corrections are required
% numsp2    - number of spikes for data2. supply only if finite size corrections are required
%
% Outputs: 
%          confC - confidence level for C - only for err(1)>=1
%          phistd - theoretical or jackknife standard deviation for phi for err(1)=1 and err(1)=2 
%                   respectively. returns zero if coherence is 1
%          Cerr  - Jacknife error bars for C  - only for err(1)=2

if nargin < 5; error('Need at least 5 input arguments'); end;
if err(1)==0; error('Need err=[1 p] or [2 p] for error bar calculation'); end;
if nargout==4  && err(1)==1; error('Cerr contains Jackknife errors: only computed when err(1) is 2'); end;
[nf,K,Ch]=size(J1);
errchk=err(1);
p=err(2);
pp=1-p/2;
%
% Find the number of degrees of freedom
%
if trialave;
   dim=K*Ch;
   dof=2*dim;
   dof1=dof;
   dof2=dof;
   Ch=1;
   if nargin>=6 && ~isempty(numsp1) 
      totspikes1=sum(numsp1);
      dof1=fix(2*totspikes1*dof/(2*totspikes1+dof));
   end
   if nargin==7 && ~isempty(numsp2); 
      totspikes2=sum(numsp2);
      dof2=fix(2*totspikes2*dof/(2*totspikes2+dof));
   end;
   dof=min(dof1,dof2);
   J1=reshape(J1,nf,dim);
   J2=reshape(J2,nf,dim);
else
   dim=K;
   dof=2*dim;
   dof1=dof;
   dof2=dof;
   for ch=1:Ch;
      if nargin>=6 && ~isempty(numsp1);
         totspikes1=numsp1(ch); 
        dof1=fix(2*totspikes1*dof/(2*totspikes1+dof));
      end;
      if nargin==7 && ~isempty(numsp2);
         totspikes2=numsp2(ch);
        dof2=fix(2*totspikes2*dof/(2*totspikes2+dof));
      end;
      dof(ch)=min(dof1,dof2);
   end;
end;
%
% variance of the phase
%
%
% Old code is the next few lines - new code is in the if statement below
% beginning line 87
%
% if isempty(find((C-1).^2 < 10^-16));
%    phierr = sqrt((2./dof(ones(nf,1),:)).*(1./(C.^2) - 1));  
% else
%    phierr = zeros(nf,Ch);
% end  

%
% theoretical, asymptotic confidence level
%
if dof <= 2
   confC = 1;
else     
   df = 1./((dof/2)-1);
   confC = sqrt(1 - p.^df);
end;
%
% Phase standard deviation (theoretical and jackknife) and jackknife
% confidence intervals for C
%
if errchk==1;
   totnum=nf*Ch;
   phistd=zeros(totnum,1); 
   CC=reshape(C,[totnum,1]); 
   indx=find(abs(CC-1)>=1.e-16);
   dof=repmat(dof,[nf,1]);
   dof=reshape(dof,[totnum 1]); 
   phistd(indx)= sqrt((2./dof(indx).*(1./(C(indx).^2) - 1))); 
   phistd=reshape(phistd,[nf Ch]);
elseif errchk==2;
    tcrit=tinv(pp,dof-1);
    for k=1:dim;
        indxk=setdiff(1:dim,k);
        J1k=J1(:,indxk,:);
        J2k=J2(:,indxk,:);
        eJ1k=squeeze(sum(J1k.*conj(J1k),2));
        eJ2k=squeeze(sum(J2k.*conj(J2k),2));
        eJ12k=squeeze(sum(conj(J1k).*J2k,2)); 
        Cxyk=eJ12k./sqrt(eJ1k.*eJ2k);
        absCxyk=abs(Cxyk);
        atanhCxyk(k,:,:)=sqrt(2*dim-2)*atanh(absCxyk);
        phasefactorxyk(k,:,:)=Cxyk./absCxyk;
%         indxk=setdiff(1:dim,k);
%         J1jk=J1(:,indxk,:);
%         J2jk=J2(:,indxk,:);
%         eJ1jk=squeeze(sum(J1jk.*conj(J1jk),2));
%         eJ2jk=squeeze(sum(J2jk.*conj(J2jk),2));
%         eJ12jk=squeeze(sum(conj(J1jk).*J2jk,2)); 
%         atanhCxyjk(k,:,:)=sqrt(2*dim-2)*atanh(abs(eJ12jk)./sqrt(eJ1jk.*eJ2jk));
    end; 
    atanhC=sqrt(2*dim-2)*atanh(C);
    sigma12=sqrt(dim-1)*squeeze(std(atanhCxyk,1,1));
%     sigma12=sqrt(dim-1)*squeeze(std(atanhCxyjk,1,1));
    if Ch==1; sigma12=sigma12'; end;
    Cu=atanhC+tcrit(ones(nf,1),:).*sigma12;
    Cl=atanhC-tcrit(ones(nf,1),:).*sigma12;
    Cerr(1,:,:) = tanh(Cl/sqrt(2*dim-2));
    Cerr(2,:,:) = tanh(Cu/sqrt(2*dim-2));
    phistd=(2*dim-2)*(1-abs(squeeze(mean(phasefactorxyk))));
    if trialave; phistd=phistd'; end;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S,f,Serr]=mtspectrumc_FS(data,params)
% Multi-taper spectrum - continuous process
%
% Usage:
%
% [S,f,Serr]=mtspectrumc(data,params)
% Input: 
% Note units have to be consistent. See chronux.m for more information.
%       data (in form samples x channels/trials) -- required
%       params: structure with fields tapers, pad, Fs, fpass, err, trialave
%       -optional
%           tapers : precalculated tapers from dpss or in the one of the following
%                    forms: 
%                    (1) A numeric vector [TW K] where TW is the
%                        time-bandwidth product and K is the number of
%                        tapers to be used (less than or equal to
%                        2TW-1). 
%                    (2) A numeric vector [W T p] where W is the
%                        bandwidth, T is the duration of the data and p 
%                        is an integer such that 2TW-p tapers are used. In
%                        this form there is no default i.e. to specify
%                        the bandwidth, you have to specify T and p as
%                        well. Note that the units of W and T have to be
%                        consistent: if W is in Hz, T must be in seconds
%                        and vice versa. Note that these units must also
%                        be consistent with the units of params.Fs: W can
%                        be in Hz if and only if params.Fs is in Hz.
%                        The default is to use form 1 with TW=3 and K=5
%
%	        pad		    (padding factor for the FFT) - optional (can take values -1,0,1,2...). 
%                    -1 corresponds to no padding, 0 corresponds to padding
%                    to the next highest power of 2 etc.
%			      	 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
%			      	 to 512 points, if pad=1, we pad to 1024 points etc.
%			      	 Defaults to 0.
%           Fs   (sampling frequency) - optional. Default 1.
%           fpass    (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%           err  (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
%                                   [0 p] or 0 - no error bars) - optional. Default 0.
%           trialave (average over trials/channels when 1, don't average when 0) - optional. Default 0
% Output:
%       S       (spectrum in form frequency x channels/trials if trialave=0; 
%               in the form frequency if trialave=1)
%       f       (frequencies)
%       Serr    (error bars) only for err(1)>=1

if nargin < 1; error('Need data'); end;
if nargin < 2; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams_FS(params);
if nargout > 2 && err(1)==0; 
%   Cannot compute error bars with err(1)=0. Change params and run again. 
    error('When Serr is desired, err(1) has to be non-zero.');
end;
data=change_row_to_column_FS(data);
N=size(data,1);
nfft=max(2^(nextpow2(N)+pad),N);
[f,findx]=getfgrid_FS(Fs,nfft,fpass); 
tapers=dpsschk_FS(tapers,N,Fs); % check tapers
J=mtfftc_FS(data,tapers,nfft,Fs);
J=J(findx,:,:);
S=permute(mean(conj(J).*J,2),[1 3 2]);
if trialave; S=squeeze(mean(S,2));else S=squeeze(S);end;
if nargout==3; 
   Serr=specerr_FS(S,J,err,trialave);
end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Serr=specerr_FS(S,J,err,trialave,numsp)
% Function to compute lower and upper confidence intervals on the spectrum 
% Usage: Serr=specerr(S,J,err,trialave,numsp)
% Outputs: Serr (Serr(1,...) - lower confidence level, Serr(2,...) upper confidence level)
%
% Inputs:
% S - spectrum
% J - tapered fourier transforms 
% err - [errtype p] (errtype=1 - asymptotic estimates; errchk=2 - Jackknife estimates; 
%                   p - p value for error estimates)
% trialave - 0: no averaging over trials/channels
%            1 : perform trial averaging
% numsp    - number of spikes in each channel. specify only when finite
%            size correction required (and of course, only for point
%            process data)
%
% Outputs:
% Serr - error estimates. Only for err(1)>=1. If err=[1 p] or [2 p] Serr(...,1) and Serr(...,2)
% contain the lower and upper error bars with the specified method. 
if nargin < 4; error('Need at least 4 input arguments'); end;
if err(1)==0; error('Need err=[1 p] or [2 p] for error bar calculation. Make sure you are not asking for the output of Serr'); end;
[nf,K,C]=size(J);
errchk=err(1);
p=err(2);
pp=1-p/2;
qq=1-pp;

if trialave
   dim=K*C;
   C=1;
   dof=2*dim;
   if nargin==5; dof = fix(1/(1/dof + 1/(2*sum(numsp)))); end
   J=reshape(J,nf,dim);
else
   dim=K;
   dof=2*dim*ones(1,C);
   for ch=1:C;
     if nargin==5; dof(ch) = fix(1/(1/dof + 1/(2*numsp(ch)))); end 
   end;
end;
Serr=zeros(2,nf,C);
if errchk==1;
   Qp=chi2inv(pp,dof);
   Qq=chi2inv(qq,dof);
   Serr(1,:,:)=dof(ones(nf,1),:).*S./Qp(ones(nf,1),:);
   Serr(2,:,:)=dof(ones(nf,1),:).*S./Qq(ones(nf,1),:);
elseif errchk==2;
   tcrit=tinv(pp,dim-1);
   for k=1:dim;
       indices=setdiff(1:dim,k);
       Jjk=J(:,indices,:); % 1-drop projection
       eJjk=squeeze(sum(Jjk.*conj(Jjk),2));
       Sjk(k,:,:)=eJjk/(dim-1); % 1-drop spectrum
   end;
   sigma=sqrt(dim-1)*squeeze(std(log(Sjk),1,1)); if C==1; sigma=sigma'; end; 
   conf=repmat(tcrit,nf,C).*sigma;
   conf=squeeze(conf); 
   Serr(1,:,:)=S.*exp(-conf); Serr(2,:,:)=S.*exp(conf);
end;
Serr=squeeze(Serr);
end