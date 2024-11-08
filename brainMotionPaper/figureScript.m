%% Script to generate all figures for paper

clear;
close all;

%% Process data sets

combinedMovementDataBrain_221207_008 = combineMotionTracking_FS('221207_008_Layer1_',1:3);
combinedMovementDataSkull_221207_008 = combineMotionTracking_FS('221207_008_Layer2_',1:3);

combinedMovementDataBrain_221207_010 = combineMotionTracking_FS('221207_010_Layer1_',1:3);
combinedMovementDataSkull_221207_010 = combineMotionTracking_FS('221207_010_Layer2_',1:3);

combinedMovementDataBrain_221208_006 = combineMotionTracking_FS('221208_006_Layer1_',1:3);
combinedMovementDataSkull_221208_006 = combineMotionTracking_FS('221208_006_Layer2_',1:3);

combinedMovementData_210218_006 = combineMotionTracking_FS('210218_006_',1:8);

combinedMovementDataBrain_211216_001 = combineMotionTracking_FS('211216_001_Layer1_',1:4);
combinedMovementDataSkull_211216_001 = combineMotionTracking_FS('211216_001_Layer2_',1:4);

combinedMovementDataBrain_211216_002 = combineMotionTracking_FS('211216_002_Layer1_',1:5);
combinedMovementDataSkull_211216_002 = combineMotionTracking_FS('211216_002_Layer2_',1:5);

combinedMovementDataBrain_220816_004 = combineMotionTracking_FS('220816_004_Layer1_',1:3);
combinedMovementDataSkull_220816_004 = combineMotionTracking_FS('220816_004_Layer2_',1:3);

combinedMovementDataBrain_240611_001 = combineMotionTracking_FS('240611_001_Layer1_',1:4);
combinedMovementDataSkull_240611_001 = combineMotionTracking_FS('240611_001_Layer2_',1:4);

combinedMovementDataBrain_240612_002 = combineMotionTracking_FS('240612_002_Layer1_',1:3);
combinedMovementDataSkull_240612_002 = combineMotionTracking_FS('240612_002_Layer2_',1:3);

combinedMovementDataBrain_240612_020 = combineMotionTracking_FS('240612_020_Layer1_',1:3);
combinedMovementDataSkull_240612_020 = combineMotionTracking_FS('240612_020_Layer2_',1:3);

combinedMovementDataBrain_220813_003 = combineMotionTracking_FS('220813_003_Layer1_',1:3);
combinedMovementDataSkull_220813_003 = combineMotionTracking_FS('220813_003_Layer2_',1:3);

combinedMovementDataSkull_240814_006 = combineMotionTracking_FS('240814_006_',1:3);

combinedMovementDataBrain_221205_004 = combineMotionTracking_FS('221205_004_Layer1_',1:3);
combinedMovementDataDura_221205_004 = combineMotionTracking_FS('221205_004_Layer2_',1:3);
combinedMovementDataSkull_221205_004 = combineMotionTracking_FS('221205_004_Layer3_',1:3);

load('211216_002_rawData.mat');

%% figure 1e

plotMotionTrackingBrainAndSkull_FS(combinedMovementDataBrain_221207_010,combinedMovementDataSkull_221207_010);

%% figure 2a

plotMotionTracking2P2LPCA_FS(combinedMovementDataBrain_211216_002,combinedMovementDataSkull_211216_002);

%% figure 2b

plotMovementQuiver_FS

%% figure 2c

plotThermoCoherence_FS(combinedMovementDataBrain_211216_002,combinedMovementDataSkull_211216_002,rawData)

%% figure 2d

plotLocomotionXCorr_FS(combinedMovementDataBrain_211216_002,combinedMovementDataSkull_211216_002)

%% figure 2e

[LTATimeData,LTADispData] = plotLocomotionTriggeredAvg_FS;

%% figure 3b

plotLocomotionTriggeredAvgEMGSingleTrial_FS(combinedMovementDataBrain_211216_002)

%% figure 3c

plotLocomotionTriggeredAvgEMG_FS

%% figure 3d

plotMovementBrainInSkull_FS(combinedMovementDataBrain_211216_002,combinedMovementDataSkull_211216_002,rawData)

%% figure 3e

plotEMGMovement_FS(combinedMovementDataBrain_211216_002,combinedMovementDataSkull_211216_002)

%% figure 3f

plotEMGXCorr_FS(combinedMovementDataBrain_211216_002,combinedMovementDataSkull_211216_002)

%% figure 3g

[ETATimeData,ETADispData] = plotEMGTriggeredAvg_FS;

%% figure 5b

plotMotionTracking2P2LPCASqueeze_FS(combinedMovementDataBrain_240612_020,combinedMovementDataSkull_240612_020,'Figure 5b');

%% figure 5c

plotMotionTrackingBrainAndSkullSqueeze_FS(combinedMovementDataBrain_240612_020,combinedMovementDataSkull_240612_020,'Figure 5c');

%% figure 5d

plotSqueezeQuiver_FS

%% figure 5e

[STATimeData,STADispData] = plotSqueezeTriggeredAvg_FS;

%% figure 5f

plotSqueezeTriggeredAvgSkull_FS

%% supplementary figure 2b-c

plotCalibration2D_FS

%% supplementary figure 3c

plotPSF_FS

%% supplementary figure 3f-g

plotCalibrationETL_FS

%% supplementary figure 4b

[RTATimeData,RTADispData] = plotRespTriggeredAvg_FS;

%% supplementary figure 4c

plotRespTriggeredAvgSkull_FS

%% supplementary figure 4d

plotMotionTrackingBrainAndSkullResp_FS(combinedMovementDataBrain_240611_001,combinedMovementDataSkull_240611_001,'Supplementary Figure 4d_1');

%% supplementary figure 4e

plotMotionTrackingBrainAndSkullRespComp_FS(combinedMovementDataBrain_240612_002,combinedMovementDataSkull_240612_002);

%% supplementary figure 6a

plotMotionTrackingSkullBregma_FS(combinedMovementDataSkull_240814_006);

%% supplementary figure 6b

plotSkullQuiver_FS

%% supplementary figure 6c

plotLocomotionTriggeredAvgSkull_FS

%% supplementary figure 6d

plotEMGTriggeredAvgSkull_FS

%% supplementary figure 8a-c

plotMouseStats_FS

%% supplementary figure 9b

plotMotionTrackingBrainOnly_FS(combinedMovementData_210218_006,'Supplementary Figure 9b')

%% supplementary figure 10a-d

plotLocoEMGHist_FS(LTATimeData,LTADispData,ETATimeData,ETADispData)

%% supplementary figure 11a-d

plotSqueezeRespHist_FS(STATimeData,STADispData,RTATimeData,RTADispData)

%% supplementary figure 14a

plotMotionTrackingBrainAndSkullOlf_FS(combinedMovementDataBrain_220816_004,combinedMovementDataSkull_220816_004)

%% supplementary figure 14b

plotLocomotionTriggeredAvgOlf_FS

%% supplementary figure 15b

plotAbdominalCompressionPressure_FS

%% supplementary figure 16a-d

plotXCorrData_FS

%% animation movie

plotMovementBrainInSkullMovie_FS(combinedMovementDataBrain_211216_001,combinedMovementDataSkull_211216_001)

%% olfactory movie

plotMovementBrainInSkullOlfMovie_FS(combinedMovementDataBrain_220813_003,combinedMovementDataSkull_220813_003)

%% 8 locations movie

plotMotionTrackingBrainOnly_FS(combinedMovementData_210218_006,'8 Locations Movie')

%% 3 layer movie

plotMovementBrainInSkull3LayerMovie_FS(combinedMovementDataBrain_221205_004,combinedMovementDataDura_221205_004,combinedMovementDataSkull_221205_004)

%% motion no locomotion movie

plotMovementBrainInSkullNoLocoMovie_FS(combinedMovementDataBrain_221207_008,combinedMovementDataSkull_221207_008)
plotMotionTracking2P2LPCANoLoco_FS(combinedMovementDataBrain_221207_008,combinedMovementDataSkull_221207_008)

%% squeeze movie

plotMotionTrackingBrainAndSkullSqueeze_FS(combinedMovementDataBrain_240612_020,combinedMovementDataSkull_240612_020,'Squeeze Movie 1')
plotMotionTracking2P2LPCASqueeze_FS(combinedMovementDataBrain_240612_020,combinedMovementDataSkull_240612_020,'Squeeze Movie 2')

%% respiration movie

plotMotionTrackingBrainAndSkullResp_FS(combinedMovementDataBrain_240611_001,combinedMovementDataSkull_240611_001,'Respiration Movie')

%% hemisphere comparison movie

plotMovementBrainInSkullHemisphereMovie_FS(combinedMovementDataBrain_221207_010,combinedMovementDataSkull_221207_010,combinedMovementDataBrain_221208_006,combinedMovementDataSkull_221208_006)
plotMotionTracking2P2LPCAHemisphere_FS(combinedMovementDataBrain_221207_010,combinedMovementDataSkull_221207_010,combinedMovementDataBrain_221208_006,combinedMovementDataSkull_221208_006)

%% vectorize figures for publication

% convert all figures to render for vectorization
% h =  findobj('type','figure');
% i = length(h);
% for n = 1:i
%     figure(n)
%     set(gcf,'renderer','Painters')
% end

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
h(6) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
subplot(3,1,1)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
hold off
% title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
xlabel('Time (s)')
ylabel('X Position (\mum)')
grid on
axis([148 398 -3 6])
text(148,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(148,-3,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(3,1,2)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
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
ylabel('Treadmill Velocity (m/s)')
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

function plotMotionTrackingBrainOnly_FS(movementData,figureTitle)

subtitle = 'Figure 1e';
yLimitMotion = [-2 15];
yLimitLocomotion = [0 0.3];

h(1) = figure('Color','White','Name',figureTitle,'NumberTitle','off');

x1 = subplot(3,1,1);
plot((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame,movementData.targetPosition(:,1),'r')
hold on;
f = fill([(1:size(movementData.targetPosition,1))*movementData.secondsPerFrame flip((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame)],movementData.cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
hold off
% title(['\fontsize{20pt}\bf{Figure 1e}'])
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
f = fill([(1:size(movementData.targetPosition,1))*movementData.secondsPerFrame flip((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame)],-1*movementData.cIntFillPtsY,[.5 .5 .5],'Linestyle','none');
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
ylabel('Treadmill Velocity m/s')
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
h(2) = figure('Color','White','Name','Figure 2a','NumberTitle','off');
s = scatter(targetPositionInSkull(:,1),targetPositionInSkull(:,2),10,'g','filled');
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
% title(['\fontsize{20pt}\bf{Figure 2a}'])
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

function plotMotionTracking2P2LPCANoLoco_FS(movementData,stationaryData)
targetPositionInSkull = combineBrainSkullMovement_FS(movementData,stationaryData);
motionVec = pcaMotionAnalysis_FS(targetPositionInSkull);
h(2) = figure('Color','White','Name','Motion No Locomotion Movie 2','NumberTitle','off');
s = scatter(targetPositionInSkull(:,1),targetPositionInSkull(:,2),10,'g','filled');
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
% title(['\fontsize{20pt}\bf{Figure 2a}'])
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

function plotMotionTracking2P2LPCAHemisphere_FS(movementData_1,stationaryData_1,movementData_2,stationaryData_2)
targetPositionInSkull = combineBrainSkullMovement_FS(movementData_1,stationaryData_1);
motionVec = pcaMotionAnalysis_FS(targetPositionInSkull);
h(2) = figure('Color','White','Name','Hemisphere Movie 3','NumberTitle','off');
s = scatter(targetPositionInSkull(:,1),targetPositionInSkull(:,2),10,'g','filled');
s.MarkerFaceAlpha = .1;
hold on
drawArrow_FS([0;0],[motionVec(1);motionVec(2)]);
hold off
axis equal square
axis([-8 8 -8 8])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
% title(['\fontsize{20pt}\bf{Position of Brain in Skull with PCA Vector}'])
% title(['\fontsize{20pt}\bf{Figure 2a}'])
xlabel('\mum')
ylabel('\mum')
if movementData_1.hemisphere == 1
    text(8,0,'Lateral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(-8,0,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    text(0,8,'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(0,-8,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
else
    text(8,0,'Medial','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(-8,0,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    text(0,8,'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(0,-8,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
end

targetPositionInSkull = combineBrainSkullMovement_FS(movementData_2,stationaryData_2);
motionVec = pcaMotionAnalysis_FS(targetPositionInSkull);
h(2) = figure('Color','White','Name','Hemisphere Movie 4','NumberTitle','off');
s = scatter(targetPositionInSkull(:,1),targetPositionInSkull(:,2),10,'g','filled');
s.MarkerFaceAlpha = .1;
hold on
drawArrow_FS([0;0],[-1*motionVec(1);-1*motionVec(2)]);
hold off
axis equal square
axis([-5 5 -5 5])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
% title(['\fontsize{20pt}\bf{Position of Brain in Skull with PCA Vector}'])
% title(['\fontsize{20pt}\bf{Figure 2a}'])
xlabel('\mum')
ylabel('\mum')
if movementData_2.hemisphere == 1
    text(5,0,'Lateral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(-5,0,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    text(0,5,'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(0,-5,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
else
    text(5,0,'Medial','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(-5,0,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    text(0,5,'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(0,-5,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
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

% if ~exist('swapTF','var')
%     swapTF = false;
% end
% if ~exist('vecReverseTF','var')
%     vecReverseTF = false;
% end
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
for n = 1:size(uniqueLocs,1)
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
h(3) = figure('Color','White','Name','Figure 2b','NumberTitle','off');
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
%     plot(moveDataMat{n,4},-moveDataMat{n,5},[arrowColor 'o'],'MarkerSize',5)
    drawArrowQuiver_FS([moveDataMat{n,4},-moveDataMat{n,5}],[moveDataMat{n,4}+(moveDataMat{n,7}(1)*100),-moveDataMat{n,5}+(moveDataMat{n,7}(2)*100)]);
%     quiver(moveDataMat{n,4},-moveDataMat{n,5},moveDataMat{n,7}(1),moveDataMat{n,7}(2),150,'Color',arrowColor,'LineWidth',1.5,'MaxHeadSize',20)
    skullMeanX{catIdx+1}(end+1) = moveDataMat{n,8};
    skullMeanY{catIdx+1}(end+1) = moveDataMat{n,9};
    skullStdX{catIdx+1}(end+1) = moveDataMat{n,10};
    skullStdY{catIdx+1}(end+1) = moveDataMat{n,11};
    brainMeanX{catIdx+1}(end+1) = moveDataMat{n,12};
    brainMeanY{catIdx+1}(end+1) = moveDataMat{n,13};
    brainStdX{catIdx+1}(end+1) = moveDataMat{n,14};
    brainStdY{catIdx+1}(end+1) = moveDataMat{n,15};
end
xlim([-5000 5000])
ylim([-5000 5000])
axis square
plot(0,0,'ko','MarkerSize',5)
plot(0,-2600,'ko','MarkerSize',5)
% quiver(-2500,2500,3,0,150,'LineWidth',1.5,'MaxHeadSize',20)
drawArrowQuiver_FS([-2500,2500],[-2500+(3*100),2500]);
% title(['Figure 2b' 10 'Rostral, \mum'])
ylabel('\mum')
xlabel('\mum')
text(55,55,'Bregma')
text(55,-2500,'Lambda')
text(-2500,2300,'3 \mum')
text(0,5000,'Rostral')
text(0,-4500,'Caudal')
text(5000,0,'Right')
text(-5000,0,'Left')
rectangle('Position',[-2650 2000 800 800])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotSkullQuiver_FS
load('movementDataLogSkull_FS.mat')
moveDataMat = moveDataMatSkull;
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
for n = 1:size(uniqueLocs,1)
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
h(3) = figure('Color','White','Name','Supplementary Figure 6b','NumberTitle','off');
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
%     plot(moveDataMat{n,4},-moveDataMat{n,5},[arrowColor 'o'],'MarkerSize',5)
    drawArrowQuiver_FS([moveDataMat{n,4},-moveDataMat{n,5}],[moveDataMat{n,4}+(moveDataMat{n,7}(1)*100),-moveDataMat{n,5}+(moveDataMat{n,7}(2)*100)]);
%     quiver(moveDataMat{n,4},-moveDataMat{n,5},moveDataMat{n,7}(1),moveDataMat{n,7}(2),150,'Color',arrowColor,'LineWidth',1.5,'MaxHeadSize',20)
    skullMeanX{catIdx+1}(end+1) = moveDataMat{n,8};
    skullMeanY{catIdx+1}(end+1) = moveDataMat{n,9};
    skullStdX{catIdx+1}(end+1) = moveDataMat{n,10};
    skullStdY{catIdx+1}(end+1) = moveDataMat{n,11};
    brainMeanX{catIdx+1}(end+1) = moveDataMat{n,12};
    brainMeanY{catIdx+1}(end+1) = moveDataMat{n,13};
    brainStdX{catIdx+1}(end+1) = moveDataMat{n,14};
    brainStdY{catIdx+1}(end+1) = moveDataMat{n,15};
end
xlim([-5000 5000])
ylim([-5000 5000])
axis square
plot(0,0,'ko','MarkerSize',5)
plot(0,-2600,'ko','MarkerSize',5)
drawArrowQuiver_FS([-2500,2500],[-2500+(3*100),2500]);
% quiver(-2500,2500,3,0,150,'LineWidth',1.5,'MaxHeadSize',20)
% title(['Figure 2b' 10 'Rostral, \mum'])
ylabel('\mum')
xlabel('\mum')
text(55,55,'Bregma')
text(55,-2500,'Lambda')
text(-2500,2300,'3 \mum')
text(0,5000,'Rostral')
text(0,-4500,'Caudal')
text(5000,0,'Right')
text(-5000,0,'Left')
rectangle('Position',[-2650 2000 800 800])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotSqueezeQuiver_FS
load('movementDataLogSqueeze_FS.mat')
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
for n = 1:size(uniqueLocs,1)
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
h(3) = figure('Color','White','Name','Figure 5d','NumberTitle','off');
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
%     plot(moveDataMat{n,4},-moveDataMat{n,5},[arrowColor 'o'],'MarkerSize',5)
    drawArrowQuiver_FS([moveDataMat{n,4},-moveDataMat{n,5}],[moveDataMat{n,4}+(moveDataMat{n,7}(1)*100),-moveDataMat{n,5}+(moveDataMat{n,7}(2)*100)]);
%     quiver(moveDataMat{n,4},-moveDataMat{n,5},moveDataMat{n,7}(1),moveDataMat{n,7}(2),150,'Color',arrowColor,'LineWidth',1.5,'MaxHeadSize',20)
    skullMeanX{catIdx+1}(end+1) = moveDataMat{n,8};
    skullMeanY{catIdx+1}(end+1) = moveDataMat{n,9};
    skullStdX{catIdx+1}(end+1) = moveDataMat{n,10};
    skullStdY{catIdx+1}(end+1) = moveDataMat{n,11};
    brainMeanX{catIdx+1}(end+1) = moveDataMat{n,12};
    brainMeanY{catIdx+1}(end+1) = moveDataMat{n,13};
    brainStdX{catIdx+1}(end+1) = moveDataMat{n,14};
    brainStdY{catIdx+1}(end+1) = moveDataMat{n,15};
end
xlim([-5000 5000])
ylim([-5000 5000])
axis square
plot(0,0,'ko','MarkerSize',5)
plot(0,-2600,'ko','MarkerSize',5)
% quiver(-2500,2500,3,0,150,'LineWidth',1.5,'MaxHeadSize',20)
drawArrowQuiver_FS([-2500,2500],[-2500+(3*100),2500]);
% title(['Figure 2b' 10 'Rostral, \mum'])
ylabel('\mum')
xlabel('\mum')
text(55,55,'Bregma')
text(55,-2500,'Lambda')
text(-2500,2300,'3 \mum')
text(0,5000,'Rostral')
text(0,-4500,'Caudal')
text(5000,0,'Right')
text(-5000,0,'Left')
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
% % f = fill([HzX flip(HzX)],[ErrorX(1,:) flip(ErrorX(2,:))],[.5 .5 .5],'Linestyle','none');
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
% % f = fill([HzT flip(HzT)],[ErrorT(1,:) flip(ErrorT(2,:))],[.5 .5 .5],'Linestyle','none');
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
% % f = fill([fx flip(fx)],[cErrx(1,:) flip(cErrx(2,:))],[.5 .5 .5],'Linestyle','none');
% % set(f,'facea',[.2]);
% yline(confCx);
% hold off
% text(.1,.8,['Conf = ' num2str(confCx)])
% title('Coherence')

h(4) = figure('Color','White','Name','Figure 2c','NumberTitle','off');
subplot(3,1,1)
PowerY = downsample(PowerY,3);
HzY = downsample(HzY,3);
semilogx(HzY,PowerY,'b')
xlabel('Frequency (Hz)')
ylabel('A.U.')
title(['Brain Motion (R/C)'])
% hold on
% f = fill([HzY flip(HzY)],[ErrorY(1,:) flip(ErrorY(2,:))],[.5 .5 .5],'Linestyle','none');
% set(f,'facea',[.2]);
% hold off
subplot(3,1,2)
semilogx(HzT,PowerT*100,'Color',[0.8549 0.1098 0.3607])
xlabel('Frequency (Hz)')
ylabel('A.U.')
title('Thermocouple (Respiration)')
% hold on
% f = fill([HzT flip(HzT)],[ErrorT(1,:) flip(ErrorT(2,:))],[.5 .5 .5],'Linestyle','none');
% set(f,'facea',[.2]);
% hold off
subplot(3,1,3)
Cy = downsample(Cy,3);
fy = downsample(fy,3);
semilogx(fy,Cy,'k');
ylim([0 1])
xlabel('Frequency (Hz)')
ylabel('Coherence')
title('Thermocouple Motion Coherence')
hold on
% f = fill([fy flip(fy)],[cErry(1,:) flip(cErry(2,:))],[.5 .5 .5],'Linestyle','none');
% set(f,'facea',[.2]);
yline(confCx);
hold off
% text(.1,.8,['Conf = ' num2str(confCy)])
title('Thermocouple Motion Coherence')
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
    h(1) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
    semilogy(HzX,PowerX,'k')
    hold on
    f = fill([HzX flip(HzX)],[ErrorX(1,:) flip(ErrorX(2,:))],[.5 .5 .5],'Linestyle','none');
    set(f,'facea',[.2]);
    hold off
    title(['\fontsize{20pt}\bf{X Position Frequency Domain}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Frequency (Hz)')
    ylabel('Power')
    
    h(2) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
    semilogy(HzY,PowerY,'k')
    hold on
    f = fill([HzY flip(HzY)],[ErrorY(1,:) flip(ErrorY(2,:))],[.5 .5 .5],'Linestyle','none');
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
movementData.locDataInterp=zeros(size(movementData.targetPosition,1));
movementData.locDataInterp(:,1)=movement_time;
movementData.locDataInterp(:,2)=interp1(movementData.ballData(:,1),abs(movementData.ballData(:,2)),movement_time,'linear');

h(5) = figure('Color','White','Name','Figure 2d','NumberTitle','off');
maxlag=500;
xc_1=xcorr(detrend(movementData.locDataInterp(1:(end-100),2))', detrend(movementData.targetPosition(1:(end-100),1))',maxlag,'coeff');
xc_2=xcorr(detrend(movementData.locDataInterp(1:(end-100),2))', detrend(movementData.targetPosition(1:(end-100),2))',maxlag,'coeff');
plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_1,'r')
hold on
plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_2,'b')
hold off
title('Brain Motion and Locomotion Cross-Correlation')
ylabel('Noramlized Cross-Correlation')
xlabel('Lags (s)')
legend({'Lateral-Medial','Rostral-Caudal'})
ylim([0 1])
axes('Position',[.2 .7 .2 .2])
box on
plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_1,'r')
hold on
plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_2,'b')
hold off
axis([-10 10 .4 .7])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMovementBrainInSkull_FS(movementData,stationaryData,rawData)
targetPositionInSkull = combineBrainSkullMovement_FS(movementData,stationaryData);
movementData.secondsPerFrame = movementData.secondsPerFrame/2;
h(6) = figure('Color','White','Name','Figure 3d','NumberTitle','off');
subplot(4,1,1)
plot([1:size(targetPositionInSkull,1)]*movementData.secondsPerFrame,targetPositionInSkull(:,1),'r')
% title(['Figure 3c' 10 '\fontsize{20pt}\bf{Position of Brain in Skull}'])
xlabel('Time (s)')
ylabel('Brain Shift (\mum)')
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
plot([1:size(targetPositionInSkull,1)]*movementData.secondsPerFrame,targetPositionInSkull(:,2),'b')
xlabel('Time (s)')
ylabel('Brain Shift (\mum)')
grid on
axis([0 600 -2 6])
text(0,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(0,-2,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(4,1,3)
if all(movementData.emgData(:,2) == 0)
    title('\fontsize{20pt}\bf{No EMG Data}')
else
    semilogy(movementData.emgData(:,1),movementData.emgData(:,2),'Color',[1 0.5 0])
    xlabel('Time (s)')
    ylabel('EMG Power')
    grid on
    xlim([0 600])
end
subplot(4,1,4)
plot(rawData(:,1),rawData(:,3),'Color',[1 0.5 0])
axis([0 600 -1 1])
xlabel('Time (s)')
ylabel('EMG (mV)')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMovementBrainInSkullMovie_FS(movementData,stationaryData)
targetPositionInSkull = combineBrainSkullMovement_FS(movementData,stationaryData);
movementData.secondsPerFrame = movementData.secondsPerFrame/2;
h(6) = figure('Color','White','Name','Animation Movie','NumberTitle','off');
subplot(5,1,1)
plot([1:size(targetPositionInSkull,1)]*movementData.secondsPerFrame,targetPositionInSkull(:,1),'r')
% title(['Figure 3c' 10 '\fontsize{20pt}\bf{Position of Brain in Skull}'])
xlabel('Time (s)')
ylabel('Brain Shift (\mum)')
grid on
axis([200 220 -4 6])
if movementData.hemisphere == 1
    text(200,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(200,-4,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
else
    text(200,6,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(200,-4,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
end
subplot(5,1,2)
plot([1:size(targetPositionInSkull,1)]*movementData.secondsPerFrame,targetPositionInSkull(:,2),'b')
xlabel('Time (s)')
ylabel('Brain Shift (\mum)')
grid on
axis([200 220 -4 6])
text(200,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(200,-4,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(5,1,3)
if all(movementData.emgData(:,2) == 0)
    title('\fontsize{20pt}\bf{No EMG Data}')
else
    semilogy(movementData.emgData(:,1),movementData.emgData(:,2),'Color',[1 0.5 0])
    xlabel('Time (s)')
    ylabel('EMG Power')
    grid on
    xlim([200 220])
end
subplot(5,1,4)
plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)*2*pi*.06),'k')
axis([200 220 0 .25])
xlabel('Time (s)')
ylabel('Treadmill Velocity (m/s)')
subplot(5,1,5)
plot(movementData.thermoData(:,1),movementData.thermoData(:,2),'Color',[0.8549 0.1098 0.3607])
axis([200 220 -.5 .5])
xlabel('Time (s)')
ylabel('Respiration Thermocouple (au)')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMovementBrainInSkullOlfMovie_FS(movementData,stationaryData)
[movementData.targetPosition,stationaryData.targetPosition] = interpBrainSkullMovement_FS(movementData,stationaryData);
if movementData.hemisphere == 2
    movementData.targetPosition(:,1) = movementData.targetPosition(:,1)*-1;
    stationaryData.targetPosition(:,1) = stationaryData.targetPosition(:,1)*-1;
end
movementData.targetPosition(:,2) = movementData.targetPosition(:,2)*-1;
stationaryData.targetPosition(:,2) = stationaryData.targetPosition(:,2)*-1;
% movementData.targetPosition(:,1) = movementData.targetPosition(:,1)-mean(movementData.targetPosition(2480:2530,1));
% movementData.targetPosition(:,2) = movementData.targetPosition(:,2)-mean(movementData.targetPosition(2480:2530,2));
% stationaryData.targetPosition(:,1) = stationaryData.targetPosition(:,1)-mean(stationaryData.targetPosition(2480:2530,1));
% stationaryData.targetPosition(:,2) = stationaryData.targetPosition(:,2)-mean(stationaryData.targetPosition(2480:2530,2));
movementData.secondsPerFrame = movementData.secondsPerFrame/2;
h(6) = figure('Color','White','Name','Olfactory Movie','NumberTitle','off');
subplot(3,1,1)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
hold off
% title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
xlabel('Time (s)')
ylabel('Medial-Lateral Shift (\mum)')
grid on
axis([100 150 -4 4])
text(100,4,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(100,-4,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(3,1,2)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
hold off
xlabel('Time (s)')
ylabel('Rostral-Caudal Shift (\mum)')
grid on
axis([100 150 -4 4])
text(100,4,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(100,-4,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(4,1,3)
% plot(movementData.emgData(:,1),movementData.emgData(:,2),'k')
% xlabel('Time (s)')
% ylabel('Abdominal EMG (au)')
% grid on
% axis([62 272 0.5 3])
subplot(3,1,3)
plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)*2*pi*.06),'k')
xlabel('Time (s)')
ylabel('Treadmill Velocity (m/s)')
grid on
axis([100 150 0 0.2])

% h(6) = figure('Color','White','Name','Animation Movie','NumberTitle','off');
% subplot(3,1,1)
% plot([1:size(targetPositionInSkull,1)]*movementData.secondsPerFrame,targetPositionInSkull(:,1),'r')
% % title(['Figure 3c' 10 '\fontsize{20pt}\bf{Position of Brain in Skull}'])
% xlabel('Time (s)')
% ylabel('Brain Shift (\mum)')
% grid on
% axis([100 150 -4 4])
% if movementData.hemisphere == 1
%     text(200,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
%     text(200,-4,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% else
%     text(200,6,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
%     text(200,-4,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% end
% subplot(3,1,2)
% plot([1:size(targetPositionInSkull,1)]*movementData.secondsPerFrame,targetPositionInSkull(:,2),'b')
% xlabel('Time (s)')
% ylabel('Brain Shift (\mum)')
% grid on
% axis([100 150 -4 4])
% text(200,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(200,-4,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(3,1,3)
% plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)*2*pi*.06),'k')
% axis([100 150 0 .2])
% xlabel('Time (s)')
% ylabel('Treadmill Velocity (m/s)')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMovementBrainInSkull3LayerMovie_FS(movementData,duraData,stationaryData)
% [movementData.targetPosition,stationaryData.targetPosition] = interpBrainSkullMovement_FS(movementData,stationaryData);
    movementData.targetPosition(:,1) = movementData.targetPosition(:,1)*-1;
    stationaryData.targetPosition(:,1) = stationaryData.targetPosition(:,1)*-1;
    duraData.targetPosition(:,1) = duraData.targetPosition(:,1)*-1;
movementData.targetPosition(:,2) = movementData.targetPosition(:,2)*-1;
stationaryData.targetPosition(:,2) = stationaryData.targetPosition(:,2)*-1;
duraData.targetPosition(:,2) = duraData.targetPosition(:,2)*-1;
% movementData.targetPosition(:,1) = movementData.targetPosition(:,1)-mean(movementData.targetPosition(2480:2530,1));
% movementData.targetPosition(:,2) = movementData.targetPosition(:,2)-mean(movementData.targetPosition(2480:2530,2));
% stationaryData.targetPosition(:,1) = stationaryData.targetPosition(:,1)-mean(stationaryData.targetPosition(2480:2530,1));
% stationaryData.targetPosition(:,2) = stationaryData.targetPosition(:,2)-mean(stationaryData.targetPosition(2480:2530,2));
% movementData.secondsPerFrame = movementData.secondsPerFrame/2;
h(6) = figure('Color','White','Name','3 Layer Movie','NumberTitle','off');
subplot(3,1,1)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
plot([1:size(duraData.targetPosition,1)]*movementData.secondsPerFrame,duraData.targetPosition(:,1),'Color',[.9 .9 .9])
hold off
% title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
xlabel('Time (s)')
ylabel('Medial-Lateral Shift (\mum)')
grid on
axis([75 165 -1 5])
text(75,5,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(75,-1,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(3,1,2)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
plot([1:size(duraData.targetPosition,1)]*movementData.secondsPerFrame,duraData.targetPosition(:,2),'Color',[.9 .9 .9])
hold off
xlabel('Time (s)')
ylabel('Rostral-Caudal Shift (\mum)')
grid on
axis([75 165 -1 5])
text(75,5,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(75,-1,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(4,1,3)
% plot(movementData.emgData(:,1),movementData.emgData(:,2),'k')
% xlabel('Time (s)')
% ylabel('Abdominal EMG (au)')
% grid on
% axis([62 272 0.5 3])
subplot(3,1,3)
plot(movementData.emgData(:,1),movementData.emgData(:,2),'Color',[1 0.5 0])
xlabel('Time (s)')
ylabel('EMG Power')
grid on
axis([75 165 0 4])

% h(6) = figure('Color','White','Name','Animation Movie','NumberTitle','off');
% subplot(3,1,1)
% plot([1:size(targetPositionInSkull,1)]*movementData.secondsPerFrame,targetPositionInSkull(:,1),'r')
% % title(['Figure 3c' 10 '\fontsize{20pt}\bf{Position of Brain in Skull}'])
% xlabel('Time (s)')
% ylabel('Brain Shift (\mum)')
% grid on
% axis([100 150 -4 4])
% if movementData.hemisphere == 1
%     text(200,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
%     text(200,-4,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% else
%     text(200,6,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
%     text(200,-4,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% end
% subplot(3,1,2)
% plot([1:size(targetPositionInSkull,1)]*movementData.secondsPerFrame,targetPositionInSkull(:,2),'b')
% xlabel('Time (s)')
% ylabel('Brain Shift (\mum)')
% grid on
% axis([100 150 -4 4])
% text(200,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(200,-4,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(3,1,3)
% plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)*2*pi*.06),'k')
% axis([100 150 0 .2])
% xlabel('Time (s)')
% ylabel('Treadmill Velocity (m/s)')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMovementBrainInSkullHemisphereMovie_FS(movementData_1,stationaryData_1,movementData_2,stationaryData_2)
% [movementData.targetPosition,stationaryData.targetPosition] = interpBrainSkullMovement_FS(movementData,stationaryData);
%     movementData_1.targetPosition(:,1) = movementData_1.targetPosition(:,1)*-1;
%     stationaryData_1.targetPosition(:,1) = stationaryData_1.targetPosition(:,1)*-1;
movementData_1.targetPosition(:,2) = movementData_1.targetPosition(:,2)*-1;
stationaryData_1.targetPosition(:,2) = stationaryData_1.targetPosition(:,2)*-1;
% movementData_1.targetPosition(:,1) = movementData_1.targetPosition(:,1)-mean(movementData_1.targetPosition(2480:2530,1));
% movementData_1.targetPosition(:,2) = movementData_1.targetPosition(:,2)-mean(movementData_1.targetPosition(2480:2530,2));
% stationaryData_1.targetPosition(:,1) = stationaryData_1.targetPosition(:,1)-mean(stationaryData_1.targetPosition(2480:2530,1));
% stationaryData_1.targetPosition(:,2) = stationaryData_1.targetPosition(:,2)-mean(stationaryData_1.targetPosition(2480:2530,2));
% movementData_1.secondsPerFrame = movementData_1.secondsPerFrame/2;
h(6) = figure('Color','White','Name','Hemisphere Movie 1','NumberTitle','off');
subplot(3,1,1)
plot([1:size(movementData_1.targetPosition,1)]*movementData_1.secondsPerFrame,movementData_1.targetPosition(:,1),'g')
hold on
plot([1:size(stationaryData_1.targetPosition,1)]*movementData_1.secondsPerFrame,stationaryData_1.targetPosition(:,1),'m')
hold off
% title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
xlabel('Time (s)')
ylabel('Medial-Lateral Shift (\mum)')
grid on
axis([150 400 -5 5])
text(150,5,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(150,-5,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(3,1,2)
plot([1:size(movementData_1.targetPosition,1)]*movementData_1.secondsPerFrame,movementData_1.targetPosition(:,2),'g')
hold on
plot([1:size(stationaryData_1.targetPosition,1)]*movementData_1.secondsPerFrame,stationaryData_1.targetPosition(:,2),'m')
hold off
xlabel('Time (s)')
ylabel('Rostral-Caudal Shift (\mum)')
grid on
axis([150 400 -5 5])
text(150,5,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(150,-5,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(4,1,3)
% plot(movementData_1.emgData(:,1),movementData_1.emgData(:,2),'k')
% xlabel('Time (s)')
% ylabel('Abdominal EMG (au)')
% grid on
% axis([62 272 0.5 3])
subplot(3,1,3)
plot(movementData_1.emgData(:,1),movementData_1.emgData(:,2),'Color',[1 0.5 0])
xlabel('Time (s)')
ylabel('EMG Power')
grid on
axis([150 400 0 4])

% [movementData.targetPosition,stationaryData.targetPosition] = interpBrainSkullMovement_FS(movementData,stationaryData);
% movementData_2.targetPosition(:,1) = movementData_2.targetPosition(:,1)*-1;
% stationaryData_2.targetPosition(:,1) = stationaryData_2.targetPosition(:,1)*-1;
movementData_2.targetPosition(:,2) = movementData_2.targetPosition(:,2)*-1;
stationaryData_2.targetPosition(:,2) = stationaryData_2.targetPosition(:,2)*-1;
% movementData_2.targetPosition(:,1) = movementData_2.targetPosition(:,1)-mean(movementData_2.targetPosition(2480:2530,1));
% movementData_2.targetPosition(:,2) = movementData_2.targetPosition(:,2)-mean(movementData_2.targetPosition(2480:2530,2));
% stationaryData_2.targetPosition(:,1) = stationaryData_2.targetPosition(:,1)-mean(stationaryData_2.targetPosition(2480:2530,1));
% stationaryData_2.targetPosition(:,2) = stationaryData_2.targetPosition(:,2)-mean(stationaryData_2.targetPosition(2480:2530,2));
% movementData_2.secondsPerFrame = movementData_2.secondsPerFrame/2;
h(6) = figure('Color','White','Name','Hemisphere Movie 2','NumberTitle','off');
subplot(3,1,1)
plot([1:size(movementData_2.targetPosition,1)]*movementData_2.secondsPerFrame,movementData_2.targetPosition(:,1),'g')
hold on
plot([1:size(stationaryData_2.targetPosition,1)]*movementData_2.secondsPerFrame,stationaryData_2.targetPosition(:,1),'m')
hold off
% title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
xlabel('Time (s)')
ylabel('Medial-Lateral Shift (\mum)')
grid on
axis([90 340 -5 5])
text(90,5,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(90,-5,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(3,1,2)
plot([1:size(movementData_2.targetPosition,1)]*movementData_2.secondsPerFrame,movementData_2.targetPosition(:,2),'g')
hold on
plot([1:size(stationaryData_2.targetPosition,1)]*movementData_2.secondsPerFrame,stationaryData_2.targetPosition(:,2),'m')
hold off
xlabel('Time (s)')
ylabel('Rostral-Caudal Shift (\mum)')
grid on
axis([90 340 -5 5])
text(90,5,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(90,-5,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(4,1,3)
% plot(movementData_2.emgData(:,1),movementData_2.emgData(:,2),'k')
% xlabel('Time (s)')
% ylabel('Abdominal EMG (au)')
% grid on
% axis([62 272 0.5 3])
subplot(3,1,3)
plot(movementData_2.emgData(:,1),movementData_2.emgData(:,2),'Color',[1 0.5 0])
xlabel('Time (s)')
ylabel('EMG Power')
grid on
axis([90 340 0 4])

% h(6) = figure('Color','White','Name','Animation Movie','NumberTitle','off');
% subplot(3,1,1)
% plot([1:size(targetPositionInSkull,1)]*movementData.secondsPerFrame,targetPositionInSkull(:,1),'r')
% % title(['Figure 3c' 10 '\fontsize{20pt}\bf{Position of Brain in Skull}'])
% xlabel('Time (s)')
% ylabel('Brain Shift (\mum)')
% grid on
% axis([100 150 -4 4])
% if movementData.hemisphere == 1
%     text(200,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
%     text(200,-4,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% else
%     text(200,6,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
%     text(200,-4,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% end
% subplot(3,1,2)
% plot([1:size(targetPositionInSkull,1)]*movementData.secondsPerFrame,targetPositionInSkull(:,2),'b')
% xlabel('Time (s)')
% ylabel('Brain Shift (\mum)')
% grid on
% axis([100 150 -4 4])
% text(200,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(200,-4,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(3,1,3)
% plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)*2*pi*.06),'k')
% axis([100 150 0 .2])
% xlabel('Time (s)')
% ylabel('Treadmill Velocity (m/s)')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMovementBrainInSkullNoLocoMovie_FS(movementData,stationaryData)
[movementData.targetPosition,stationaryData.targetPosition] = interpBrainSkullMovement_FS(movementData,stationaryData);
if movementData.hemisphere == 2
    movementData.targetPosition(:,1) = movementData.targetPosition(:,1)*-1;
    stationaryData.targetPosition(:,1) = stationaryData.targetPosition(:,1)*-1;
end
movementData.targetPosition(:,2) = movementData.targetPosition(:,2)*-1;
stationaryData.targetPosition(:,2) = stationaryData.targetPosition(:,2)*-1;
% movementData.targetPosition(:,1) = movementData.targetPosition(:,1)-mean(movementData.targetPosition(2480:2530,1));
% movementData.targetPosition(:,2) = movementData.targetPosition(:,2)-mean(movementData.targetPosition(2480:2530,2));
% stationaryData.targetPosition(:,1) = stationaryData.targetPosition(:,1)-mean(stationaryData.targetPosition(2480:2530,1));
% stationaryData.targetPosition(:,2) = stationaryData.targetPosition(:,2)-mean(stationaryData.targetPosition(2480:2530,2));
movementData.secondsPerFrame = movementData.secondsPerFrame/2;
h(6) = figure('Color','White','Name','Motion No Locomotion Movie 1','NumberTitle','off');
subplot(4,1,1)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
hold off
% title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
xlabel('Time (s)')
ylabel('Medial-Lateral Shift (\mum)')
grid on
axis([280 350 -1 6])
text(100,4,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(100,-4,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(4,1,2)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
hold off
xlabel('Time (s)')
ylabel('Rostral-Caudal Shift (\mum)')
grid on
axis([280 350 -1 6])
text(100,4,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(100,-4,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(4,1,3)
plot(movementData.emgData(:,1),movementData.emgData(:,2),'Color',[1 0.5 0])
xlabel('Time (s)')
ylabel('Abdominal EMG (au)')
grid on
axis([280 350 1 3])
subplot(4,1,4)
plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)*2*pi*.06),'k')
xlabel('Time (s)')
ylabel('Treadmill Velocity (m/s)')
grid on
axis([280 350 0 .4])

% h(6) = figure('Color','White','Name','Animation Movie','NumberTitle','off');
% subplot(3,1,1)
% plot([1:size(targetPositionInSkull,1)]*movementData.secondsPerFrame,targetPositionInSkull(:,1),'r')
% % title(['Figure 3c' 10 '\fontsize{20pt}\bf{Position of Brain in Skull}'])
% xlabel('Time (s)')
% ylabel('Brain Shift (\mum)')
% grid on
% axis([100 150 -4 4])
% if movementData.hemisphere == 1
%     text(200,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
%     text(200,-4,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% else
%     text(200,6,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
%     text(200,-4,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% end
% subplot(3,1,2)
% plot([1:size(targetPositionInSkull,1)]*movementData.secondsPerFrame,targetPositionInSkull(:,2),'b')
% xlabel('Time (s)')
% ylabel('Brain Shift (\mum)')
% grid on
% axis([100 150 -4 4])
% text(200,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(200,-4,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(3,1,3)
% plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)*2*pi*.06),'k')
% axis([100 150 0 .2])
% xlabel('Time (s)')
% ylabel('Treadmill Velocity (m/s)')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotEMGMovement_FS(movementData,stationaryData)
movementData.targetPosition = combineBrainSkullMovement_FS(movementData,stationaryData);
movementData.secondsPerFrame = movementData.secondsPerFrame/2;
best_lag=round(1/movementData.secondsPerFrame);
movement_time=movementData.secondsPerFrame*(1:length(movementData.targetPosition));
movementData.emgDataInterp=zeros(size(movementData.targetPosition,1));
movementData.emgDataInterp(:,1)=movement_time;
movementData.emgDataInterp(:,2)=interp1(movementData.emgData(:,1),movementData.emgData(:,2),movement_time,'linear');
% emg_bins=.5:.1:3.5;
% pixel_bins=-3:.5:7;
emg_bins=.5:.025:3.5;
pixel_bins=-3:.1:7;
h(7) = figure('Color','White','Name','Figure 3e','NumberTitle','off');
% % % subplot(2,2,2)
% % % scatter(movementData.emgDataInterp(1:(end-best_lag),2),movementData.targetPosition(best_lag:end-1,2))
% % % title('Figure 3d')
% % % xlabel('Log EMG Power')
% % % ylabel('Brain Shift, \mum')
% % % text(3.5,6,'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
% % % text(3.5,-2,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
% % % ylim([-2 6])
% % % xlim([0.5 3.5])

% % % subplot(2,2,1)
% % % scatter(movementData.emgDataInterp(1:(end-best_lag),2),movementData.targetPosition(best_lag:end-1,1))
% % % title('Figure 3d')
% % % xlabel('Log EMG Power')
% % % ylabel('Brain Shift, \mum')
% % % text(3.5,6,'Lateral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
% % % text(3.5,-2,'Medial','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
% % % ylim([-2 6])
% % % xlim([0.5 3.5])

subplot(2,1,1)
hh2=histogram2(movementData.emgDataInterp(1:(end-best_lag),2),movementData.targetPosition(best_lag:end-1,1), emg_bins, pixel_bins,...
    'DisplayStyle','tile','ShowEmptyBins','on');
imagesc(emg_bins, pixel_bins,log(hh2.Values'+1))
axis xy
colormap hot
xlabel('EMG Power')
ylabel('Brain Shift, \mum')
axis square
ylim([-2 6])
xlim([0.5 3.5])
text(.5,5.5,'Lateral','Color','w')
text(.5,-1.5,'Medial','Color','w')

subplot(2,1,2)
hh2=histogram2(movementData.emgDataInterp(1:(end-best_lag),2),movementData.targetPosition(best_lag:end-1,2), emg_bins, pixel_bins,...
    'DisplayStyle','tile','ShowEmptyBins','on');
imagesc(emg_bins, pixel_bins,log(hh2.Values'+1))
axis xy
colormap hot
xlabel('EMG Power')
ylabel('Brain Shift, \mum')
axis square
ylim([-2 6])
xlim([0.5 3.5])
text(.5,5.5,'Rostral','Color','w')
text(.5,-1.5,'Caudal','Color','w')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotEMGXCorr_FS(movementData,stationaryData)
targetPositionInSkull = combineBrainSkullMovement_FS(movementData,stationaryData);
movementData.secondsPerFrame = movementData.secondsPerFrame/2;
movementData.targetPosition = targetPositionInSkull;
movement_time=movementData.secondsPerFrame*(1:length(movementData.targetPosition));
movementData.emgDataInterp=zeros(size(movementData.targetPosition,1));
movementData.emgDataInterp(:,1)=movement_time;
movementData.emgDataInterp(:,2)=interp1(movementData.emgData(:,1),abs(movementData.emgData(:,2)),movement_time,'linear');

h(8) = figure('Color','White','Name','Figure 3f','NumberTitle','off');
maxlag=500;
xc_1=xcorr(detrend(movementData.emgDataInterp(1:(end-100),2))', detrend(movementData.targetPosition(1:(end-100),1))',maxlag,'coeff');
xc_2=xcorr(detrend(movementData.emgDataInterp(1:(end-100),2))', detrend(movementData.targetPosition(1:(end-100),2))',maxlag,'coeff');
plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_1,'r')
hold on
plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_2,'b')
hold off
title('Brain Motion and EMG Cross-Correlation')
ylabel('Noramlized Cross-Correlation')
xlabel('Lags (s)')
legend({'Medial-Lateral','Rostral-Caudal'})
ylim([0 1])
axes('Position',[.2 .7 .2 .2])
box on
plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_1,'r')
hold on
plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_2,'b')
hold off
axis([-10 10 .5 .8])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [LTATimeData, LTADispData] = plotLocomotionTriggeredAvg_FS
motionEventsLocationsX = [];
motionEventsLocationsY = [];
timeToThreshX = [];
timeToThreshY = [];
dispTimeThreshX = [];
dispTimeThreshY = [];
moveThresh = .75;
timeThresh = 1.5;
binEdgesTime = -2:.25:3;
binEdgesDisp = -4:.25:4;
load('LTADataCell_FS.mat')
for n = 1:size(locDataCell,1)
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

h(9) = figure('Color','White','Name','Figure 2e','NumberTitle','off');
subplot(2,2,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
f = fill([3 0 0 3],[2.9 2.9 -2.9 -2.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecX,meanX,'k')
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-3 3],'r')
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[1,0,0,0.1])
end
% plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean Motion During Locomotion Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Brian Shift (\mum)')
ylim([-3 3])
xlim([-2 3])
grid on

subplot(2,2,3)
f = fill([3 0 0 3],[2.9 2.9 -2.9 -2.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecY,meanY,'k')
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-3 3],'r')
for n = 1:size(motionEventsLocationsY,1)
    plot(timeVecY,-1*motionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
end
% plot(timeVecY(brainMotionStart),meanY(brainMotionStart),'rx')
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
for n = 1:size(locDataCell,1)
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
f = fill([-2 0 0 -2],[2.9 2.9 -2.9 -2.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecX,meanX,'k')
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-3 3],'r')
for n = 1:size(stopMotionEventsLocationsX,1)
    plot(timeVecX,stopMotionEventsLocationsX(n,:),'Color',[1,0,0,0.1])
end
hold off
text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% title(['\fontsize{20pt}\bf{Mean Motion During Stopping Locomotion Events, n = ' num2str(size(stopMotionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Brian Shift (\mum)')
ylim([-3 3])
xlim([-2 3])
grid on

subplot(2,2,4)
f = fill([-2 0 0 -2],[2.9 2.9 -2.9 -2.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecY,meanY,'k')
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-3 3],'r')
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

% h(10) = figure('Color','White','Name','Figure 2e','NumberTitle','off');
% subplot(2,1,1)
% histogram(timeToThreshX,binEdgesTime);
% hold on
[pdfXVals,pdfYVals] = findKernelPDF(timeToThreshX,binEdgesTime);
% plot(pdfXVals,pdfYVals*20,'r','LineWidth',2)
% title(['Time for brain to displace laterally ' num2str(moveThresh) ' micrometers following locomotion trigger'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 60])
% mu = mean(timeToThreshX);
% sig = std(timeToThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
LTATimeData.timeToThreshX = timeToThreshX;
LTATimeData.binEdgesTime = binEdgesTime;
LTATimeData.X_pdfXVals = pdfXVals;
LTATimeData.X_pdfYVals = pdfYVals*20;


% subplot(2,1,2)
% histogram(timeToThreshY,binEdgesTime);
% hold on
[pdfXVals,pdfYVals] = findKernelPDF(timeToThreshY,binEdgesTime);
% plot(pdfXVals,pdfYVals*40,'r','LineWidth',2)
% title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following locomotion trigger'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 60])
% mu = mean(timeToThreshY);
% sig = std(timeToThreshY);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
LTATimeData.timeToThreshY = timeToThreshY;
LTATimeData.Y_pdfXVals = pdfXVals;
LTATimeData.Y_pdfYVals = pdfYVals*40;

% h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
% histogram(dispTimeThreshX,binEdgesDisp);
% hold on
[pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshX,binEdgesDisp);
% plot(pdfXVals,pdfYVals*40,'r','LineWidth',2)
% title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following locomotion trigger'])
% xlabel('Displacement (\mum)')
% xlim([-4 4])
% ylim([0 60])
% mu = mean(dispTimeThreshX);
% sig = std(dispTimeThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
LTADispData.dispTimeThreshX = dispTimeThreshX;
LTADispData.binEdgesDisp = binEdgesDisp;
LTADispData.X_pdfXVals = pdfXVals;
LTADispData.X_pdfYVals = pdfYVals*40;

% subplot(2,1,2)
% histogram(dispTimeThreshY,binEdgesDisp);
% hold on
[pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshY,binEdgesDisp);
% plot(pdfXVals,pdfYVals*40,'r','LineWidth',2)
% title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following locomotion trigger'])
% xlabel('Displacement (\mum)')
% xlim([-4 4])
% ylim([0 60])
% mu = mean(dispTimeThreshY);
% sig = std(dispTimeThreshY);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
LTADispData.dispTimeThreshY = dispTimeThreshY;
LTADispData.Y_pdfXVals = pdfXVals;
LTADispData.Y_pdfYVals = pdfYVals*40;

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
% h(13) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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
% h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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

function [ETATimeData, ETADispData] = plotEMGTriggeredAvg_FS
motionEventsLocationsX = [];
motionEventsLocationsY = [];
timeToThreshX = [];
timeToThreshY = [];
dispTimeThreshX = [];
dispTimeThreshY = [];
moveThresh = .75;
timeThresh = 1.5;
binEdgesTime = -2:.25:3;
binEdgesDisp = -4:.25:4;
load('ETADataCell_FS.mat')
for n = 1:size(EMGDataCell,1)
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

h(9) = figure('Color','White','Name','Figure 3g','NumberTitle','off');
subplot(2,2,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
f = fill([3 0 0 3],[2.9 2.9 -2.9 -2.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecX,meanX,'k')
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
% plot([0 0],[-3 3],'r')
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[1,0,0,0.1])
end
% plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean Motion During EMG Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Brian Shift (\mum)')
ylim([-3 3])
xlim([-2 3])
grid on

subplot(2,2,3)
f = fill([3 0 0 3],[2.9 2.9 -2.9 -2.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecY,meanY,'k')
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-3 3],'r')
for n = 1:size(motionEventsLocationsY,1)
    plot(timeVecY,-1*motionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
end
% plot(timeVecY(brainMotionStart),meanY(brainMotionStart),'rx')
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
for n = 1:size(EMGDataCell,1)
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
f = fill([-2 0 0 -2],[2.9 2.9 -2.9 -2.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecX,meanX,'k')
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-3 3],'r')
for n = 1:size(stopMotionEventsLocationsX,1)
    plot(timeVecX,stopMotionEventsLocationsX(n,:),'Color',[1,0,0,0.1])
end
hold off
text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% title(['\fontsize{20pt}\bf{Mean Motion During Stopping EMG Events, n = ' num2str(size(stopMotionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Brian Shift (\mum)')
ylim([-3 3])
xlim([-2 3])
grid on

subplot(2,2,4)
f = fill([-2 0 0 -2],[2.9 2.9 -2.9 -2.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecY,meanY,'k')
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-3 3],'r')
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

% h(10) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
% histogram(timeToThreshX,binEdgesTime);
% hold on
[pdfXVals,pdfYVals] = findKernelPDF(timeToThreshX,binEdgesTime);
% plot(pdfXVals,pdfYVals*25,'r','LineWidth',2)
% title(['Time for brain to displace laterally ' num2str(moveThresh) ' micrometers following EMG trigger'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 60])
% mu = mean(timeToThreshX);
% sig = std(timeToThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
ETATimeData.timeToThreshX = timeToThreshX;
ETATimeData.binEdgesTime = binEdgesTime;
ETATimeData.X_pdfXVals = pdfXVals;
ETATimeData.X_pdfYVals = pdfYVals*25;
% 
% subplot(2,1,2)
% histogram(timeToThreshY,binEdgesDisp);
% hold on
[pdfXVals,pdfYVals] = findKernelPDF(timeToThreshY,binEdgesTime);
% plot(pdfXVals,pdfYVals*25,'r','LineWidth',2)
% title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following EMG trigger'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 60])
% mu = mean(timeToThreshY);
% sig = std(timeToThreshY);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
ETATimeData.timeToThreshY = timeToThreshY;
ETATimeData.Y_pdfXVals = pdfXVals;
ETATimeData.Y_pdfYVals = pdfYVals*25;
% 
% h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
% histogram(dispTimeThreshX,binEdgesDisp);
% hold on
[pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshX,binEdgesDisp);
% plot(pdfXVals,pdfYVals*35,'r','LineWidth',2)
% title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following EMG trigger'])
% xlabel('Displacement (\mum)')
% xlim([-4 4])
% ylim([0 60])
% mu = mean(dispTimeThreshX);
% sig = std(dispTimeThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
ETADispData.dispTimeThreshX = dispTimeThreshX;
ETADispData.binEdgesDisp = binEdgesDisp;
ETADispData.X_pdfXVals = pdfXVals;
ETADispData.X_pdfYVals = pdfYVals*35;

% subplot(2,1,2)
% histogram(dispTimeThreshY,binEdgesDisp);
% hold on
[pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshY,binEdgesDisp);
% plot(pdfXVals,pdfYVals*35,'r','LineWidth',2)
% title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following EMG trigger'])
% xlabel('Displacement (\mum)')
% xlim([-4 4])
% ylim([0 60])
% mu = mean(dispTimeThreshY);
% sig = std(dispTimeThreshY);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
ETADispData.dispTimeThreshY = dispTimeThreshY;
ETADispData.Y_pdfXVals = pdfXVals;
ETADispData.Y_pdfYVals = pdfYVals*35;

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
% h(13) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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
% h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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

function plotLocoEMGHist_FS(LTATimeData,LTADispData,ETATimeData,ETADispData)

h(10) = figure('Color','White','Name','Supplementary Figure 10a','NumberTitle','off');
subplot(1,2,1)
histogram(LTATimeData.timeToThreshX,LTATimeData.binEdgesTime,'FaceColor',[.25 .25 .25],'FaceAlpha',0.5);
hold on
histogram(ETATimeData.timeToThreshX,ETATimeData.binEdgesTime,'FaceColor',[1 .5 0],'FaceAlpha',0.5);
hold off
% [pdfXVals,pdfYVals] = findKernelPDF(LTATimeData.timeToThreshX,LTATimeData.binEdgesTime);
% plot(LTATimeData.pdfXVals,LTATimeData.pdfYVals,'r','LineWidth',2)
title('Time for brain to displace laterally 0.75 micrometers following trigger')
xlabel('Time (s)')
xlim([-2 3])
ylim([0 60])
text(-2,55,'Locomotion Trigger','Color',[.25 .25 .25],'FontSize',12)
text(-2,50,'EMG Trigger','Color',[1 .5 0],'FontSize',12)
% mu = mean(LTATimeData.timeToThreshX);
% sig = std(LTATimeData.timeToThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% LTATimeData.timeToThreshX = timeToThreshX;
% LTATimeData.binEdgesTime = binEdgesTime;
% LTATimeData.X_pdfXVals = pdfXVals;
% LTATimeData.X_pdfYVals = pdfYVals*20;

subplot(1,2,2)
plot(LTATimeData.X_pdfXVals,LTATimeData.X_pdfYVals,'Color',[.25 .25 .25],'LineWidth',2)
hold on
plot(ETATimeData.X_pdfXVals,ETATimeData.X_pdfYVals,'Color',[1 .5 0],'LineWidth',2)
% [pdfXVals,pdfYVals] = findKernelPDF(timeToThreshX,binEdgesTime);
% plot(pdfXVals,pdfYVals*25,'r','LineWidth',2)
title(['Time for brain to displace laterally 0.75 micrometers following trigger'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 60])
text(-2,55,'Locomotion Trigger','Color',[.25 .25 .25],'FontSize',12)
text(-2,50,'EMG Trigger','Color',[1 .5 0],'FontSize',12)
mu = mean(LTATimeData.timeToThreshX);
sig = std(LTATimeData.timeToThreshX);
plot([mu mu],[0 60],'Color',[.25 .25 .25],'LineWidth',2);
plot([mu+sig mu+sig],[0 60],'Color',[.25 .25 .25],'LineStyle','--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'Color',[.25 .25 .25],'LineStyle','--','LineWidth',2);
text(-2,15,['displacement thresh = 0.75, n = ' num2str(length(LTATimeData.timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)],'Color',[.25 .25 .25])
mu = mean(ETATimeData.timeToThreshX);
sig = std(ETATimeData.timeToThreshX);
plot([mu mu],[0 60],'Color',[1 .5 0],'LineWidth',2);
plot([mu+sig mu+sig],[0 60],'Color',[1 .5 0],'LineStyle','--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'Color',[1 .5 0],'LineStyle','--','LineWidth',2);
text(-2,10,['displacement thresh = 0.75, n = ' num2str(length(ETATimeData.timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)],'Color',[1 .5 0])
hold off
% ETATimeData.timeToThreshX = timeToThreshX;
% ETATimeData.binEdgesTime = binEdgesTime;
% ETATimeData.X_pdfXVals = pdfXVals;
% ETATimeData.X_pdfYVals = pdfYVals*25;

h(10) = figure('Color','White','Name','Supplementary Figure 10b','NumberTitle','off');
subplot(1,2,1)
histogram(LTATimeData.timeToThreshY,LTATimeData.binEdgesTime,'FaceColor',[.25 .25 .25],'FaceAlpha',0.5);
hold on
histogram(ETATimeData.timeToThreshY,ETATimeData.binEdgesTime,'FaceColor',[1 .5 0],'FaceAlpha',0.5);
hold off
% [pdfXVals,pdfYVals] = findKernelPDF(LTATimeData.timeToThreshX,LTATimeData.binEdgesTime);
% plot(LTATimeData.pdfXVals,LTATimeData.pdfYVals,'r','LineWidth',2)
title('Time for brain to displace rostrally 0.75 micrometers following trigger')
xlabel('Time (s)')
xlim([-2 3])
ylim([0 60])
text(-2,55,'Locomotion Trigger','Color',[.25 .25 .25],'FontSize',12)
text(-2,50,'EMG Trigger','Color',[1 .5 0],'FontSize',12)
% mu = mean(LTATimeData.timeToThreshX);
% sig = std(LTATimeData.timeToThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% LTATimeData.timeToThreshX = timeToThreshX;
% LTATimeData.binEdgesTime = binEdgesTime;
% LTATimeData.X_pdfXVals = pdfXVals;
% LTATimeData.X_pdfYVals = pdfYVals*20;

subplot(1,2,2)
plot(LTATimeData.Y_pdfXVals,LTATimeData.Y_pdfYVals,'Color',[.25 .25 .25],'LineWidth',2)
hold on
plot(ETATimeData.Y_pdfXVals,ETATimeData.Y_pdfYVals,'Color',[1 .5 0],'LineWidth',2)
% [pdfXVals,pdfYVals] = findKernelPDF(timeToThreshX,binEdgesTime);
% plot(pdfXVals,pdfYVals*25,'r','LineWidth',2)
title(['Time for brain to displace rostrally 0.75 micrometers following trigger'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 60])
text(-2,55,'Locomotion Trigger','Color',[.25 .25 .25],'FontSize',12)
text(-2,50,'EMG Trigger','Color',[1 .5 0],'FontSize',12)
mu = mean(LTATimeData.timeToThreshY);
sig = std(LTATimeData.timeToThreshY);
plot([mu mu],[0 60],'Color',[.25 .25 .25],'LineWidth',2);
plot([mu+sig mu+sig],[0 60],'Color',[.25 .25 .25],'LineStyle','--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'Color',[.25 .25 .25],'LineStyle','--','LineWidth',2);
text(-2,15,['displacement thresh = 0.75, n = ' num2str(length(LTATimeData.timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)],'Color',[.25 .25 .25])
mu = mean(ETATimeData.timeToThreshY);
sig = std(ETATimeData.timeToThreshY);
plot([mu mu],[0 60],'Color',[1 .5 0],'LineWidth',2);
plot([mu+sig mu+sig],[0 60],'Color',[1 .5 0],'LineStyle','--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'Color',[1 .5 0],'LineStyle','--','LineWidth',2);
text(-2,10,['displacement thresh = 0.75, n = ' num2str(length(ETATimeData.timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)],'Color',[1 .5 0])
hold off
% ETATimeData.timeToThreshX = timeToThreshX;
% ETATimeData.binEdgesTime = binEdgesTime;
% ETATimeData.X_pdfXVals = pdfXVals;
% ETATimeData.X_pdfYVals = pdfYVals*25;

h(50) = figure('Color','White','Name','Supplementary Figure 10c','NumberTitle','off');
subplot(1,2,1)
histogram(LTADispData.dispTimeThreshX,LTADispData.binEdgesDisp,'FaceColor',[.25 .25 .25],'FaceAlpha',0.5);
hold on
histogram(ETADispData.dispTimeThreshX,ETADispData.binEdgesDisp,'FaceColor',[1 .5 0],'FaceAlpha',0.5);
% [pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshX,binEdgesDisp);
% plot(pdfXVals,pdfYVals*40,'r','LineWidth',2)
title('Lateral displacement of brain after 1.5 seconds following trigger')
xlabel('Displacement (\mum)')
xlim([-4 4])
ylim([0 60])
text(-4,55,'Locomotion Trigger','Color',[.25 .25 .25],'FontSize',12)
text(-4,50,'EMG Trigger','Color',[1 .5 0],'FontSize',12)
% mu = mean(dispTimeThreshX);
% sig = std(dispTimeThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(1,2,2)
plot(LTADispData.X_pdfXVals,LTADispData.X_pdfYVals,'Color',[.25 .25 .25],'LineWidth',2)
hold on
plot(ETADispData.X_pdfXVals,ETADispData.X_pdfYVals,'Color',[1 .5 0],'LineWidth',2)
% [pdfXVals,pdfYVals] = findKernelPDF(timeToThreshX,binEdgesTime);
% plot(pdfXVals,pdfYVals*25,'r','LineWidth',2)
title('Lateral displacement of brain after 1.5 seconds following trigger')
xlabel('Displacement (\mum)')
xlim([-4 4])
ylim([0 60])
text(-4,55,'Locomotion Trigger','Color',[.25 .25 .25],'FontSize',12)
text(-4,50,'EMG Trigger','Color',[1 .5 0],'FontSize',12)
mu = mean(LTADispData.dispTimeThreshX);
sig = std(LTADispData.dispTimeThreshX);
plot([mu mu],[0 60],'Color',[.25 .25 .25],'LineWidth',2);
plot([mu+sig mu+sig],[0 60],'Color',[.25 .25 .25],'LineStyle','--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'Color',[.25 .25 .25],'LineStyle','--','LineWidth',2);
text(-2,15,['displacement thresh = 0.75, n = ' num2str(length(LTADispData.dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)],'Color',[.25 .25 .25])
mu = mean(ETADispData.dispTimeThreshX);
sig = std(ETADispData.dispTimeThreshX);
plot([mu mu],[0 60],'Color',[1 .5 0],'LineWidth',2);
plot([mu+sig mu+sig],[0 60],'Color',[1 .5 0],'LineStyle','--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'Color',[1 .5 0],'LineStyle','--','LineWidth',2);
text(-2,10,['displacement thresh = 0.75, n = ' num2str(length(ETADispData.dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)],'Color',[1 .5 0])
hold off

h(50) = figure('Color','White','Name','Supplementary Figure 10d','NumberTitle','off');
subplot(1,2,1)
histogram(LTADispData.dispTimeThreshY,LTADispData.binEdgesDisp,'FaceColor',[.25 .25 .25],'FaceAlpha',0.5);
hold on
histogram(ETADispData.dispTimeThreshY,ETADispData.binEdgesDisp,'FaceColor',[1 .5 0],'FaceAlpha',0.5);
% [pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshX,binEdgesDisp);
% plot(pdfXVals,pdfYVals*40,'r','LineWidth',2)
title('Rostral displacement of brain after 1.5 seconds following trigger')
xlabel('Displacement (\mum)')
xlim([-4 4])
ylim([0 60])
text(-4,55,'Locomotion Trigger','Color',[.25 .25 .25],'FontSize',12)
text(-4,50,'EMG Trigger','Color',[1 .5 0],'FontSize',12)
% mu = mean(dispTimeThreshX);
% sig = std(dispTimeThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(1,2,2)
plot(LTADispData.Y_pdfXVals,LTADispData.Y_pdfYVals,'Color',[.25 .25 .25],'LineWidth',2)
hold on
plot(ETADispData.Y_pdfXVals,ETADispData.Y_pdfYVals,'Color',[1 .5 0],'LineWidth',2)
% [pdfXVals,pdfYVals] = findKernelPDF(timeToThreshX,binEdgesTime);
% plot(pdfXVals,pdfYVals*25,'r','LineWidth',2)
title('Rostral displacement of brain after 1.5 seconds following trigger')
xlabel('Displacement (\mum)')
xlim([-4 4])
ylim([0 60])
text(-4,55,'Locomotion Trigger','Color',[.25 .25 .25],'FontSize',12)
text(-4,50,'EMG Trigger','Color',[1 .5 0],'FontSize',12)
mu = mean(LTADispData.dispTimeThreshY);
sig = std(LTADispData.dispTimeThreshY);
plot([mu mu],[0 60],'Color',[.25 .25 .25],'LineWidth',2);
plot([mu+sig mu+sig],[0 60],'Color',[.25 .25 .25],'LineStyle','--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'Color',[.25 .25 .25],'LineStyle','--','LineWidth',2);
text(-2,15,['displacement thresh = 0.75, n = ' num2str(length(LTADispData.dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)],'Color',[.25 .25 .25])
mu = mean(ETADispData.dispTimeThreshY);
sig = std(ETADispData.dispTimeThreshY);
plot([mu mu],[0 60],'Color',[1 .5 0],'LineWidth',2);
plot([mu+sig mu+sig],[0 60],'Color',[1 .5 0],'LineStyle','--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'Color',[1 .5 0],'LineStyle','--','LineWidth',2);
text(-2,10,['displacement thresh = 0.75, n = ' num2str(length(ETADispData.dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)],'Color',[1 .5 0])
hold off

% histogram(timeToThreshY,binEdgesTime);
% hold on
% [pdfXVals,pdfYVals] = findKernelPDF(timeToThreshY,binEdgesTime);
% plot(pdfXVals,pdfYVals*40,'r','LineWidth',2)
% title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following locomotion trigger'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 60])
% mu = mean(timeToThreshY);
% sig = std(timeToThreshY);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% % LTATimeData.timeToThreshY = timeToThreshY;
% % LTATimeData.Y_pdfXVals = pdfXVals;
% % LTATimeData.Y_pdfYVals = pdfYVals*40;
% 
% h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
% histogram(dispTimeThreshX,binEdgesDisp);
% hold on
% [pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshX,binEdgesDisp);
% plot(pdfXVals,pdfYVals*40,'r','LineWidth',2)
% title(['Lateral displacement of brain after 1.5 seconds following locomotion trigger'])
% xlabel('Displacement (\mum)')
% xlim([-4 4])
% ylim([0 60])
% mu = mean(dispTimeThreshX);
% sig = std(dispTimeThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% % LTADispData.dispTimeThreshX = dispTimeThreshX;
% % LTADispData.binEdgesDisp = binEdgesDisp;
% % LTADispData.X_pdfXVals = pdfXVals;
% % LTADispData.X_pdfYVals = pdfYVals*40;
% 
% subplot(2,1,2)
% histogram(dispTimeThreshY,binEdgesDisp);
% hold on
% [pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshY,binEdgesDisp);
% plot(pdfXVals,pdfYVals*40,'r','LineWidth',2)
% title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following locomotion trigger'])
% xlabel('Displacement (\mum)')
% xlim([-4 4])
% ylim([0 60])
% mu = mean(dispTimeThreshY);
% sig = std(dispTimeThreshY);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% % LTADispData.dispTimeThreshY = dispTimeThreshY;
% % LTADispData.Y_pdfXVals = pdfXVals;
% % LTADispData.Y_pdfYVals = pdfYVals*40;
% 
% subplot(2,1,2)
% histogram(timeToThreshY,binEdgesDisp);
% hold on
% [pdfXVals,pdfYVals] = findKernelPDF(timeToThreshY,binEdgesTime);
% plot(pdfXVals,pdfYVals*25,'r','LineWidth',2)
% title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following EMG trigger'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 60])
% mu = mean(timeToThreshY);
% sig = std(timeToThreshY);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% % ETATimeData.timeToThreshY = timeToThreshY;
% % ETATimeData.Y_pdfXVals = pdfXVals;
% % ETATimeData.Y_pdfYVals = pdfYVals*25;
% 
% h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
% histogram(dispTimeThreshX,binEdgesDisp);
% hold on
% [pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshX,binEdgesDisp);
% plot(pdfXVals,pdfYVals*35,'r','LineWidth',2)
% title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following EMG trigger'])
% xlabel('Displacement (\mum)')
% xlim([-4 4])
% ylim([0 60])
% mu = mean(dispTimeThreshX);
% sig = std(dispTimeThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% % ETADispData.dispTimeThreshX = dispTimeThreshX;
% % ETADispData.binEdgesDisp = binEdgesDisp;
% % ETADispData.X_pdfXVals = pdfXVals;
% % ETADispData.X_pdfYVals = pdfYVals*35;
% 
% subplot(2,1,2)
% histogram(dispTimeThreshY,binEdgesDisp);
% hold on
% [pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshY,binEdgesDisp);
% plot(pdfXVals,pdfYVals*35,'r','LineWidth',2)
% title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following EMG trigger'])
% xlabel('Displacement (\mum)')
% xlim([-4 4])
% ylim([0 60])
% mu = mean(dispTimeThreshY);
% sig = std(dispTimeThreshY);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% % ETADispData.dispTimeThreshY = dispTimeThreshY;
% % ETADispData.Y_pdfXVals = pdfXVals;
% % ETADispData.Y_pdfYVals = pdfYVals*35;
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
for n = 1:size(locDataCell,1)
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

h(9) = figure('Color','White','Name','Supplementary Figure 6c','NumberTitle','off');
subplot(2,2,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
f = fill([3 0 0 3],[2.9 2.9 -2.9 -2.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecX,meanX,'k')
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-3 3],'r')
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[1,0,0,0.1])
end
% plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean Skull Motion During Locomotion Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Skull Shift (\mum)')
ylim([-3 3])
xlim([-2 3])
grid on

subplot(2,2,3)
f = fill([3 0 0 3],[2.9 2.9 -2.9 -2.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecY,meanY,'k')
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-3 3],'r')
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
for n = 1:size(locDataCell,1)
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
f = fill([-2 0 0 -2],[2.9 2.9 -2.9 -2.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecX,meanX,'k')
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-3 3],'r')
for n = 1:size(stopMotionEventsLocationsX,1)
    plot(timeVecX,stopMotionEventsLocationsX(n,:),'Color',[1,0,0,0.1])
end
hold off
text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% title(['\fontsize{20pt}\bf{Mean Skull Motion During Stopping Locomotion Events, n = ' num2str(size(stopMotionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Skull Shift (\mum)')
ylim([-3 3])
xlim([-2 3])
grid on

subplot(2,2,4)
f = fill([-2 0 0 -2],[2.9 2.9 -2.9 -2.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecY,meanY,'k')
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-3 3],'r')
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

% h(10) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
% histfit(timeToThreshX,numBins,'kernel')
% title(['Time for skull to displace laterally ' num2str(moveThresh) ' micrometers following locomotion trigger'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 30])
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
% title(['Time for skull to displace rostrally ' num2str(moveThresh) ' micrometers following locomotion trigger'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 30])
% hold on
% mu = mean(timeToThreshY);
% sig = std(timeToThreshY);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
% histfit(dispTimeThreshX,numBins,'kernel')
% title(['Lateral displacement of skull after ' num2str(timeThresh) ' s following locomotion trigger'])
% xlabel('Displacement (\mum)')
% xlim([-2 3])
% ylim([0 60])
% hold on
% mu = mean(dispTimeThreshX);
% sig = std(dispTimeThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% subplot(2,1,2)
% histfit(dispTimeThreshY,numBins,'kernel')
% title(['Rostral displacement of skull after ' num2str(timeThresh) ' s following locomotion trigger'])
% xlabel('Displacement (\mum)')
% xlim([-2 3])
% ylim([0 60])
% hold on
% mu = mean(dispTimeThreshY);
% sig = std(dispTimeThreshY);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
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
% h(13) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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
% h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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
for n = 1:size(EMGDataCell,1)
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

h(9) = figure('Color','White','Name','Supplementary Figure 6d','NumberTitle','off');
subplot(2,2,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
f = fill([3 0 0 3],[2.9 2.9 -2.9 -2.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecX,meanX,'k')
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-3 3],'r')
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[1,0,0,0.1])
end
% plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean Skull Motion During EMG Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Skull Shift (\mum)')
ylim([-3 3])
xlim([-2 3])
grid on

subplot(2,2,3)
f = fill([3 0 0 3],[2.9 2.9 -2.9 -2.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecY,meanY,'k')
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-3 3],'r')
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
for n = 1:size(EMGDataCell,1)
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
f = fill([-2 0 0 -2],[2.9 2.9 -2.9 -2.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecX,meanX,'k')
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-3 3],'r')
for n = 1:size(stopMotionEventsLocationsX,1)
    plot(timeVecX,stopMotionEventsLocationsX(n,:),'Color',[1,0,0,0.1])
end
hold off
text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% title(['\fontsize{20pt}\bf{Mean Skull Motion During Stopping EMG Events, n = ' num2str(size(stopMotionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Skull Shift (\mum)')
ylim([-3 3])
xlim([-2 3])
grid on

subplot(2,2,4)
f = fill([-2 0 0 -2],[2.9 2.9 -2.9 -2.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecY,meanY,'k')
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-3 3],'r')
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

% h(10) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
% histfit(timeToThreshX,numBins,'kernel')
% title(['Time for skull to displace laterally ' num2str(moveThresh) ' micrometers following EMG trigger'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 30])
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
% title(['Time for skull to displace rostrally ' num2str(moveThresh) ' micrometers following EMG trigger'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 30])
% hold on
% mu = mean(timeToThreshY);
% sig = std(timeToThreshY);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
% histfit(dispTimeThreshX,numBins,'kernel')
% title(['Lateral displacement of skull after ' num2str(timeThresh) ' s following EMG trigger'])
% xlabel('Displacement (\mum)')
% xlim([-3 4])
% ylim([0 60])
% hold on
% mu = mean(dispTimeThreshX);
% sig = std(dispTimeThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% subplot(2,1,2)
% histfit(dispTimeThreshY,numBins,'kernel')
% title(['Rostral displacement of skull after ' num2str(timeThresh) ' s following EMG trigger'])
% xlabel('Displacement (\mum)')
% xlim([-3 4])
% ylim([0 60])
% hold on
% mu = mean(dispTimeThreshY);
% sig = std(dispTimeThreshY);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

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
% h(13) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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
% h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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

function plotLocomotionTriggeredAvgOlf_FS
motionEventsLocationsX = [];
motionEventsLocationsY = [];
timeToThreshX = [];
timeToThreshY = [];
dispTimeThreshX = [];
dispTimeThreshY = [];
moveThresh = .75;
timeThresh = 1.5;
binEdgesTime = -2:.25:3;
binEdgesDisp = -4:.25:4;
load('LTADataCellOlf_FS.mat')
for n = 1:size(locDataCellOlf,1)
    if isnan(locDataCellOlf{n,3})
        continue
    end
    motionVectorX = locDataCellOlf{n,3}(2,:);
    motionVectorY = locDataCellOlf{n,4}(2,:);
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
    
    singleTimeVecX = linspace(round(locDataCellOlf{1,2}(1,1)-locDataCellOlf{1,2}(1,2)),round(locDataCellOlf{1,2}(1,3)-locDataCellOlf{1,2}(1,2)),length(motionVectorX));
    dispTimeThreshX(end+1) = motionVectorX((find(singleTimeVecX>timeThresh,1))) - motionVectorX((find(singleTimeVecX>0,1)));
    idxToThreshXSingle = find((motionVectorX - motionVectorX(find(singleTimeVecX>0,1))) > moveThresh & singleTimeVecX>0 & singleTimeVecX<=3,1);
    if ~isempty(idxToThreshXSingle)
        timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
    end
    singleTimeVecY = linspace(round(locDataCellOlf{1,2}(1,1)-locDataCellOlf{1,2}(1,2)),round(locDataCellOlf{1,2}(1,3)-locDataCellOlf{1,2}(1,2)),length(motionVectorY));
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
timeVecX = linspace(round(locDataCellOlf{1,2}(1,1)-locDataCellOlf{1,2}(1,2)),round(locDataCellOlf{1,2}(1,3)-locDataCellOlf{1,2}(1,2)),length(meanX));
timeVecY = linspace(round(locDataCellOlf{1,2}(1,1)-locDataCellOlf{1,2}(1,2)),round(locDataCellOlf{1,2}(1,3)-locDataCellOlf{1,2}(1,2)),length(meanY));

stdWindowSize = 5;
for n = stdWindowSize+1:length(meanY)
    if std(meanY(n-stdWindowSize:n)) > .01
        brainMotionStart = n;
        break
    end 
end

h(9) = figure('Color','White','Name','Supplementary Figure 14b','NumberTitle','off');
subplot(2,2,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
f = fill([3 0 0 3],[2.9 2.9 -2.9 -2.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecX,meanX,'k')
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-3 3],'r')
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[1,0,0,0.1])
end
% plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean Olfactory Motion During Locomotion Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Brian Shift (\mum)')
ylim([-3 3])
xlim([-2 3])
grid on

subplot(2,2,3)
f = fill([3 0 0 3],[2.9 2.9 -2.9 -2.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecY,meanY,'k')
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-3 3],'r')
for n = 1:size(motionEventsLocationsY,1)
    plot(timeVecY,-1*motionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
end
% plot(timeVecY(brainMotionStart),meanY(brainMotionStart),'rx')
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
for n = 1:size(locDataCellOlf,1)
    if isnan(locDataCellOlf{n,6})
        continue
    end
    stopMotionVectorX = locDataCellOlf{n,6}(2,:);
    stopMotionVectorY = locDataCellOlf{n,7}(2,:);
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
timeVecX = linspace(round(locDataCellOlf{1,5}(1,1)-locDataCellOlf{1,5}(1,2)),round(locDataCellOlf{1,5}(1,3)-locDataCellOlf{1,5}(1,2)),length(meanX));
timeVecY = linspace(round(locDataCellOlf{1,5}(1,1)-locDataCellOlf{1,5}(1,2)),round(locDataCellOlf{1,5}(1,3)-locDataCellOlf{1,5}(1,2)),length(meanY));

subplot(2,2,2)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
f = fill([-2 0 0 -2],[2.9 2.9 -2.9 -2.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecX,meanX,'k')
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-3 3],'r')
for n = 1:size(stopMotionEventsLocationsX,1)
    plot(timeVecX,stopMotionEventsLocationsX(n,:),'Color',[1,0,0,0.1])
end
hold off
text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% title(['\fontsize{20pt}\bf{Mean Motion During Stopping Locomotion Events, n = ' num2str(size(stopMotionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Brian Shift (\mum)')
ylim([-3 3])
xlim([-2 3])
grid on

subplot(2,2,4)
f = fill([-2 0 0 -2],[2.9 2.9 -2.9 -2.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecY,meanY,'k')
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[-3 3],'r')
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
% 
% h(10) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
% histogram(timeToThreshX,binEdgesTime);
% hold on
% [pdfXVals,pdfYVals] = findKernelPDF(timeToThreshX,binEdgesTime);
% plot(pdfXVals,pdfYVals*20,'r','LineWidth',2)
% title(['Time for brain to displace laterally ' num2str(moveThresh) ' micrometers following locomotion trigger'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 60])
% mu = mean(timeToThreshX);
% sig = std(timeToThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% subplot(2,1,2)
% histogram(timeToThreshY,binEdgesTime);
% hold on
% [pdfXVals,pdfYVals] = findKernelPDF(timeToThreshY,binEdgesTime);
% plot(pdfXVals,pdfYVals*40,'r','LineWidth',2)
% title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following locomotion trigger'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 60])
% mu = mean(timeToThreshY);
% sig = std(timeToThreshY);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
% histogram(dispTimeThreshX,binEdgesDisp);
% hold on
% [pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshX,binEdgesDisp);
% plot(pdfXVals,pdfYVals*40,'r','LineWidth',2)
% title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following locomotion trigger'])
% xlabel('Displacement (\mum)')
% xlim([-4 4])
% ylim([0 60])
% mu = mean(dispTimeThreshX);
% sig = std(dispTimeThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% subplot(2,1,2)
% histogram(dispTimeThreshY,binEdgesDisp);
% hold on
% [pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshY,binEdgesDisp);
% plot(pdfXVals,pdfYVals*40,'r','LineWidth',2)
% title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following locomotion trigger'])
% xlabel('Displacement (\mum)')
% xlim([-4 4])
% ylim([0 60])
% mu = mean(dispTimeThreshY);
% sig = std(dispTimeThreshY);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
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
% h(13) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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
% h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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

function plotMotionTrackingBrainAndSkullOlf_FS(movementData,stationaryData)
[movementData.targetPosition,stationaryData.targetPosition] = interpBrainSkullMovement_FS(movementData,stationaryData);
if movementData.hemisphere == 2
    movementData.targetPosition(:,1) = movementData.targetPosition(:,1)*-1;
    stationaryData.targetPosition(:,1) = stationaryData.targetPosition(:,1)*-1;
end
movementData.targetPosition(:,2) = movementData.targetPosition(:,2)*-1;
stationaryData.targetPosition(:,2) = stationaryData.targetPosition(:,2)*-1;
% movementData.targetPosition(:,1) = movementData.targetPosition(:,1)-mean(movementData.targetPosition(2480:2530,1));
% movementData.targetPosition(:,2) = movementData.targetPosition(:,2)-mean(movementData.targetPosition(2480:2530,2));
% stationaryData.targetPosition(:,1) = stationaryData.targetPosition(:,1)-mean(stationaryData.targetPosition(2480:2530,1));
% stationaryData.targetPosition(:,2) = stationaryData.targetPosition(:,2)-mean(stationaryData.targetPosition(2480:2530,2));
movementData.secondsPerFrame = movementData.secondsPerFrame/2;
h(6) = figure('Color','White','Name','Supplementary Figure 14a','NumberTitle','off');
subplot(3,1,1)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
hold off
% title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
xlabel('Time (s)')
ylabel('Medial-Lateral Shift (\mum)')
grid on
axis([295 595 -5 5])
text(295,5,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(295,-5,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(3,1,2)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
hold off
xlabel('Time (s)')
ylabel('Rostral-Caudal Shift (\mum)')
grid on
axis([295 595 -5 5])
text(295,5,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(295,-5,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(4,1,3)
% plot(movementData.emgData(:,1),movementData.emgData(:,2),'k')
% xlabel('Time (s)')
% ylabel('Abdominal EMG (au)')
% grid on
% axis([62 272 0.5 3])
subplot(3,1,3)
plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)*2*pi*.06),'k')
xlabel('Time (s)')
ylabel('Treadmill Velocity (m/s)')
grid on
axis([295 595 0 0.3])

% h(6) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(4,1,1)
% plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
% hold on
% plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
% hold off
% title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
% xlabel('Time (s)')
% ylabel('X Position (\mum)')
% grid on
% axis([66 76 -6 6])
% text(62,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(62,-6,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(4,1,2)
% plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
% hold on
% plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
% hold off
% xlabel('Time (s)')
% ylabel('Y Position (\mum)')
% grid on
% axis([66 76 -6 6])
% text(62,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(62,-6,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(4,1,3)
% plot(movementData.emgData(:,1),movementData.emgData(:,2),'k')
% xlabel('Time (s)')
% ylabel('Abdominal EMG (au)')
% grid on
% axis([66 76 0.5 3])
% subplot(4,1,4)
% plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)*2*pi*.06),'r')
% xlabel('Time (s)')
% ylabel('Locomotion (m/s)')
% grid on
% axis([66 76 0 0.2])
% 
% h(6) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(4,1,1)
% plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
% hold on
% plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
% hold off
% title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
% xlabel('Time (s)')
% ylabel('X Position (\mum)')
% grid on
% axis([150 160 -6 6])
% text(62,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(62,-6,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(4,1,2)
% plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
% hold on
% plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
% hold off
% xlabel('Time (s)')
% ylabel('Y Position (\mum)')
% grid on
% axis([150 160 -6 6])
% text(62,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(62,-6,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(4,1,3)
% plot(movementData.emgData(:,1),movementData.emgData(:,2),'k')
% xlabel('Time (s)')
% ylabel('Abdominal EMG (au)')
% grid on
% axis([150 160 0.5 3])
% subplot(4,1,4)
% plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)*2*pi*.06),'r')
% xlabel('Time (s)')
% ylabel('Locomotion (m/s)')
% grid on
% axis([150 160 0 0.2])

% h(6) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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
% h(6) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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

function [STATimeData,STADispData] = plotSqueezeTriggeredAvg_FS
motionEventsLocationsX = [];
motionEventsLocationsY = [];
timeToThreshX = [];
timeToThreshY = [];
dispTimeThreshX = [];
dispTimeThreshY = [];
moveThresh = .75;
timeThresh = 1.5;
binEdgesTime = -2:.25:3;
binEdgesDisp = -4:.25:4;
load('squeezeDataCell_FS.mat')
for n = 1:size(squeezeDataCell,1)
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

h(9) = figure('Color','White','Name','Figure 5e','NumberTitle','off');
subplot(2,1,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
f = fill([2 0 0 2],[4.9 4.9 -.9 -.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecX,meanX,'k')
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[1,0,0,0.1])
end
% plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
text(6,-1,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(6,5,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean Motion During Squeeze Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Brian Shift (\mum)')
ylim([-1 5])
xlim([-2 6])
grid on

subplot(2,1,2)
f = fill([2 0 0 2],[4.9 4.9 -.9 -.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecY,meanY,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
for n = 1:size(motionEventsLocationsY,1)
    plot(timeVecY,-1*motionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
end
% plot(timeVecY(brainMotionStart),meanY(brainMotionStart),'rx')
hold off
text(6,-1,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(6,5,'Rostral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
xlabel('Time (s)')
ylabel('\Delta Brian Shift (\mum)')
ylim([-1 5])
xlim([-2 6])
grid on
clear movementData

% h(10) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
% histogram(timeToThreshX,binEdgesTime);
% hold on
[pdfXVals,pdfYVals] = findKernelPDF(timeToThreshX,binEdgesTime);
% plot(pdfXVals,pdfYVals*18,'r','LineWidth',2)
% title(['Time for brain to displace laterally ' num2str(moveThresh) ' micrometers following squeeze'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 60])
% mu = mean(timeToThreshX);
% sig = std(timeToThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
STATimeData.timeToThreshX = timeToThreshX;
STATimeData.binEdgesTime = binEdgesTime;
STATimeData.X_pdfXVals = pdfXVals;
STATimeData.X_pdfYVals = pdfYVals*18;
% 
% subplot(2,1,2)
% histogram(timeToThreshY,binEdgesTime);
% hold on
[pdfXVals,pdfYVals] = findKernelPDF(timeToThreshY,binEdgesTime);
% plot(pdfXVals,pdfYVals*30,'r','LineWidth',2)
% title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following squeeze'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 60])
% mu = mean(timeToThreshY);
% sig = std(timeToThreshY);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
STATimeData.timeToThreshY = timeToThreshY;
STATimeData.Y_pdfXVals = pdfXVals;
STATimeData.Y_pdfYVals = pdfYVals*30;
% 
% h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
% histogram(dispTimeThreshX,binEdgesDisp);
% hold on
[pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshX,binEdgesDisp);
% plot(pdfXVals,pdfYVals*32,'r','LineWidth',2)
% title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following squeeze'])
% xlabel('Displacement (\mum)')
% xlim([-4 4])
% ylim([0 60])
% mu = mean(dispTimeThreshX);
% sig = std(dispTimeThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
STADispData.dispTimeThreshX = dispTimeThreshX;
STADispData.binEdgesDisp = binEdgesDisp;
STADispData.X_pdfXVals = pdfXVals;
STADispData.X_pdfYVals = pdfYVals*32;
% 
% subplot(2,1,2)
% histogram(dispTimeThreshY,binEdgesDisp);
% hold on
[pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshY,binEdgesDisp);
% plot(pdfXVals,pdfYVals*32,'r','LineWidth',2)
% title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following squeeze'])
% xlabel('Displacement (\mum)')
% xlim([-4 4])
% ylim([0 60])
% mu = mean(dispTimeThreshY);
% sig = std(dispTimeThreshY);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
STADispData.dispTimeThreshY = dispTimeThreshY;
STADispData.Y_pdfXVals = pdfXVals;
STADispData.Y_pdfYVals = pdfYVals*32;

% timeToThreshX = [];
% timeToThreshY = [];
% dispTimeThreshX = [];
% dispTimeThreshY = [];
% for n = 1:size(squeezeDataCell)
%     if isnan(squeezeDataCell{n,3})
%         continue
%     end
%     motionVectorX = squeezeDataCell{n,3}(2,:);
%     motionVectorY = squeezeDataCell{n,4}(2,:);
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
%     singleTimeVecX = linspace(round(squeezeDataCell{1,2}(1,1)-squeezeDataCell{1,2}(1,2)),round(squeezeDataCell{1,2}(1,3)-squeezeDataCell{1,2}(1,2)),length(motionVectorX));
%     timeThreshIdx = find(singleTimeVecX>(singleTimeVecX(brainMotionStart)+timeThresh),1);
%     dispTimeThreshX(end+1) = motionVectorX(timeThreshIdx) - motionVectorX(brainMotionStart);
%     idxToThreshXSingle = find((motionVectorX - motionVectorX(brainMotionStart)) > moveThresh & singleTimeVecX>singleTimeVecX(brainMotionStart) & singleTimeVecX<=3,1);
%     if ~isempty(idxToThreshXSingle)
%         timeToThreshX(end+1) = singleTimeVecX(idxToThreshXSingle);
%     end
%     singleTimeVecY = linspace(round(squeezeDataCell{1,2}(1,1)-squeezeDataCell{1,2}(1,2)),round(squeezeDataCell{1,2}(1,3)-squeezeDataCell{1,2}(1,2)),length(motionVectorY));
%     dispTimeThreshY(end+1) = (motionVectorY(timeThreshIdx) - motionVectorY(brainMotionStart))*-1;
%     idxToThreshYSingle = find((motionVectorY - motionVectorY(brainMotionStart))*-1 > moveThresh & singleTimeVecY>singleTimeVecY(brainMotionStart) & singleTimeVecY<=3,1);
%     if ~isempty(idxToThreshYSingle)
%         timeToThreshY(end+1) = singleTimeVecY(idxToThreshYSingle);
%     end
% end
% h(13) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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
% h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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

function plotSqueezeRespHist_FS(STATimeData,STADispData,RTATimeData,RTADispData)
timeThresh = 1.5;
moveThresh = 0.75;
h(10) = figure('Color','White','Name','Supplementary Figure 11a','NumberTitle','off');
subplot(1,2,1)
histogram(STATimeData.timeToThreshX,STATimeData.binEdgesTime,'FaceColor',[1 0 0]);
hold on
% [pdfXVals,pdfYVals] = findKernelPDF(timeToThreshX,binEdgesTime);
plot(STATimeData.X_pdfXVals,STATimeData.X_pdfYVals,'r','LineWidth',2)
title(['Time for brain to displace laterally ' num2str(moveThresh) ' micrometers following squeeze'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 60])
mu = mean(STATimeData.timeToThreshX);
sig = std(STATimeData.timeToThreshX);
plot([mu mu],[0 60],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(STATimeData.timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(1,2,2)
histogram(STADispData.dispTimeThreshX,STADispData.binEdgesDisp,'FaceColor',[1 0 0]);
hold on
% [pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshX,binEdgesDisp);
plot(STADispData.X_pdfXVals,STADispData.X_pdfYVals,'r','LineWidth',2)
title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following squeeze'])
xlabel('Displacement (\mum)')
xlim([-4 4])
ylim([0 60])
mu = mean(STADispData.dispTimeThreshX);
sig = std(STADispData.dispTimeThreshX);
plot([mu mu],[0 60],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(STADispData.dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

h(50) = figure('Color','White','Name','Supplementary Figure 11b','NumberTitle','off');
subplot(1,2,1)
histogram(STATimeData.timeToThreshY,STATimeData.binEdgesTime,'FaceColor',[0 0 1]);
hold on
% [pdfXVals,pdfYVals] = findKernelPDF(timeToThreshY,binEdgesTime);
plot(STATimeData.Y_pdfXVals,STATimeData.Y_pdfYVals,'r','LineWidth',2)
title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following squeeze'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 60])
mu = mean(STATimeData.timeToThreshY);
sig = std(STATimeData.timeToThreshY);
plot([mu mu],[0 60],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(STATimeData.timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(1,2,2)
histogram(STADispData.dispTimeThreshY,STADispData.binEdgesDisp,'FaceColor',[0 0 1]);
hold on
% [pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshY,binEdgesDisp);
plot(STADispData.Y_pdfXVals,STADispData.Y_pdfYVals,'r','LineWidth',2)
title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following squeeze'])
xlabel('Displacement (\mum)')
xlim([-4 4])
ylim([0 60])
mu = mean(STADispData.dispTimeThreshY);
sig = std(STADispData.dispTimeThreshY);
plot([mu mu],[0 60],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(STADispData.dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

timeThresh = 0.4;
moveThresh = 0.25;
h(10) = figure('Color','White','Name','Supplementary Figure 11c','NumberTitle','off');
subplot(1,2,1)
histogram(RTATimeData.timeToThreshX,RTATimeData.binEdgesTime,'FaceColor',[1 0 0]);
hold on
% [pdfXVals,pdfYVals] = findKernelPDF(timeToThreshX,binEdgesTime);
plot(RTATimeData.X_pdfXVals,RTATimeData.X_pdfYVals,'r','LineWidth',2)
title(['Time for brain to displace laterally ' num2str(moveThresh) ' micrometers following squeeze'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 60])
mu = mean(RTATimeData.timeToThreshX);
sig = std(RTATimeData.timeToThreshX);
plot([mu mu],[0 60],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(RTATimeData.timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(1,2,2)
histogram(RTADispData.dispTimeThreshX,RTADispData.binEdgesDisp,'FaceColor',[1 0 0]);
hold on
% [pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshX,binEdgesDisp);
plot(RTADispData.X_pdfXVals,RTADispData.X_pdfYVals,'r','LineWidth',2)
title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following squeeze'])
xlabel('Displacement (\mum)')
xlim([-4 4])
ylim([0 60])
mu = mean(RTADispData.dispTimeThreshX);
sig = std(RTADispData.dispTimeThreshX);
plot([mu mu],[0 60],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(RTADispData.dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

h(50) = figure('Color','White','Name','Supplementary Figure 11d','NumberTitle','off');
subplot(1,2,1)
histogram(RTATimeData.timeToThreshY,RTATimeData.binEdgesTime,'FaceColor',[0 0 1]);
hold on
% [pdfXVals,pdfYVals] = findKernelPDF(timeToThreshY,binEdgesTime);
plot(RTATimeData.Y_pdfXVals,RTATimeData.Y_pdfYVals,'r','LineWidth',2)
title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following squeeze'])
xlabel('Time (s)')
xlim([-2 3])
ylim([0 60])
mu = mean(RTATimeData.timeToThreshY);
sig = std(RTATimeData.timeToThreshY);
plot([mu mu],[0 60],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
hold off
text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(RTATimeData.timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

subplot(1,2,2)
histogram(RTADispData.dispTimeThreshY,RTADispData.binEdgesDisp,'FaceColor',[0 0 1]);
hold on
% [pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshY,binEdgesDisp);
plot(RTADispData.Y_pdfXVals,RTADispData.Y_pdfYVals,'r','LineWidth',2)
title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following squeeze'])
xlabel('Displacement (\mum)')
xlim([-4 4])
ylim([0 60])
mu = mean(RTADispData.dispTimeThreshY);
sig = std(RTADispData.dispTimeThreshY);
plot([mu mu],[0 60],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
hold off
text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(RTADispData.dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
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
for n = 1:size(squeezeDataCell,1)
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

h(9) = figure('Color','White','Name','Figure 5f','NumberTitle','off');
subplot(2,1,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
f = fill([2 0 0 2],[4.9 4.9 -.9 -.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecX,meanX,'k')
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[1,0,0,0.1])
end
% plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
text(6,-1,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(6,5,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean Skull Motion During Squeeze Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Skull Shift (\mum)')
ylim([-1 5])
xlim([-2 6])
grid on

subplot(2,1,2)
f = fill([2 0 0 2],[4.9 4.9 -.9 -.9],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecY,meanY,'k')
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
for n = 1:size(motionEventsLocationsY,1)
    plot(timeVecY,-1*motionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
end
% plot(timeVecY(brainMotionStart),meanY(brainMotionStart),'rx')
hold off
text(6,-1,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(6,5,'Rostral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
xlabel('Time (s)')
ylabel('\Delta Skull Shift (\mum)')
ylim([-1 5])
xlim([-2 6])
grid on
clear movementData

% h(10) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
% histfit(timeToThreshX,numBins,'kernel')
% title(['Time for skull to displace laterally ' num2str(moveThresh) ' micrometers following squeeze'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 30])
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
% title(['Time for skull to displace rostrally ' num2str(moveThresh) ' micrometers following squeeze'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 30])
% hold on
% mu = mean(timeToThreshY);
% sig = std(timeToThreshY);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
% histfit(dispTimeThreshX,numBins,'kernel')
% title(['Lateral displacement of skull after ' num2str(timeThresh) ' s following squeeze'])
% xlabel('Displacement (\mum)')
% xlim([-2 3])
% ylim([0 60])
% hold on
% mu = mean(dispTimeThreshX);
% sig = std(dispTimeThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% subplot(2,1,2)
% histfit(dispTimeThreshY,numBins,'kernel')
% title(['Rostral displacement of skull after ' num2str(timeThresh) ' s following squeeze'])
% xlabel('Displacement (\mum)')
% xlim([-2 3])
% ylim([0 60])
% hold on
% mu = mean(dispTimeThreshY);
% sig = std(dispTimeThreshY);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [RTATimeData,RTADispData] = plotRespTriggeredAvg_FS
motionEventsLocationsX = [];
motionEventsLocationsY = [];
timeToThreshX = [];
timeToThreshY = [];
dispTimeThreshX = [];
dispTimeThreshY = [];
moveThresh = .25;
timeThresh = .4;
binEdgesTime = -.25:.05:1;
binEdgesDisp = -2.5:.1:2.5;
load('respDataCell_FS.mat')
for n = 1:size(respDataCell,1)
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

h(9) = figure('Color','White','Name','Supplementary Figure 4b','NumberTitle','off');
subplot(2,1,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
plot(timeVecX,meanX,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
plot([0 0],[-1 3],'r')
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[1,0,0,0.1])
end
% plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
text(5,-1,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(5,5,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean Motion During Respiration Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Brian Shift (\mum)')
ylim([-1 3])
xlim([-.25 1])
grid on

subplot(2,1,2)
plot(timeVecY,meanY,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,[.5 .5 .5],'Linestyle','none');
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
ylabel('\Delta Brian Shift (\mum)')
ylim([-1 3])
xlim([-.25 1])
grid on
clear movementData

% h(10) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
% histogram(timeToThreshX,binEdgesTime);
% hold on
[pdfXVals,pdfYVals] = findKernelPDF(timeToThreshX,binEdgesTime);
% plot(pdfXVals,pdfYVals*5,'r','LineWidth',2)
% title(['Time for brain to displace laterally ' num2str(moveThresh) ' micrometers following respiration'])
% xlabel('Time (s)')
% xlim([-.25 1])
% ylim([0 60])
% mu = mean(timeToThreshX);
% sig = std(timeToThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
RTATimeData.timeToThreshX = timeToThreshX;
RTATimeData.binEdgesTime = binEdgesTime;
RTATimeData.X_pdfXVals = pdfXVals;
RTATimeData.X_pdfYVals = pdfYVals*5;
% 
% subplot(2,1,2)
% histogram(timeToThreshY,binEdgesTime);
% hold on
[pdfXVals,pdfYVals] = findKernelPDF(timeToThreshY,binEdgesTime);
% plot(pdfXVals,pdfYVals*5,'r','LineWidth',2)
% title(['Time for brain to displace rostrally ' num2str(moveThresh) ' micrometers following respiration'])
% xlabel('Time (s)')
% xlim([-.25 1])
% ylim([0 60])
% mu = mean(timeToThreshY);
% sig = std(timeToThreshY);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
RTATimeData.timeToThreshY = timeToThreshY;
RTATimeData.Y_pdfXVals = pdfXVals;
RTATimeData.Y_pdfYVals = pdfYVals*5;
% 
% h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
% histogram(dispTimeThreshX,binEdgesDisp);
% hold on
[pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshX,binEdgesDisp);
% plot(pdfXVals,pdfYVals*12,'r','LineWidth',2)
% title(['Lateral displacement of brain after ' num2str(timeThresh) ' s following respiration'])
% xlabel('Displacement (\mum)')
% xlim([-2.5 2.5])
% ylim([0 60])
% mu = mean(dispTimeThreshX);
% sig = std(dispTimeThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
RTADispData.dispTimeThreshX = dispTimeThreshX;
RTADispData.binEdgesDisp = binEdgesDisp;
RTADispData.X_pdfXVals = pdfXVals;
RTADispData.X_pdfYVals = pdfYVals*12;
% 
% subplot(2,1,2)
% histogram(dispTimeThreshY,binEdgesDisp);
% hold on
[pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshY,binEdgesDisp);
% plot(pdfXVals,pdfYVals*12,'r','LineWidth',2)
% title(['Rostral displacement of brain after ' num2str(timeThresh) ' s following respiration'])
% xlabel('Displacement (\mum)')
% xlim([-2.5 2.5])
% ylim([0 60])
% mu = mean(dispTimeThreshY);
% sig = std(dispTimeThreshY);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
RTADispData.dispTimeThreshY = dispTimeThreshY;
RTADispData.Y_pdfXVals = pdfXVals;
RTADispData.Y_pdfYVals = pdfYVals*12;
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
for n = 1:size(respDataCell,1)
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

h(9) = figure('Color','White','Name','Supplementary Figure 4c','NumberTitle','off');
subplot(2,1,1)
maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
plot(timeVecX,meanX,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
plot([0 0],[-1 3],'r')
set(f,'facea',[.2]);
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[1,0,0,0.1])
end
% plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
text(5,-1,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(5,5,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean Skull Motion During Respiration Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('\Delta Skull Shift (\mum)')
ylim([-1 3])
xlim([-.25 1])
grid on

subplot(2,1,2)
plot(timeVecY,meanY,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,[.5 .5 .5],'Linestyle','none');
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

% h(10) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
% histfit(timeToThreshX,numBins,'kernel')
% title(['Time for skull to displace laterally ' num2str(moveThresh) ' micrometers following respiration'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 30])
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
% title(['Time for skull to displace rostrally ' num2str(moveThresh) ' micrometers following respiration'])
% xlabel('Time (s)')
% xlim([-2 3])
% ylim([0 30])
% hold on
% mu = mean(timeToThreshY);
% sig = std(timeToThreshY);
% plot([mu mu],[0 40],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 40],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 40],'k--','LineWidth',2);
% hold off
% text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
% histfit(dispTimeThreshX,numBins,'kernel')
% title(['Lateral displacement of skull after ' num2str(timeThresh) ' s following respiration'])
% xlabel('Displacement (\mum)')
% xlim([-2 3])
% ylim([0 60])
% hold on
% mu = mean(dispTimeThreshX);
% sig = std(dispTimeThreshX);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
% 
% subplot(2,1,2)
% histfit(dispTimeThreshY,numBins,'kernel')
% title(['Rostral displacement of skull after ' num2str(timeThresh) ' s following respiration'])
% xlabel('Displacement (\mum)')
% xlim([-2 3])
% ylim([0 60])
% hold on
% mu = mean(dispTimeThreshY);
% sig = std(dispTimeThreshY);
% plot([mu mu],[0 60],'k','LineWidth',2);
% plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% hold off
% text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshY)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotLocomotionTriggeredAvgEMG_FS
motionEventsLocationsX = [];
% motionEventsLocationsY = [];
timeToThreshX = [];
% timeToThreshY = [];
dispTimeThreshX = [];
% dispTimeThreshY = [];
moveThresh = 1.5;
timeThresh = -.25;
binEdgesTime = -2:.1:3;
binEdgesDisp = 0:.1:4;
load('LTADataCellEMG_FS.mat')
for n = 1:size(locDataCell,1)
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
%     dispTimeThreshX(end+1) = motionVectorX((find(singleTimeVecX>timeThresh,1))) - motionVectorX((find(singleTimeVecX>0,1)));
    dispTimeThreshX(end+1) = motionVectorX((find(singleTimeVecX>timeThresh,1)));
%     idxToThreshXSingle = find((motionVectorX - motionVectorX(find(singleTimeVecX>0,1))) > moveThresh & singleTimeVecX>0 & singleTimeVecX<=3,1);
    idxToThreshXSingle = find(motionVectorX > moveThresh & singleTimeVecX<=3,1);
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

h(9) = figure('Color','White','Name','Figure 3c','NumberTitle','off');
subplot(1,2,1)
% maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
f = fill([3 0 0 3],[3.9 3.9 .1 .1],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecX,meanX,'k')
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
% plot([0 0],[0 4],'r')
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[1,0.5,0,0.1])
end
% plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
% text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean EMG During Locomotion Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('EMG Power (au)')
ylim([0 4])
xlim([-2 3])
grid on

subplot(1,2,2)
% maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
f = fill([.25 0 0 .25],[3.9 3.9 .1 .1],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecX,meanX,'k')
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
% plot([0 0],[0 4],'r')
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[1,0.5,0,0.1])
end
% plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
% text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean EMG During Locomotion Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('EMG Power (au)')
ylim([0 4])
xlim([-1 .25])
grid on

% subplot(2,2,3)
% plot(timeVecY,meanY,'k')
% hold on
% f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,[.5 .5 .5],'Linestyle','none');
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
for n = 1:size(locDataCell,1)
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
% f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
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
% f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,[.5 .5 .5],'Linestyle','none');
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

% % % h(10) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% % % % subplot(2,1,1)
% % % histogram(timeToThreshX,binEdgesTime);
% % % hold on
% % % [pdfXVals,pdfYVals] = findKernelPDF(timeToThreshX,binEdgesTime);
% % % plot(pdfXVals,pdfYVals*10,'r','LineWidth',2)
% % % title('Time for EMG power increase one-half order of magnitude in relation to locomotion trigger')
% % % xlabel('Time (s)')
% % % xlim([-2 3])
% % % ylim([0 60])
% % % mu = mean(timeToThreshX);
% % % sig = std(timeToThreshX);
% % % plot([mu mu],[0 60],'k','LineWidth',2);
% % % plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% % % plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% % % hold off
% % % text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

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

% % % h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% % % % subplot(2,1,1)
% % % histogram(dispTimeThreshX,binEdgesDisp);
% % % hold on
% % % [pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshX,binEdgesDisp);
% % % plot(pdfXVals,pdfYVals*12,'r','LineWidth',2)
% % % title('EMG power 0.25s before locomotion trigger')
% % % xlabel('Power (au)')
% % % xlim([0 4])
% % % ylim([0 60])
% % % hold on
% % % mu = mean(dispTimeThreshX);
% % % sig = std(dispTimeThreshX);
% % % plot([mu mu],[0 60],'k','LineWidth',2);
% % % plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% % % plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% % % hold off
% % % text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

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
% h(13) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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
% h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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
moveThresh = 1.5;
timeThresh = -.25;
binEdgesTime = -2:.1:3;
binEdgesDisp = 0:.1:4;
load('LTADataCellEMG_FS.mat')
for n = 1:size(emgEventsLocationsX,1)
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
%     dispTimeThreshX(end+1) = motionVectorX((find(singleTimeVecX>timeThresh,1))) - motionVectorX((find(singleTimeVecX>0,1)));
    dispTimeThreshX(end+1) = motionVectorX((find(singleTimeVecX>timeThresh,1)));
%     idxToThreshXSingle = find((motionVectorX - motionVectorX(find(singleTimeVecX>0,1))) > moveThresh & singleTimeVecX>0 & singleTimeVecX<=3,1);
    idxToThreshXSingle = find(motionVectorX > moveThresh & singleTimeVecX<=3,1);
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

h(9) = figure('Color','White','Name','Figure 3b','NumberTitle','off');
% subplot(2,1,1)
% maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
f = fill([3 0 0 3],[3.9 3.9 .1 .1],[.9 .9 .9],'Linestyle','none','FaceAlpha',0.5);
hold on
plot(timeVecX,meanX,'k')
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
% plot([0 0],[0 4],'r')
for n = 1:size(motionEventsLocationsX,1)
    plot(timeVecX,motionEventsLocationsX(n,:),'Color',[1,0.5,0,0.1])
end
% plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
hold off
% text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean EMG During Locomotion Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
xlabel('Time (s)')
ylabel('EMG Power (au)')
ylim([0 4])
xlim([-2 3])
grid on

% subplot(2,2,3)
% plot(timeVecY,meanY,'k')
% hold on
% f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,[.5 .5 .5],'Linestyle','none');
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
for n = 1:size(locDataCell,1)
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
% f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
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
% f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,[.5 .5 .5],'Linestyle','none');
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

% % % h(10) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% % % % subplot(2,1,1)
% % % histogram(timeToThreshX,binEdgesTime);
% % % hold on
% % % [pdfXVals,pdfYVals] = findKernelPDF(timeToThreshX,binEdgesTime);
% % % plot(pdfXVals,pdfYVals*2,'r','LineWidth',2)
% % % title('Time for EMG power increase one-half order of magnitude in relation to locomotion trigger')
% % % xlabel('Time (s)')
% % % xlim([-2 3])
% % % ylim([0 20])
% % % mu = mean(timeToThreshX);
% % % sig = std(timeToThreshX);
% % % plot([mu mu],[0 60],'k','LineWidth',2);
% % % plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% % % plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% % % hold off
% % % text(-2,15,['displacement thresh = ' num2str(moveThresh) ', n = ' num2str(length(timeToThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

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

% % % h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% % % % subplot(2,1,1)
% % % histogram(dispTimeThreshX,binEdgesDisp);
% % % hold on
% % % [pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshX,binEdgesDisp);
% % % plot(pdfXVals,pdfYVals*2,'r','LineWidth',2)
% % % title('EMG power 0.25s before locomotion trigger')
% % % xlabel('Power (au)')
% % % xlim([0 4])
% % % ylim([0 20])
% % % hold on
% % % mu = mean(dispTimeThreshX);
% % % sig = std(dispTimeThreshX);
% % % plot([mu mu],[0 60],'k','LineWidth',2);
% % % plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
% % % plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
% % % hold off
% % % text(-2,30,['time thresh = ' num2str(timeThresh) ', n = ' num2str(length(dispTimeThreshX)) 10 'mean = ' num2str(mu) ', std = ' num2str(sig)])

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
% h(13) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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
% h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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

function plotLocomotionTriggeredAvgSingleTrial_FS(movementData)
load('squeezeDataCell_FS.mat')
%     locIdx = find(strcmp(squeezeDataCell(:,1),'D:/21-12-16_MouseExp/211216_002_processe_2layerBrainInSkullDataFinal.mat'));
locIdx = 1;
%     targetPositionInSkull = combineBrainSkullMovement_FS(movementData,stationaryData);
% movementData.secondsPerFrame = movementData.secondsPerFrame/2;
%     posL2 = [movementData.emgData(:,1), movementData.emgData(:,2)];
%     meanPosValL2 = [.5*(posL2(1:end-1,1) + posL2(2:end,1)),.5*(posL2(1:end-1,2) + posL2(2:end,2))];
%     posL2 = zipperVecs(posL2,meanPosValL2);
    emgEventsLocationsX = [];
%     motionEventsLocationsY = [];
%     deleteRows = [];
% movementData.secondsPerFrame = 0.0252*2;
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
        locTriggerEMGIdx = find(min(abs(locTriggerTime - timeVec)) == abs(locTriggerTime - timeVec));
%         locDataCell{i,2}(n,4:6) = [locTriggerEMGIdx-60 locTriggerEMGIdx locTriggerEMGIdx+90];
        emgVector = movementData.targetPosition(locTriggerEMGIdx-round(2/movementData.secondsPerFrame):locTriggerEMGIdx+round(3/movementData.secondsPerFrame),2);
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
moveThresh = 1.5;
timeThresh = -.25;
binEdgesTime = -2:.1:3;
binEdgesDisp = 0:.1:4;
load('LTADataCellEMG_FS.mat')
for n = 1:size(emgEventsLocationsX,1)
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
%     dispTimeThreshX(end+1) = motionVectorX((find(singleTimeVecX>timeThresh,1))) - motionVectorX((find(singleTimeVecX>0,1)));
    dispTimeThreshX(end+1) = motionVectorX((find(singleTimeVecX>timeThresh,1)));
%     idxToThreshXSingle = find((motionVectorX - motionVectorX(find(singleTimeVecX>0,1))) > moveThresh & singleTimeVecX>0 & singleTimeVecX<=3,1);
    idxToThreshXSingle = find(motionVectorX > moveThresh & singleTimeVecX<=3,1);
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

h(9) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
% maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
plot(timeVecX,meanX,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
%plot([0 0],[0 4],'r')
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
% f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,[.5 .5 .5],'Linestyle','none');
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
for n = 1:size(locDataCell,1)
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
% f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,[.5 .5 .5],'Linestyle','none');
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
% f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,[.5 .5 .5],'Linestyle','none');
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

h(10) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
histogram(timeToThreshX,binEdgesTime);
hold on
[pdfXVals,pdfYVals] = findKernelPDF(timeToThreshX,binEdgesTime);
plot(pdfXVals,pdfYVals*2,'r','LineWidth',2)
title('Time for EMG power increase one-half order of magnitude in relation to locomotion trigger')
xlabel('Time (s)')
xlim([-2 3])
ylim([0 20])
mu = mean(timeToThreshX);
sig = std(timeToThreshX);
plot([mu mu],[0 60],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
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

h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(2,1,1)
histogram(dispTimeThreshX,binEdgesDisp);
hold on
[pdfXVals,pdfYVals] = findKernelPDF(dispTimeThreshX,binEdgesDisp);
plot(pdfXVals,pdfYVals*2,'r','LineWidth',2)
title('EMG power 0.25s before locomotion trigger')
xlabel('Power (au)')
xlim([0 4])
ylim([0 20])
hold on
mu = mean(dispTimeThreshX);
sig = std(dispTimeThreshX);
plot([mu mu],[0 60],'k','LineWidth',2);
plot([mu+sig mu+sig],[0 60],'k--','LineWidth',2);
plot([mu-sig mu-sig],[0 60],'k--','LineWidth',2);
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
% h(13) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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
% h(50) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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

h(1) = figure('Color','White','Name','Supplementary Figure 2b','NumberTitle','off');
subplot(1,2,1)
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

subplot(1,2,2)
% h(2) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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
h(3) = figure('Color','White','Name','Supplementary Figure 2c','NumberTitle','off');
subplot(1,2,1)
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
% h(4) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
subplot(1,2,2)
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

function plotPSF_FS
load('PSFData_FS.mat')
h(1) = figure('Color','White','Name','Supplementary Figure 3c','NumberTitle','off');
subplot(1,2,1)
plot(PSFData.X,PSFData.XIntensity,'k')
hold on
scatter(PSFData.XPoints,PSFData.XPointsIntensity,'kx')
hold off
xlabel('\mum')
ylabel('Normalized Pixel Intensity')
title('X')
xlim([-2 2])
ylim([0 1])

subplot(1,2,2)
plot(PSFData.ZIntensity,PSFData.Z,'k')
hold on
scatter(PSFData.ZPointsIntensity,PSFData.ZPoints,'kx')
hold off
ylabel('\mum')
xlabel('Normalized Pixel Intensity')
title('Z')
xlim([0 1])
ylim([-7 7])
set(gca,'YDir','reverse')
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

h(1) = figure('Color','White','Name','Supplementary Figure 3f','NumberTitle','off');
subplot(2,1,1)
hold on
for n = 1:numel(fn)
    plot(calibrationValues.(fn{n}).diopterVals,calibrationValues.(fn{n}).zMatchMicrons - yZeroDiff,'-*');
end
f = fill([stdX, fliplr(stdX)], [stdYPlus, fliplr(stdYMinus)], [.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
hold off
xlabel('Diopter Input (m^{-1})')
ylabel('Z (\mum)')
% title('Controller Diopter Values vs. Focal Plane Position in Z')
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
xlabel('Diopter Input (m^{-1})')
ylabel('Z (\mum)')
% title('Controller Diopter Values vs. Focal Plane Position in Z')
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

h(2) = figure('Color','White','Name','Supplementary Figure 3g','NumberTitle','off');
subplot(3,1,1)
plot(diopterVals,widthVals,'mo')
xlabel('Diopter Input (m^{-1})')
ylabel('\mum/Pixel in X')
xlim([-1.5 2])
ylim([.7 .9])
lsline
subplot(3,1,2)
plot(diopterVals,heightVals,'go')
xlabel('Diopter Input (m^{-1})')
ylabel('\mum/Pixel in Y')
xlim([-1.5 2])
ylim([.5 .7])
lsline
subplot(3,1,3)
plot(diopterValsPower,powerVals,'ko')
xlabel('Diopter Input (m^{-1})')
ylabel('Laser Power Through Objective (mW)')
xlim([-1.5 2])
ylim([46 52])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMotionTrackingBrainAndSkullResp_FS(movementData,stationaryData,figureTitle)
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
h(6) = figure('Color','White','Name',figureTitle,'NumberTitle','off');
subplot(4,1,1)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
hold off
% title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
xlabel('Time (s)')
ylabel('Medial-Lateral Shift (\mum)')
grid on
axis([15 225 -3 3])
text(15,3,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(15,-3,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(4,1,2)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
hold off
xlabel('Time (s)')
ylabel('Rostral-Caudal Shift (\mum)')
grid on
axis([15 225 -3 3])
text(15,3,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(15,-3,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(4,1,3)
plot(movementData.emgData(:,1),movementData.emgData(:,2),'Color',[1 0.5 0])
xlabel('Time (s)')
ylabel('Abdominal EMG (au)')
grid on
axis([15 225 0.5 2])
subplot(4,1,4)
plot(movementData.videoRespiration(:,1)-1.7,movementData.videoRespiration(:,2),'Color',[0.8549 0.1098 0.3607])
xlabel('Time (s)')
ylabel('Respiration (mean pixel intensity)')
grid on
axis([15 225 50 100])

if contains(figureTitle,'Supplementary Figure')
    h(6) = figure('Color','White','Name','Supplementary Figure 4d_2','NumberTitle','off');
    subplot(4,1,1)
    plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
    hold on
    plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
    hold off
    % title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
    xlabel('Time (s)')
    ylabel('Medial-Lateral Shift (\mum)')
    grid on
    axis([15 25 -3 3])
    text(15,3,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(15,-3,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    subplot(4,1,2)
    plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
    hold on
    plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
    hold off
    xlabel('Time (s)')
    ylabel('Rostal-Caudal Shift (\mum)')
    grid on
    axis([15 25 -3 3])
    text(15,3,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(15,-3,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    subplot(4,1,3)
    plot(movementData.emgData(:,1),movementData.emgData(:,2),'Color',[1 .5 0])
    xlabel('Time (s)')
    ylabel('Abdominal EMG (au)')
    grid on
    axis([15 25 0.5 2])
    subplot(4,1,4)
    plot(movementData.videoRespiration(:,1)-1.7,movementData.videoRespiration(:,2),'Color',[0.8549 0.1098 0.3607])
    xlabel('Time (s)')
    ylabel('Respiration (mean pixel intensity)')
    grid on
    axis([15 25 50 100])
    
    h(6) = figure('Color','White','Name','Supplementary Figure 4d_3','NumberTitle','off');
    subplot(4,1,1)
    plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
    hold on
    plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
    hold off
    % title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
    xlabel('Time (s)')
    ylabel('Medial-Lateral Shift (\mum)')
    grid on
    axis([135 145 -3 3])
    text(135,3,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(135,-3,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    subplot(4,1,2)
    plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
    hold on
    plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
    hold off
    xlabel('Time (s)')
    ylabel('Rostral-Caudal Shift (\mum)')
    grid on
    axis([135 145 -3 3])
    text(135,3,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(135,-3,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    subplot(4,1,3)
    plot(movementData.emgData(:,1),movementData.emgData(:,2),'Color',[1 .5 0])
    xlabel('Time (s)')
    ylabel('Abdominal EMG (au)')
    grid on
    axis([135 145 0.5 2])
    subplot(4,1,4)
    plot(movementData.videoRespiration(:,1)-1.7,movementData.videoRespiration(:,2),'Color',[0.8549 0.1098 0.3607])
    xlabel('Time (s)')
    ylabel('Respiration (mean pixel intensity)')
    grid on
    axis([135 145 50 100])
end
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
h(6) = figure('Color','White','Name','Supplementary Figure 4e_1','NumberTitle','off');
subplot(4,1,1)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
hold off
% title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
xlabel('Time (s)')
ylabel('Lateral-Medial Shift (\mum)')
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
ylabel('Rostral-Caudal Shift (\mum)')
grid on
axis([62 272 -6 6])
text(62,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(62,-6,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(4,1,3)
plot(movementData.emgData(:,1),movementData.emgData(:,2),'Color',[1 0.5 0])
xlabel('Time (s)')
ylabel('Abdominal EMG (au)')
grid on
axis([62 272 0.5 3])
subplot(4,1,4)
plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)*2*pi*.06),'k')
xlabel('Time (s)')
ylabel('Treadmill Velocity (m/s)')
grid on
axis([62 272 0 0.2])

h(6) = figure('Color','White','Name','Supplementary Figure 4e_2','NumberTitle','off');
subplot(4,1,1)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
hold off
% title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
xlabel('Time (s)')
ylabel('Medial-Lateral Shift (\mum)')
grid on
axis([66 76 -6 6])
text(66,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(66,-6,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(4,1,2)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
hold off
xlabel('Time (s)')
ylabel('Rostral-Caudal Shift (\mum)')
grid on
axis([66 76 -6 6])
text(66,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(66,-6,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(4,1,3)
plot(movementData.emgData(:,1),movementData.emgData(:,2),'Color',[1 .5 0])
xlabel('Time (s)')
ylabel('Abdominal EMG (au)')
grid on
axis([66 76 0.5 3])
subplot(4,1,4)
plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)*2*pi*.06),'k')
xlabel('Time (s)')
ylabel('Treadmill Velocity (m/s)')
grid on
axis([66 76 0 0.2])

h(6) = figure('Color','White','Name','Supplementary Figure 4e_3','NumberTitle','off');
subplot(4,1,1)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
hold off
% title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
xlabel('Time (s)')
ylabel('Medial-Lateral Shift (\mum)')
grid on
axis([150 160 -6 6])
text(150,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(150,-6,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(4,1,2)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
hold off
xlabel('Time (s)')
ylabel('Rostral-Caudal Shift (\mum)')
grid on
axis([150 160 -6 6])
text(150,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(150,-6,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(4,1,3)
plot(movementData.emgData(:,1),movementData.emgData(:,2),'Color',[1 .5 0])
xlabel('Time (s)')
ylabel('Abdominal EMG (au)')
grid on
axis([150 160 0.5 3])
subplot(4,1,4)
plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)*2*pi*.06),'k')
xlabel('Time (s)')
ylabel('Treadmill Velocity (m/s)')
grid on
axis([150 160 0 0.2])

% h(6) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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
% h(6) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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

function plotMotionTracking2P2LPCASqueeze_FS(movementData,stationaryData,figureTitle)
targetPositionInSkull = combineBrainSkullMovement_FS(movementData,stationaryData);
motionVec = pcaMotionAnalysis_FS(targetPositionInSkull);
h(2) = figure('Color','White','Name',figureTitle,'NumberTitle','off');
s = scatter(targetPositionInSkull(:,1),targetPositionInSkull(:,2),10,'g','filled');
s.MarkerFaceAlpha = .5;
hold on
drawArrow_FS([0;0],[motionVec(1);motionVec(2)]);
hold off
axis equal square
axis([-4 4 -4 4])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
% title(['\fontsize{20pt}\bf{Position of Brain in Skull with PCA Vector}'])
% title(['\fontsize{20pt}\bf{Figure 2a}'])
xlabel('\mum')
ylabel('\mum')
if movementData.hemisphere == 1
    text(4,0,'Lateral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(-4,0,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    text(0,4,'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(0,-4,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
else
    text(4,0,'Medial','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(-4,0,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    text(0,4,'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(0,-4,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMotionTrackingBrainAndSkullSqueeze_FS(movementData,stationaryData,figureTitle)
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
h(6) = figure('Color','White','Name',figureTitle,'NumberTitle','off');
subplot(3,1,1)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
hold off
% title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
xlabel('Time (s)')
ylabel('Lateral-Medial Shift (\mum)')
grid on
axis([0 150 -1 4])
text(0,4,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(0,-1,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(3,1,2)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
hold on
plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
hold off
xlabel('Time (s)')
ylabel('Rostral-Caudal Shift (\mum)')
grid on
axis([0 150 -1 4])
text(0,4,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(0,-1,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(3,1,3)
[~,locs] = findpeaks(movementData.ballData(:,2),movementData.ballData(:,1),'MinPeakDistance',1.5,'MinPeakHeight',max(movementData.ballData(:,2))-.2);
binarySqueeze = zeros(length(movementData.ballData(:,1)),1);
for n = 1:2:length(locs)
    startIdx = find(movementData.ballData(:,1)==locs(n));
    stopIdx = find(movementData.ballData(:,1)==locs(n+1));
    binarySqueeze(startIdx:stopIdx) = 1;
end
plot(movementData.ballData(:,1),binarySqueeze,'Color',[.5 0 .5])
xlabel('Time (s)')
ylabel('Two Second Abdominal Compression')
grid on
axis([0 150 -0.1 1])

% h(6) = figure('Color','White','Name','Figure 4e_2','NumberTitle','off');
% subplot(4,1,1)
% plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
% hold on
% plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
% hold off
% % title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
% xlabel('Time (s)')
% ylabel('Medial-Lateral Shift (\mum)')
% grid on
% axis([66 76 -6 6])
% text(62,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(62,-6,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(4,1,2)
% plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
% hold on
% plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
% hold off
% xlabel('Time (s)')
% ylabel('Rostral-Caudal Shift (\mum)')
% grid on
% axis([66 76 -6 6])
% text(62,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(62,-6,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(4,1,3)
% plot(movementData.emgData(:,1),movementData.emgData(:,2),'Color',[1 .5 0])
% xlabel('Time (s)')
% ylabel('Abdominal EMG (au)')
% grid on
% axis([66 76 0.5 3])
% subplot(4,1,4)
% plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)*2*pi*.06),'k')
% xlabel('Time (s)')
% ylabel('Treadmill Velocity (m/s)')
% grid on
% axis([66 76 0 0.2])
% 
% h(6) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(4,1,1)
% plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
% hold on
% plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
% hold off
% % title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
% xlabel('Time (s)')
% ylabel('Medial-Lateral Shift (\mum)')
% grid on
% axis([150 160 -6 6])
% text(62,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(62,-6,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(4,1,2)
% plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
% hold on
% plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
% hold off
% xlabel('Time (s)')
% ylabel('Rostral-Caudal Shift (\mum)')
% grid on
% axis([150 160 -6 6])
% text(62,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(62,-6,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(4,1,3)
% plot(movementData.emgData(:,1),movementData.emgData(:,2),'Color',[1 .5 0])
% xlabel('Time (s)')
% ylabel('Abdominal EMG (au)')
% grid on
% axis([150 160 0.5 3])
% subplot(4,1,4)
% plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)*2*pi*.06),'k')
% xlabel('Time (s)')
% ylabel('Treadmill Velocity (m/s)')
% grid on
% axis([150 160 0 0.2])

% h(6) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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
% h(6) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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

function plotMotionTrackingSkullBregma_FS(movementData)
if movementData.hemisphere == 2
    movementData.targetPosition(:,1) = movementData.targetPosition(:,1)*-1;
end
movementData.targetPosition(:,2) = movementData.targetPosition(:,2)*-1;
h(6) = figure('Color','White','Name','Supplementary Figure 6a','NumberTitle','off');
subplot(3,1,1)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'r')
% title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
xlabel('Time (s)')
ylabel('Bregma Window Marker on Skull (\mum)')
grid on
axis([40 140 -1 5])
text(40,5,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(40,-1,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(3,1,2)
plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'b')
xlabel('Time (s)')
ylabel('Bregma Window Marker on Skull (\mum)')
grid on
axis([40 140 -1 5])
text(40,5,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(40,-1,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
subplot(3,1,3)
plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)*2*pi*.06),'k')
xlabel('Time (s)')
ylabel('Treadmill Velocity (m/s)')
grid on
axis([40 140 0 0.2])

% h(6) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(4,1,1)
% plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
% hold on
% plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
% hold off
% title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
% xlabel('Time (s)')
% ylabel('X Position (\mum)')
% grid on
% axis([66 76 -6 6])
% text(62,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(62,-6,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(4,1,2)
% plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
% hold on
% plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
% hold off
% xlabel('Time (s)')
% ylabel('Y Position (\mum)')
% grid on
% axis([66 76 -6 6])
% text(62,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(62,-6,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(4,1,3)
% plot(movementData.emgData(:,1),movementData.emgData(:,2),'k')
% xlabel('Time (s)')
% ylabel('Abdominal EMG (au)')
% grid on
% axis([66 76 0.5 3])
% subplot(4,1,4)
% plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)*2*pi*.06),'r')
% xlabel('Time (s)')
% ylabel('Locomotion (m/s)')
% grid on
% axis([66 76 0 0.2])
% 
% h(6) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% subplot(4,1,1)
% plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,1),'g')
% hold on
% plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,1),'m')
% hold off
% title(['Figure 1e' 10 '\fontsize{20pt}\bf{Position of Brain and Skull}'])
% xlabel('Time (s)')
% ylabel('X Position (\mum)')
% grid on
% axis([150 160 -6 6])
% text(62,6,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(62,-6,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(4,1,2)
% plot([1:size(movementData.targetPosition,1)]*movementData.secondsPerFrame,movementData.targetPosition(:,2),'g')
% hold on
% plot([1:size(stationaryData.targetPosition,1)]*movementData.secondsPerFrame,stationaryData.targetPosition(:,2),'m')
% hold off
% xlabel('Time (s)')
% ylabel('Y Position (\mum)')
% grid on
% axis([150 160 -6 6])
% text(62,6,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(62,-6,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% subplot(4,1,3)
% plot(movementData.emgData(:,1),movementData.emgData(:,2),'k')
% xlabel('Time (s)')
% ylabel('Abdominal EMG (au)')
% grid on
% axis([150 160 0.5 3])
% subplot(4,1,4)
% plot(movementData.ballData(:,1),abs(movementData.ballData(:,2)*2*pi*.06),'r')
% xlabel('Time (s)')
% ylabel('Locomotion (m/s)')
% grid on
% axis([150 160 0 0.2])

% h(6) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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
% h(6) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
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

function plotAbdominalCompressionPressure_FS
circumMeters = 0.115;
widthMeters = 0.015;
surfArea = circumMeters*widthMeters;

forceCalArb = [1800;4400;7300;10000;12000;13700];
forceCalGrams = [50;100;150;200;250;300];
% forceCalNewtons = (forceCalGrams./1000)*9.81;
% forceCalPascals = forceCalNewtons./surfArea;
% forceCalMMHg = forceCalPascals/133.322;
% m = forceCalArb\forceCalMMHg;
% m = forceCalArb\forceCalGrams;
% plot(forceCalArb,forceCalMMHg)
% hold on
% plot([0;forceCalArb],[0;forceCalArb]*m)

h(6) = figure('Color','White','Name','Supplementary Figure 15b','NumberTitle','off');
subplot(1,4,1:3)
load('abdCompressPressureTrial.mat');
abdCompressPressureTrial(:,1) = (abdCompressPressureTrial(:,1)-abdCompressPressureTrial(1,1))/1000;
xInd = find(abdCompressPressureTrial(:,1)>150 & abdCompressPressureTrial(:,1)<300);
xStart = xInd(1);
plotLineX = abdCompressPressureTrial(xInd,1)-abdCompressPressureTrial(xStart,1);
% plotLineY = (abdCompressPressureTrial(xInd,2)-abdCompressPressureTrial(xStart,2)).*m;
pressureMax = max(abdCompressPressureTrial(xInd,2)-abdCompressPressureTrial(xStart,2));
plotLineY = (abdCompressPressureTrial(xInd,2)-abdCompressPressureTrial(xStart,2))./pressureMax;
plot(plotLineX,plotLineY)
xlim([0 150])
ylim([-0.1 1])
% title('Abdominal Compression Pressure')
xlabel('Time (s)')
ylabel('A.U.')

% h(6) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
subplot(1,4,4)
load('abdCompressPressureTrial.mat');
abdCompressPressureTrial(:,1) = (abdCompressPressureTrial(:,1)-abdCompressPressureTrial(1,1))/1000;
xInd = find(abdCompressPressureTrial(:,1)>150 & abdCompressPressureTrial(:,1)<300);
xStart = xInd(1);
plotLineX = abdCompressPressureTrial(xInd,1)-abdCompressPressureTrial(xStart,1);
% plotLineY = (abdCompressPressureTrial(xInd,2)-abdCompressPressureTrial(xStart,2)).*m;
pressureMax = max(abdCompressPressureTrial(xInd,2)-abdCompressPressureTrial(xStart,2));
plotLineY = (abdCompressPressureTrial(xInd,2)-abdCompressPressureTrial(xStart,2))./pressureMax;
plot(plotLineX,plotLineY)
xlim([101 106])
ylim([-0.1 1])
% title('Abdominal Compression Pressure')
xlabel('Time (s)')
ylabel('A.U.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMouseStats_FS
% fileName = 'C:\Workspace\Code\DrewLab\movementDataLog.mat';
load('movementDataLog_FS.mat')

moveMag = [];
weigth = [];
age = [];
sex = [];
for n = 1:size(moveDataMat,1)
    moveMag(n) = sqrt(moveDataMat{n,7}(1)^2 + moveDataMat{n,7}(2)^2);
    weight(n) = str2double(moveDataMat{n,18});
    age(n) = (datenum(moveDataMat{n,1}(1:6),'yymmdd')-datenum(moveDataMat{n,17},'mm/dd/yy'))/7;
    if strcmp(moveDataMat{n,16},'M')
        sex(n) = 1;
    else
        sex(n) = 0;
    end
end

weightX = min(weight):1:max(weight);
ageX = min(age):1:max(age);

h(6) = figure('Color','White','Name','Supplementary Figure 8a','NumberTitle','off');
% subplot(2,3,1)
s = scatter(weight(logical(sex)),moveMag(logical(sex)),50,'b','filled',"DisplayName","Male");
s.MarkerFaceAlpha = .15;
hold on
mdl = fitlm(weight(logical(sex)),moveMag(logical(sex)));
% plot(weight(logical(sex)),mdl.Fitted,"b-","LineWidth",2,"DisplayName","Male Fit") ;
plot(weightX,mdl.Coefficients{1,1}+(weightX*mdl.Coefficients{2,1}),"b-","LineWidth",2,"DisplayName","Male Fit") ;
% [~,yCI]=predict(mdl,[min(weight(logical(sex))):1:max(weight(logical(sex)))]',"Prediction","curve");
[~,yCI]=predict(mdl,weightX',"Prediction","curve");
% [~,yci_obs]=predict(mdl,x,"Prediction","observation");
% plot([min(weight(logical(sex))):1:max(weight(logical(sex)))]',yCI,"b--","DisplayName","Male CI") ;
plot(weightX',yCI,"b--","DisplayName","Male CI") ;
text(26,13,['p: ' num2str(mdl.Coefficients{2,4}) ', R^{2}: ' num2str(mdl.Rsquared.Ordinary)],'Color','b')
% plot(x,yci_obs,"ro","DisplayName","Prediction bounds") ;
s = scatter(weight(~logical(sex)),moveMag(~logical(sex)),50,'r','filled',"DisplayName","Female");
s.MarkerFaceAlpha = .15;
mdl = fitlm(weight(~logical(sex)),moveMag(~logical(sex)));
% plot(weight(~logical(sex)),mdl.Fitted,"r-","LineWidth",2,"DisplayName","Female Fit") ;
plot(weightX,mdl.Coefficients{1,1}+(weightX.*mdl.Coefficients{2,1}),"r-","LineWidth",2,"DisplayName","Female Fit") ;
% [~,yCI]=predict(mdl,[min(weight(logical(sex))):1:max(weight(logical(sex)))]',"Prediction","curve");
[~,yCI]=predict(mdl,weightX',"Prediction","curve");
% [~,yci_obs]=predict(mdl,x,"Prediction","observation");
% plot([min(weight(logical(sex))):1:max(weight(logical(sex)))]',yCI,"r--","DisplayName","Female CI") ;
plot(weightX',yCI,"r--","DisplayName","Female CI") ;
text(26,12,['p: ' num2str(mdl.Coefficients{2,4}) ', R^{2}: ' num2str(mdl.Rsquared.Ordinary)],'Color','r')
% plot(x,yci_obs,"ro","DisplayName","Prediction bounds") ;
hold off
% title('Mouse Weight vs Brain Shift')
xlabel('Mouse Weight (g)')
ylabel('Brain Shift (\mum)')
xlim([25 65])
ylim([0 15])
legend

h(6) = figure('Color','White','Name','Supplementary Figure 8b','NumberTitle','off');
% subplot(2,3,2)
s = scatter(age(logical(sex)),moveMag(logical(sex)),50,'b','filled',"DisplayName","Male");
s.MarkerFaceAlpha = .15;
hold on
mdl = fitlm(age(logical(sex)),moveMag(logical(sex)));
% plot(age(logical(sex)),mdl.Fitted,"b-","LineWidth",2,"DisplayName","Male Fit") ;
plot(ageX,mdl.Coefficients{1,1}+(ageX*mdl.Coefficients{2,1}),"b-","LineWidth",2,"DisplayName","Male Fit") ;
% [~,yCI]=predict(mdl,[min(age(logical(sex))):1:max(age(logical(sex)))]',"Prediction","curve");
[~,yCI]=predict(mdl,ageX',"Prediction","curve");
% [~,yci_obs]=predict(mdl,x,"Prediction","observation");
% plot([min(age(logical(sex))):1:max(age(logical(sex)))]',yCI,"b--","DisplayName","Male CI") ;
plot(ageX',yCI,"b--","DisplayName","Male CI") ;
text(1,13,['p: ' num2str(mdl.Coefficients{2,4}) ', R^{2}: ' num2str(mdl.Rsquared.Ordinary)],'Color','b')
% plot(x,yci_obs,"ro","DisplayName","Prediction bounds") ;
s = scatter(age(~logical(sex)),moveMag(~logical(sex)),50,'r','filled',"DisplayName","Female");
s.MarkerFaceAlpha = .15;
mdl = fitlm(age(~logical(sex)),moveMag(~logical(sex)));
% plot(age(~logical(sex)),mdl.Fitted,"r-","LineWidth",2,"DisplayName","Female Fit") ;
plot(ageX,mdl.Coefficients{1,1}+(ageX.*mdl.Coefficients{2,1}),"r-","LineWidth",2,"DisplayName","Female Fit") ;
% [~,yCI]=predict(mdl,[min(age(logical(sex))):1:max(age(logical(sex)))]',"Prediction","curve");
[~,yCI]=predict(mdl,ageX',"Prediction","curve");
% [~,yci_obs]=predict(mdl,x,"Prediction","observation");
% plot([min(age(logical(sex))):1:max(age(logical(sex)))]',yCI,"r--","DisplayName","Female CI") ;
plot(ageX',yCI,"r--","DisplayName","Female CI") ;
text(1,12,['p: ' num2str(mdl.Coefficients{2,4}) ', R^{2}: ' num2str(mdl.Rsquared.Ordinary)],'Color','r')
% plot(x,yci_obs,"ro","DisplayName","Prediction bounds") ;
hold off
% title('Mouse Age vs Brain Shift')
xlabel('Mouse Age (Weeks)')
ylabel('Brain Shift (\mum)')
xlim([0 70])
ylim([0 15])
legend

h(6) = figure('Color','White','Name','Supplementary Figure 8c','NumberTitle','off');
% subplot(2,3,3)
plotSpread_FS({moveMag(logical(sex)),moveMag(~logical(sex))},'categoryIdx',[zeros(length(moveMag(logical(sex))),1);ones(length(moveMag(~logical(sex))),1)],'categoryColors',{'b','r'},'categoryLabels',{'Male','Female'},'xNames',{'Male','Female'},'showMM',5);
meanMale = mean(moveMag(logical(sex)));
stdMale = std(moveMag(logical(sex)));
meanFemale = mean(moveMag(~logical(sex)));
stdFemale = std(moveMag(~logical(sex)));
% text(1,12,['$\bar{x}_{Male}$: ' num2str(meanMale)],'Color','b','Interpreter','latex')
% text(1,11,['\pm\sigma:' num2str(stdMale)],'Color','b')
% text(2,12,['$\bar{x}_{Male}$: ' num2str(meanFemale)],'Color','r','Interpreter','latex')
% text(2,11,['\pm\sigma:' num2str(stdFemale)],'Color','r')
text(.3,12,['Mean: ' num2str(meanMale) ', \pm\sigma: ' num2str(stdMale)],'Color','b')
text(1.6,12,['Mean: ' num2str(meanFemale) ', \pm\sigma: ' num2str(stdFemale)],'Color','r')
[h,p] = kstest2(moveMag(logical(sex)),moveMag(~logical(sex)));
text(1.2,13,['KS Test: ' num2str(h) ', p: ' num2str(p)])
% title('Mouse Sex vs Brain Shift')
ylabel('Brain Shift (\mum)')
xlabel('Mouse Sex')

% h(6) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% % subplot(2,3,4)
% removeDup = unique([age',weight'],'rows');
% scatter(removeDup(:,1),removeDup(:,2))
% title('Mouse Age vs Weight')
% xlabel('Mouse Age (Weeks)')
% ylabel('Mouse Weight (g)')
% xlim([0 70])
% 
% h(6) = figure('Color','White','Name','Figure 1e','NumberTitle','off');
% % subplot(2,3,5)
% removeDup = unique([sex',weight'],'rows');
% sexRem = removeDup(:,1);
% weightRem = removeDup(:,2);
% plotSpread_FS({weightRem(~logical(sexRem)),weightRem(logical(sexRem))},'categoryIdx',[zeros(length(weightRem(~logical(sexRem))),1);ones(length(weightRem(logical(sexRem))),1)],'categoryColors',{'r','b'},'categoryLabels',{'F','M'});
% title('Mouse Sex vs Weight')
% ylabel('Mouse Weight (g)')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotXCorrData_FS
load('xCorrData_FS.mat')
xCorrXLoc = [];
xCorrYLoc = [];
xCorrXEMG = [];
xCorrYEMG = [];
xLoc = [];
yLoc = [];
mouseNum = [];
% xCorrData = {};
maxlag = 500;
for i = 1:size(xCorrData,1)
%     load([fileName{i}(1:31) '.txt'])
%     load(fileName{i});
    xLoc(i) = xCorrData{i,6};
    yLoc(i) = xCorrData{i,7};
    mouseNum(i) = xCorrData{i,8};
    motionVectorX = xCorrData{i,2};
    motionVectorY = xCorrData{i,3};
%     if movementData.hemisphere == 2
%         movementData.targetPosition(:,1) = movementData.targetPosition(:,1) * -1;
%     end
%     movementData.secondsPerFrame = movementData.secondsPerFrame/2;
locDataInterp = xCorrData{i,4};
    emgDataInterp = xCorrData{i,5};
%     movement_time=movementData.secondsPerFrame*(1:length(motionVectorX));
%     locDataInterp=zeros(size(movementData.targetPosition));
%     locDataInterp(:,1)=movement_time;
%     locDataInterp(:,2)=interp1(movementData.ballData(:,1),abs(movementData.ballData(:,2)),movement_time,'linear');
%     emgDataInterp=zeros(size(movementData.targetPosition));
%     emgDataInterp(:,1)=movement_time;
%     emgDataInterp(:,2)=interp1(movementData.emgData(:,1),abs(movementData.emgData(:,2)),movement_time,'linear');
%     xCorrData(i,:) = {fileName{i},motionVectorX,motionVectorY,locDataInterp(:,2),emgDataInterp(:,2),xLoc(i),yLoc(i),mouseNum(i)};
%     xc_1=xcorr(detrend(locDataInterp(1:(end-100),2)-mean(locDataInterp(1:(end-100),2)))', detrend(movementData.targetPosition(1:(end-100),1)-mean(movementData.targetPosition(1:(end-100),1)))',maxlag,'coeff');
%     xc_2=xcorr(detrend(locDataInterp(1:(end-100),2)-mean(locDataInterp(1:(end-100),2)))', detrend(movementData.targetPosition(1:(end-100),2)-mean(movementData.targetPosition(1:(end-100),2)))',maxlag,'coeff');
%     xc_3=xcorr(detrend(emgDataInterp(1:(end-100),2)-mean(emgDataInterp(1:(end-100),2)))', detrend(movementData.targetPosition(1:(end-100),1)-mean(movementData.targetPosition(1:(end-100),1)))',maxlag,'coeff');
%     xc_4=xcorr(detrend(emgDataInterp(1:(end-100),2)-mean(emgDataInterp(1:(end-100),2)))', detrend(movementData.targetPosition(1:(end-100),2)-mean(movementData.targetPosition(1:(end-100),2)))',maxlag,'coeff');
    xc_1=xcorr(detrend(locDataInterp(1:(end-100)))', detrend(motionVectorX(1:(end-100)))',maxlag,'coeff');
    xc_2=xcorr(detrend(locDataInterp(1:(end-100)))', detrend(motionVectorY(1:(end-100)))',maxlag,'coeff');
    xc_3=xcorr(detrend(emgDataInterp(1:(end-100)))', detrend(motionVectorX(1:(end-100)))',maxlag,'coeff');
    xc_4=xcorr(detrend(emgDataInterp(1:(end-100)))', detrend(motionVectorY(1:(end-100)))',maxlag,'coeff');
%     xc_1=xcorr(detrend(locDataInterp(1:(end-100),2)).^2', detrend(movementData.targetPosition(1:(end-100),1)).^2',maxlag,'coeff');
%     xc_2=xcorr(detrend(locDataInterp(1:(end-100),2)).^2', detrend(movementData.targetPosition(1:(end-100),2)).^2',maxlag,'coeff');
%     xc_3=xcorr(detrend(emgDataInterp(1:(end-100),2)).^2', detrend(movementData.targetPosition(1:(end-100),1)).^2',maxlag,'coeff');
%     xc_4=xcorr(detrend(emgDataInterp(1:(end-100),2)).^2', detrend(movementData.targetPosition(1:(end-100),2)).^2',maxlag,'coeff');
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
%     if ~any(isnan(xc_1))
%         xCorrXLoc(end+1,:) = xc_1;
%     end
%     if ~any(isnan(xc_2))
%         xCorrYLoc(end+1,:) = xc_2;
%     end
%     if ~any(isnan(xc_3))
%         xCorrXEMG(end+1,:) = xc_3;
%     end
%     if ~any(isnan(xc_4))
%         xCorrYEMG(end+1,:) = xc_4;
%     end
xCorrXLoc(i,:) = xc_1;
xCorrYLoc(i,:) = xc_2;
xCorrXEMG(i,:) = xc_3;
xCorrYEMG(i,:) = xc_4;
% if any(xc_2<-.5)
%     disp('check')
% end
    clear movementData
end

u = 1;
yLocUnique = [];
if u == 1
    % unique locations
    [~,~,uniRows] = unique([xLoc' yLoc' mouseNum'],'rows');
    % unique mice
%     [~,~,uniRows] = unique(mouseNum');
    xCorrXLocNew = [];
    xCorrYLocNew = [];
    xCorrXEMGNew = [];
    xCorrYEMGNew = [];
    for n = 1:max(uniRows)
        xCorrXLocUnique = xCorrXLoc(uniRows == n,:);
        xCorrYLocUnique = xCorrYLoc(uniRows == n,:);
        xCorrXEMGUnique = xCorrXEMG(uniRows == n,:);
        xCorrYEMGUnique = xCorrYEMG(uniRows == n,:);
        xCorrXLocNew(n,:) = mean(xCorrXLocUnique,1);
        xCorrYLocNew(n,:) =  mean(xCorrYLocUnique,1);
        xCorrXEMGNew(n,:) =  mean(xCorrXEMGUnique,1);
        xCorrYEMGNew(n,:) =  mean(xCorrYEMGUnique,1);
%         xLocUnique = xLoc(find(uniRows == n),'first');
        yLocUnique(n) = yLoc(find(uniRows == n,1));
    end
    xCorrXLoc = xCorrXLocNew;
    xCorrYLoc = xCorrYLocNew;
    xCorrXEMG = xCorrXEMGNew;
    xCorrYEMG = xCorrYEMGNew;
end

        
% load(fileName{1});
[meanXLoc,cIntFillPtsXLoc] = getCIntMeanAndFillPts_FS(xCorrXLoc,90);
[meanYLoc,cIntFillPtsYLoc] = getCIntMeanAndFillPts_FS(xCorrYLoc,90);
[meanXEMG,cIntFillPtsXEMG] = getCIntMeanAndFillPts_FS(xCorrXEMG,90);
[meanYEMG,cIntFillPtsYEMG] = getCIntMeanAndFillPts_FS(xCorrYEMG,90);
timeVecX = 0.0253*(-maxlag:maxlag);
timeVecY = 0.0253*(-maxlag:maxlag);


% targetPositionInSkull = combineBrainSkullMovement_FS(movementData,stationaryData);
% movementData.secondsPerFrame = movementData.secondsPerFrame/2;
% movementData.targetPosition = targetPositionInSkull;
% movement_time=movementData.secondsPerFrame*(1:length(movementData.targetPosition));
% movementData.locDataInterp=zeros(size(movementData.targetPosition));
% movementData.locDataInterp(:,1)=movement_time;
% movementData.locDataInterp(:,2)=interp1(movementData.ballData(:,1),abs(movementData.ballData(:,2)),movement_time,'linear');

% % % h(1) = figure('Color','White','Name','Figure 16a','NumberTitle','off');
% % % % subplot(2,1,1)
% % % timeInd = find(timeVecX > -1,1);
% % % corrVsPosLoc = [];
% % % for n = 1:size(xCorrXLoc,1)
% % %     corrVsPosLoc(n) = xCorrXLoc(n,timeInd);
% % % end
% % % scatter(yLocUnique,corrVsPosLoc)
% % % hold on
% % % mdl = fitlm(yLocUnique,corrVsPosLoc);
% % % % plot(age(logical(sex)),mdl.Fitted,"b-","LineWidth",2,"DisplayName","Male Fit") ;
% % % plot(-6000:1:4000,mdl.Coefficients{1,1}+((-6000:1:4000)*mdl.Coefficients{2,1}),"b-","LineWidth",2,"DisplayName","Male Fit") ;
% % % % [~,yCI]=predict(mdl,[min(age(logical(sex))):1:max(age(logical(sex)))]',"Prediction","curve");
% % % [~,yCI]=predict(mdl,(-6000:1:4000)',"Prediction","curve");
% % % % [~,yci_obs]=predict(mdl,x,"Prediction","observation");
% % % % plot([min(age(logical(sex))):1:max(age(logical(sex)))]',yCI,"b--","DisplayName","Male CI") ;
% % % plot((-6000:1:4000)',yCI,"b--","DisplayName","Male CI") ;
% % % hold off
% % % ylim([-1 1])
% % % title('Locomotion xcorr')
% % % xlabel('<Rostral  Caudal>')
% % % 
% % % h(1) = figure('Color','White','Name','Figure 16b','NumberTitle','off');
% % % % subplot(2,1,2)
% % % corrVsPosEMG = [];
% % % for n = 1:size(xCorrXLoc,1)
% % %     corrVsPosEMG(n) = xCorrXEMG(n,timeInd);
% % % end
% % % scatter(yLocUnique,corrVsPosEMG)
% % % hold on
% % % mdl = fitlm(yLocUnique,corrVsPosEMG);
% % % % plot(age(logical(sex)),mdl.Fitted,"b-","LineWidth",2,"DisplayName","Male Fit") ;
% % % plot(-6000:1:4000,mdl.Coefficients{1,1}+((-6000:1:4000)*mdl.Coefficients{2,1}),"b-","LineWidth",2,"DisplayName","Male Fit") ;
% % % % [~,yCI]=predict(mdl,[min(age(logical(sex))):1:max(age(logical(sex)))]',"Prediction","curve");
% % % [~,yCI]=predict(mdl,(-6000:1:4000)',"Prediction","curve");
% % % % [~,yci_obs]=predict(mdl,x,"Prediction","observation");
% % % % plot([min(age(logical(sex))):1:max(age(logical(sex)))]',yCI,"b--","DisplayName","Male CI") ;
% % % plot((-6000:1:4000)',yCI,"b--","DisplayName","Male CI") ;
% % % hold off
% % % ylim([-1 1])
% % % title('EMG xcorr')
% % % xlabel('<Rostral  Caudal>')


h(1) = figure('Color','White','Name','Supplementary Figure 16a','NumberTitle','off');
% subplot(2,2,1)
maxMeanVal = max(abs([meanXLoc meanYLoc cIntFillPtsXLoc cIntFillPtsYLoc]));
plot(timeVecX,meanXLoc,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsXLoc,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
for n = 1:size(xCorrXLoc,1)
    plot(timeVecX,xCorrXLoc(n,:),'Color',[1,0,0,0.1])
end
peak = max(meanXLoc);
peakTime = timeVecX(meanXLoc == peak);
plot(peakTime,peak,'kx','MarkerSize',25)
peakVals = [];
for n = 1:size(xCorrXLoc,1)
    peakVals(end+1) = xCorrXLoc(n,timeVecX == peakTime);
end
stdPeak = std(peakVals);
text(-12,-.8,['Mean Peak:' num2str(peakTime) ',' num2str(peak) ', \pm\sigma:' num2str(stdPeak)])
% plot(timeVecX(brainMotionStart),meanXLoc(brainMotionStart),'rx')
hold off
title('Lateral Brain Motion and Locomotion Cross-Correlation')
ylabel('Noramlized Cross-Correlation')
xlabel('Lags (s)')
ylim([-1 1])

h(1) = figure('Color','White','Name','Supplementary Figure 16b','NumberTitle','off');
% subplot(2,2,2)
maxMeanVal = max(abs([meanXLoc meanYLoc cIntFillPtsXLoc cIntFillPtsYLoc]));
plot(timeVecY,meanYLoc,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsYLoc,[.5 .5 .5],'Linestyle','none');
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
text(-12,-.8,['Mean Peak:' num2str(peakTime) ',' num2str(peak) ', \pm\sigma:' num2str(stdPeak)])
% plot(timeVecX(brainMotionStart),meanXLoc(brainMotionStart),'rx')
hold off
title('Rostral Brain Motion and Locomotion Cross-Correlation')
ylabel('Noramlized Cross-Correlation')
xlabel('Lags (s)')
ylim([-1 1])

h(1) = figure('Color','White','Name','Supplementary Figure 16c','NumberTitle','off');
% subplot(2,2,3)
maxMeanVal = max(abs([meanXEMG meanYEMG cIntFillPtsXEMG cIntFillPtsYEMG]));
plot(timeVecX,meanXEMG,'k')
hold on
f = fill([timeVecX flip(timeVecX)],cIntFillPtsXEMG,[.5 .5 .5],'Linestyle','none');
set(f,'facea',[.2]);
for n = 1:size(xCorrXEMG,1)
    plot(timeVecX,xCorrXEMG(n,:),'Color',[1,0,0,0.1])
end
peak = max(meanXEMG);
peakTime = timeVecX(meanXEMG == peak);
plot(peakTime,peak,'kx','MarkerSize',25)
peakVals = [];
for n = 1:size(xCorrXEMG,1)
    peakVals(end+1) = xCorrXEMG(n,timeVecX == peakTime);
end
stdPeak = std(peakVals);
text(-12,-.8,['Mean Peak:' num2str(peakTime) ',' num2str(peak) ', \pm\sigma:' num2str(stdPeak)])
% plot(timeVecX(brainMotionStart),meanXEMG(brainMotionStart),'rx')
hold off
title('Lateral Brain Motion and EMG Cross-Correlation')
ylabel('Noramlized Cross-Correlation')
xlabel('Lags (s)')
ylim([-1 1])

h(1) = figure('Color','White','Name','Supplementary Figure 16d','NumberTitle','off');
% subplot(2,2,4)
maxMeanVal = max(abs([meanYEMG meanYEMG cIntFillPtsYEMG cIntFillPtsYEMG]));
plot(timeVecY,meanYEMG,'k')
hold on
f = fill([timeVecY flip(timeVecY)],cIntFillPtsYEMG,[.5 .5 .5],'Linestyle','none');
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
text(-12,-.8,['Mean Peak:' num2str(peakTime) ',' num2str(peak) ', \pm\sigma:' num2str(stdPeak)])
% plot(timeVecY(brainMotionStart),meanYEMG(brainMotionStart),'rx')
hold off
title('Rostral Brain Motion and EMG Cross-Correlation')
ylabel('Noramlized Cross-Correlation')
xlabel('Lags (s)')
ylim([-1 1])
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

function hArrow = drawArrowQuiver_FS(p0,p1,color)
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
POld = P;
multFact = 200;
P = Length*P;   %Scale
% P(1,3) = POld(1,3)*multFact;
% P(1,5) = POld(1,5)*multFact;
% P(1,2) = P(1,2)*3;
% P(1,6) = P(1,6)*3;
% P(2,1) = POld(2,1)*multFact;
% P(2,2) = POld(2,2)*multFact;
% P(2,3) = POld(2,3)*multFact;
% P(2,5) = POld(2,5)*multFact;
% P(2,6) = POld(2,6)*multFact;
% P(2,7) = POld(2,7)*multFact;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pdfXVals,pdfYVals] = findKernelPDF(distVals,binEdges)
if size(distVals,2) > 1
    distVals = distVals';
end
probDistObj = fitdist(distVals,'Kernel');
pdfXVals = binEdges(1):.05:binEdges(end);
pdfYVals = pdf(probDistObj,pdfXVals);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = plotSpread_FS(varargin)
%PLOTSPREAD plots distributions of points by spreading them around the y-axis
%
% SYNOPSIS: handles = plotSpread(data, propertyName, propertyValue, ...)
%           handles = plotSpread(ah, ...
%           deprecated:
%           handles = plotSpread(data,binWidth,spreadFcn,xNames,showMM,xValues)
%
% INPUT data: cell array of distributions or nDatapoints-by-mDistributions
%           array, or array with data that is indexed by either
%           distributionIdx or categoryIdx, or both.
%       distributionIdx: grouping variable that determines to which
%           distribution a data point belongs. Grouping is
%           resolved by calling grp2idx, and unless xNames have
%           been supplied, group names determine the x-labels.
%           If the grouping variable is numeric, group labels also
%           determine x-values, unless the parameter xValues has
%           been specified.
%       distributionColors : color identifier (string, cell array of
%           strings), or colormap, with a single color, or one color per
%           distribution (or per entry in distributionIdx). Colors the
%           distributions. Default: 'b'
%       distributionMarkers : string, or cell array of strings, with either
%           a single marker or one marker per distribution (or per entry in
%           distributionIdx). See linespec for admissible markers.
%           Default: '.'
%		categoryIdx: grouping variable that determines group membership for data
%			points across distributions. Grouping is resolved by calling
%           grp2idx.
%       categoryColors : color identifier (cell array of
%           strings), or colormap, with one color per category.
%           Colors the categories, and will override distributionColors.
%           Default is generated using distinguishable_colors by Timothy E.
%           Holy.
%       categoryMarkers : cell array of strings, with one marker per
%           category. See linespec for admissible markers. Will override
%           distributionMarkers. Default: ''
%       categoryLabels : cell array of strings with one label per category
%           (categories sorted in ascending order). Default: unique
%           category indices
%       binWidth : width of bins (along y) that control which data
%           points are considered close enough to be spread. Default: 0.1
%       spreadFcn : cell array of length 2 with {name,param}
%           if name is 'lin', the spread goes linear with the number of
%             points inside the bin, until it reaches the maximum of 0.9 at
%             n==param.
%           if name is 'xp', the spread increases as 1-exp(log(0.9)*x).
%             param is empty
%           Default {'xp',[]}
%       spreadWidth : width, along the x-axis (y-axis if flipped) that can
%           at most be covered by the points. Default:
%           median(diff(sort(xValues))); 1 if no xValues have been supplied
%       showMM : if 1, mean and median are shown as red crosses and
%                green squares, respectively. Default: 0
%                2: only mean
%                3: only median
%                4: mean +/- standard error of the mean (no median)
%                5: mean +/- standard deviation (no median)
%       xNames : cell array of length nDistributions containing x-tick names
%               (instead of the default '1,2,3')
%       xValues : list of x-values at which the data should
%                 be plotted. Default: 1,2,3...
%       xMode  : if 'auto', x-ticks are spaced automatically. If 'manual',
%                there is a tick for each distribution. If xNames is
%                provided as input, xMode is forced to 'manual'. Default:
%                'manual'.
%       xyOri  : orientation of axes. Either 'normal' (=default), or
%                'flipped'. If 'flipped', the x-and y-axes are switched, so
%                that violin plots are horizontal. Consequently,
%                axes-specific properties, such as 'yLabel' are applied to
%                the other axis.
%       yLabel : string with label for y-axis. Default : ''
%       ah  : handles of axes into which to plot
%
% OUTPUT handles: 3-by-1 cell array with handles to distributions,
%          mean/median etc, and the axes, respectively
%
% REMARKS: plotSpread is useful for distributions with a small number of
%          data points. For larger amounts of data, distributionPlot is
%          more suited.
%
% EXAMPLES: data = {randn(25,1),randn(100,1),randn(300,1)};
%           figure,plotSpread(data,[],[],{'25 pts','100 pts','300 pts'})
%
%            data = [randn(50,1);randn(50,1)+3.5]*[1 1];
%            catIdx = [ones(50,1);zeros(50,1);randi([0,1],[100,1])];
%            figure
%            plotSpread(data,'categoryIdx',catIdx,...
%                 'categoryMarkers',{'o','+'},'categoryColors',{'r','b'})
%
% END
%
% created with MATLAB ver.: 7.9.0.3470 (R2009b) on Mac OS X  Version: 10.5.7 Build: 9J61
%
% created by: jonas
% DATE: 11-Jul-2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def.binWidth = 0.1;
def.spreadFcn = {'xp',[]};
def.xNames = [];
def.showMM = false;
def.xValues = [];
def.distributionIdx = [];
def.distributionColors = 'b';
def.distributionMarkers = '.';
def.xMode = 'manual';
def.xyOri = 'normal';
def.categoryIdx = [];
def.categoryColors = [];
def.categoryMarkers = '';
def.categoryLabels = '';
def.yLabel = '';
def.spreadWidth = [];

% in development
def.individualLabels = false; % one category label across all distributions
%                               this should be smartly determined rather
%                               than hard-coded

%% CHECK INPUT

% check for axes handle
if ~iscell(varargin{1}) && length(varargin{1}) == 1 && ...
        ishandle(varargin{1}) && strcmp(get(varargin{1},'Type'),'axes')
    ah = varargin{1};
    data = varargin{2};
    varargin(1:2) = [];
    newAx = false;
else
    ah = gca;
    data = varargin{1};
    varargin(1) = [];
    % if the axes have children, it's not new (important for adjusting
    % limits below)
    newAx = isempty(get(ah,'Children'));
end

% optional arguments
parserObj = inputParser;
parserObj.FunctionName = 'plotSpread';
distributionIdx = [];distributionLabels = '';
if ~isempty(varargin) && ~ischar(varargin{1}) && ~isstruct(varargin{1})
    % old syntax
    parserObj.addOptional('binWidth',def.binWidth);
    parserObj.addOptional('spreadFcn',def.spreadFcn);
    parserObj.addOptional('xNames',def.xNames);
    parserObj.addOptional('showMM',def.showMM);
    parserObj.addOptional('xValues',def.xValues);
    
    parserObj.parse(varargin{:});
    opt = parserObj.Results;
    
    opt.distributionIdx = [];
    opt.distributionColors = def.distributionColors;
    opt.distributionMarkers = def.distributionMarkers;
    opt.xMode = def.xMode;
    opt.xyOri = def.xyOri;
    opt.categoryIdx = [];
    opt.categoryColors = def.distributionColors;
    opt.categoryMarkers = def.distributionMarkers;
    opt.categoryLabels = def.categoryLabels;
    opt.yLabel = '';
    opt.spreadWidth = def.spreadWidth;
    opt.individualLabels = false;
    
    for fn = fieldnames(def)'
        if ~isfield(opt,fn{1})
            % Manually adding the new defaults means a lot fewer bugs
            error('please add option %s to old syntax',fn{1});
        end
        if isempty(opt.(fn{1}))
            opt.(fn{1}) = def.(fn{1});
        end
    end
    
else
    % new syntax
    defNames = fieldnames(def);
    for dn = defNames(:)'
        parserObj.addParamValue(dn{1},def.(dn{1}));
    end
    
    
    parserObj.parse(varargin{:});
    opt = parserObj.Results;
end

% We want data to be a vector, so that indexing with both groupIdx and
% distributionIdx becomes straightforward, and so that we can conveniently
% eliminate NaNs that otherwise could mess up grouping.
% Consequently, if data is a cell array, we convert it, and build a
% corresponding distributionIdx (allowing a user-supplied distributionIdx
% to override, though), and then we go and take care of groupIdx. Once all
% three indices have been built, NaN can be removed.

if iscell(data)
    % make sure data is all n-by-1
    data = cellfun(@(x)x(:),data,'UniformOutput',false);
    nData = length(data);
    nn = cellfun(@numel,data);
    % make vector
    data = cat(1,data{:});
    distributionIdx = repeatEntries_FS((1:nData)',nn);
else
    % distributions in columns
    nData = size(data,2);
    distributionIdx = repeatEntries_FS((1:nData)',size(data,1));
    data = data(:);
end



% distribution groups
if ~isempty(opt.distributionIdx)
    [distributionIdx,distributionLabels,vals] = grp2idx(opt.distributionIdx);
    % convert data to cell array
    nData = length(distributionLabels);
    % if not otherwise provided, use group labels for xnames
    if isempty(opt.xNames)
        opt.xNames = distributionLabels;
        if ~iscell(opt.xNames)
            opt.xNames = num2cell(opt.xNames);
        end
    end
    if isnumeric(vals) && isempty(opt.xValues)
        opt.xValues = vals;
    end
end

if ~isempty(opt.xNames)
    opt.xMode = 'manual';
end


% distribution colors&markers
if ischar(opt.distributionColors)
    opt.distributionColors = {opt.distributionColors};
end
if iscell(opt.distributionColors)
    if length(opt.distributionColors) == 1
        % expand
        opt.distributionColors = repmat(opt.distributionColors,nData,1);
    elseif length(opt.distributionColors) ~= nData
        error('please submit one color per distribution (%i dist, %i colors)',nData,length(opt.distributionColors));
    end
    
else
    if size(opt.distributionColors,2) ~= 3
        error('please specify colormap with three columns')
    end
    if size(opt.distributionColors,1) == 1
        opt.distributionColors = repmat(opt.distributionColors,nData,1);
    elseif size(opt.distributionColors,1) ~= nData
        error('please submit one color per distribution (%i dist, %i colors)',nData,size(opt.distributionColors,1));
    end
    
    % create a cell array
    opt.distributionColors = mat2cell(opt.distributionColors,ones(nData,1),3);
end

if ischar(opt.distributionMarkers)
    opt.distributionMarkers = {opt.distributionMarkers};
end
if length(opt.distributionMarkers) == 1
    % expand
    opt.distributionMarkers = repmat(opt.distributionMarkers,nData,1);
elseif length(opt.distributionMarkers) ~= nData
    error('please submit one color per distribution (%i dist, %i colors)',nData,length(opt.distributionMarkers));
end


stdWidth = 1;
if isempty(opt.xValues)
    opt.xValues = 1:nData;
end


if isempty(opt.spreadWidth) 
    % scale width
    tmp = median(diff(sort(opt.xValues)));
    if ~isnan(tmp)
        stdWidth = tmp;
    end
else
    stdWidth = opt.spreadWidth;
end

if ~ischar(opt.xyOri) || ~any(ismember(opt.xyOri,{'normal','flipped'}))
    error('option xyOri must be either ''normal'' or ''flipped'' (is ''%s'')',opt.xyOri);
end


% check for categoryIdx/colors/markers
% If there are categories, check colors/markers individually first,
% then check whether any of them at all have been supplied, and
% if not, override distributionColors with default categoryColors

if isempty(opt.categoryIdx)
    categoryIdx = ones(size(distributionIdx));
    nCategories = 1;
    categoryLabels = '';
else
    [categoryIdx,categoryLabels] = grp2idx(opt.categoryIdx(:));
    nCategories = max(categoryIdx);
end
if ~isempty(opt.categoryLabels)
    categoryLabels = opt.categoryLabels;
elseif ~iscell(categoryLabels)
    categoryLabels = num2cell(categoryLabels);
end

% plotColors, plotMarkers, plotLabels: nDist-by-nCat arrays
plotColors = repmat(opt.distributionColors(:),1,nCategories);
plotMarkers= repmat(opt.distributionMarkers(:),1,nCategories);

if isempty(distributionLabels)
    distributionLabels = opt.xNames;
    if isempty(distributionLabels)
        distributionLabels = cellstr(num2str(opt.xValues(:)));
    end
end

if nCategories == 1
    plotLabels = distributionLabels(:);
else
    plotLabels = cell(nData,nCategories);
    for iData = 1:nData
        for iCategory = 1:nCategories
            if opt.individualLabels
            plotLabels{iData,iCategory} = ...
                sprintf('%s-%s',num2str(distributionLabels{iData}),...
                num2str(categoryLabels{iCategory}));
            else
                plotLabels{iData,iCategory} = ...
                sprintf('%s',...
                num2str(categoryLabels{iCategory}));
            end
        end
    end
    
end




categoryIsLabeled = false;
if nCategories > 1
    % if not using defaults for categoryColors: apply them
    if ~any(strcmp('categoryColors',parserObj.UsingDefaults))
        if iscell(opt.categoryColors)
            if length(opt.categoryColors) ~= nCategories
                error('please supply one category color per category')
            end
            plotColors = repmat(opt.categoryColors(:)',nData,1);
            categoryIsLabeled = true;
        else
            if all(size(opt.categoryColors) ~= [nCategories,3])
                error('please supply a #-of-categories-by-3 color array')
            end
            plotColors = repmat( mat2cell(opt.categoryColors,ones(nCategories,1),3)', nData,1);
            categoryIsLabeled = true;
        end
    end
    
    if ~any(strcmp('categoryMarkers',parserObj.UsingDefaults))
        if length(opt.categoryMarkers) ~= nCategories
            error('please supply one category marker per category')
        end
        if ~iscell(opt.categoryMarkers)
            error('please supply a list of markers as cell array')
        end
        plotMarkers = repmat(opt.categoryMarkers(:)',nData,1);
        categoryIsLabeled = true;
    end
    
    if ~categoryIsLabeled
        % use distinguishable_colors to mark categories
        
        plotColors = repmat( mat2cell(...
            distinguishable_colors(nCategories),...
            ones(nCategories,1),3)', nData,1);
        
    end
    
end


% remove NaNs from data
badData = ~isfinite(data) | ~isfinite(distributionIdx) | ~isfinite(categoryIdx);
data(badData) = [];
distributionIdx(badData) = [];
categoryIdx(badData) = [];




%% TRANSFORM DATA
% Here, I try to estimate what the aspect ratio of the data is going to be
fh = figure('Visible','off');
if ~isempty(data)
    minMax = [min(data);max(data)];
else
    minMax = [0 1];
end
switch opt.xyOri
    case 'normal'
        plot([0.5;nData+0.5],minMax,'o');
    case 'flipped'
        plot(minMax,[0.5;nData+0.5],'o');
        
end
aspectRatio = get(gca,'DataAspectRatio');
close(fh);

tFact = aspectRatio(2)/aspectRatio(1);
if strcmp(opt.xyOri,'flipped')
    tFact = 1/tFact;
end

%% SPREAD POINTS
% assign either nData, or xValues number of values, in case we're working
% with group-indices
[m,md,sem,sd] = deal(nan(max(nData,length(opt.xValues)),1));
% make sure xValues are not something weird
opt.xValues = double(opt.xValues);

    
% augment data to make n-by-2
data(:,2) = 0;
for iData = 1:nData
    currentDataIdx = distributionIdx==iData;
    currentData = data(currentDataIdx,1);
    
    if ~isempty(currentData)
        
        % transform and sort
        currentData = currentData / tFact;
        %currentData = sort(currentData);
        
        % add x
        currentData = [ones(size(currentData))*opt.xValues(iData),currentData]; %#ok<AGROW>
        
        % step through the data in 0.1 increments. If there are multiple
        % entries, spread along x
        for y = min(currentData(:,2)):opt.binWidth:max(currentData(:,2))
            % find values
            valIdx = find(currentData(:,2) >= y & currentData(:,2) < y+opt.binWidth);
            nVal = length(valIdx);
            if nVal > 1
                % spread
                switch opt.spreadFcn{1}
                    case 'xp'
                        spreadWidth = stdWidth*0.9*(1-exp(log(0.9)*(nVal-1)));
                    case 'lin'
                        spreadWidth = stdWidth*0.9*min(nVal-1,opt.spreadFcn{2})/opt.spreadFcn{2};
                end
                spreadDist = spreadWidth / (nVal - 1);
                if isEven_FS(nVal)
                    offset = spreadDist / 2;
                else
                    offset = eps;
                end
                for v = 1:nVal
                    currentData(valIdx(v),1) = opt.xValues(iData) + offset;
                    % update offset
                    offset = offset - sign(offset) * spreadDist * v;
                end
            end
        end
        
        % update data
        currentData(:,2) = data(currentDataIdx,1);
        data(currentDataIdx,:) = currentData;
        
        
        if opt.showMM > 0
            m(iData) = nanmean(currentData(:,2));
            md(iData) = nanmedian(currentData(:,2));
            sd(iData) = nanstd(currentData(:,2));
            sem(iData) = sd(iData)/sqrt(sum(isfinite(currentData(:,2))));
        end
    end % test isempty
end


%% plot
set(ah,'NextPlot','add')
ph = NaN(nData,nCategories);
for iData = 1:nData
    for iCategory = 1:nCategories
        currentIdx = distributionIdx == iData & categoryIdx == iCategory;
        if any(currentIdx)
            switch opt.xyOri
                case 'normal'
                    ph(iData,iCategory) = plot(ah,data(currentIdx,1),...
                        data(currentIdx,2),...
                        'marker',plotMarkers{iData,iCategory},...
                        'color',plotColors{iData,iCategory},...
                        'lineStyle','none',...
                        'DisplayName',plotLabels{iData,iCategory});
                case 'flipped'
                    ph(iData,iCategory) = plot(ah,data(currentIdx,2),...
                        data(currentIdx,1),...
                        'marker',plotMarkers{iData,iCategory},...
                        'color',plotColors{iData,iCategory},...
                        'lineStyle','none',...
                        'DisplayName',plotLabels{iData,iCategory});
            end
        end
    end
end



% if ~empty, use xNames
switch opt.xyOri
    case 'normal'
        switch opt.xMode
            case 'manual'
                set(ah,'XTick',opt.xValues);
                if ~isempty(opt.xNames)
                    set(ah,'XTickLabel',opt.xNames)
                end
            case 'auto'
                % no need to do anything
        end
        
        % have plot start/end properly
        minX = min(opt.xValues)-stdWidth;
        maxX = max(opt.xValues)+stdWidth;
        if ~newAx
            oldLim = xlim;
            minX = min(minX,oldLim(1));
            maxX = max(maxX,oldLim(2));
        end
        xlim([minX,maxX])
        
        ylabel(ah,opt.yLabel)
        
    case 'flipped'
        switch opt.xMode
            case 'manual'
                set(ah,'YTick',opt.xValues);
                if ~isempty(opt.xNames)
                    set(ah,'YTickLabel',opt.xNames)
                end
            case 'auto'
                % no need to do anything
        end
        
        % have plot start/end properly (for ease of copying, only switch
        % xlim to ylim
        minX = min(opt.xValues)-stdWidth;
        maxX = max(opt.xValues)+stdWidth;
        if ~newAx
            oldLim = ylim;
            minX = min(minX,oldLim(1));
            maxX = max(maxX,oldLim(2));
        end
        ylim([minX,maxX])
        
        xlabel(ah,opt.yLabel);
        
end

% ## in development
if ~opt.individualLabels
       % hack: add legend entry only once per category
       goodH = ishandle(ph);
       for iCategory = 1:nCategories
           for iData = find(goodH(:,iCategory),1,'first')+1:nData
       if goodH(iData,iCategory)
           set(get(get(ph(iData,iCategory),'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off');
       end
           end
       end
       
end


% add mean/median
mh = [];mdh=[];
if opt.showMM
    % plot mean, median. Mean is filled red circle, median is green square
    % I don't know of a very clever way to flip xy and keep everything
    % readable, thus it'll be copy-paste
    switch opt.xyOri
        case 'normal'
            if any(opt.showMM==[1,2])
                mh = plot(ah,opt.xValues,m,'+r','Color','r','MarkerSize',12);
            end
            if any(opt.showMM==[1,3])
                mdh = plot(ah,opt.xValues,md,'sg','MarkerSize',12);
            end
            if opt.showMM == 4
                mh = plot(ah,opt.xValues,m,'+r','Color','r','MarkerSize',12);
                mdh = myErrorbar_FS(ah,opt.xValues,m,sem);
            end
            if opt.showMM == 5
                mh = plot(ah,opt.xValues,m,'+r','Color','r','MarkerSize',12);
                mdh = myErrorbar_FS(ah,opt.xValues,m,sd);
            end
        case 'flipped'
            if any(opt.showMM==[1,2])
                mh = plot(ah,m,opt.xValues,'+r','Color','r','MarkerSize',12);
            end
            if any(opt.showMM==[1,3])
                mdh = plot(ah,md,opt.xValues,'sg','MarkerSize',12);
            end
            if opt.showMM == 4
                mh = plot(ah,m,opt.xValues,'+r','Color','r','MarkerSize',12);
                mdh = myErrorbar_FS(ah,m,opt.xValues,[sem,NaN(size(sem))]);
            end
            if opt.showMM == 5
                mh = plot(ah,m,opt.xValues,'+r','Color','r','MarkerSize',12);
                mdh = myErrorbar_FS(ah,m,opt.xValues,[sd,NaN(size(sd))]);
            end
    end
end

%==========================
%% CLEANUP & ASSIGN OUTPUT
%==========================

if nargout > 0
    handles{1} = ph;
    handles{2} = [mh;mdh];
    handles{3} = ah;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hh = myErrorbar_FS(varargin)
%MYERRORBAR Adds errorbars to existing plot (unlike errorbar.m, which creates a new plot, and allows only bars for y values)
%   MYERRORBAR(X,Y,L,U) adds error bars to the graph of vector X vs. vector Y with
%   error bars specified by the vectors L and U.  L and U contain the
%   lower and upper error ranges for each point in Y.  Each error bar
%   is L(i) + U(i) long and is drawn a distance of U(i) above and L(i)
%   below the points in (X,Y). If X,Y,L and U are matrices then each column
%   produces a separate line.
%   If L,U are the same size as X, Y, only error bars for Y will be plotted.
%   If L,U are twice the size of X,Y (or have twice the number of columns for
%   matrices), the first half of L, U specifies error bar lengths for X and the
%   second half specifies error bars for Y
%
%   MYERRORBAR(X,Y,E) or MYERRORBAR(Y,E) plots error bars [Y-E Y+E].
%
%   MYERRORBAR(AX,...), where AX is an axis handle, plots errorbars into
%                       axes AX
%
%   H = MYERRORBAR(...) returns a vector of line handles.
%
%   The tag of the errorbar-lines is: errorBar
%
%   For example,
%      x = 1:10;
%      y = sin(x);
%      e = std(y)*ones(size(x));
%      myErrorbar(x,y,e)
%   draws symmetric error bars of unit standard deviation for y values.
%      myErrorbar(x,y,[e,e])
%   draws symmetric error bars of unit standard deviation for x and y
%   values.
%
%   Based on the matlab-function errorbar as revised by Claude Berney
%   c: jonas, 06-03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==================
% check input
%==================

if nargin < 2
    error('not enough input arguments!')
end

% check if the first input argument is a handle
if length(varargin{1}) == 1 && ishandle(varargin{1}) && strcmpi(get(varargin{1},'Type'),'axes')
    axesH = varargin{1};
    % remove axis handle
    varargin(1) = [];
else
    axesH = gca;
end

% there could be
% y,e
% x,y,e
% x,y,l,u

switch length(varargin)
    case 2
        % y, e
        y = varargin{1};
        y = y(:);
        lengthY = length(y);
        x = [1:lengthY]';
        
        e = varargin{2};
        % check for 2 dimension errorbars
        e = e(:);
        if length(e) == 2*lengthY
            e = reshape(e,lengthY,2);
        end
        [l,u] = deal(e);
        
    case 3
        % x,y,e
        x = varargin{1};
        x = x(:);
        y = varargin{2};
        y = y(:);
        lengthY = length(y);
        
        e = varargin{3};
        % check for 2 dimension errorbars
        e = e(:);
        if length(e) == 2*lengthY
            e = reshape(e,lengthY,2);
        end
        [l,u] = deal(e);
        
    case 4
        % x,y,l,u
        % x,y,e
        x = varargin{1};
        x = x(:);
        y = varargin{2};
        y = y(:);
        lengthY = length(y);
        
        l = varargin{3};
        % check for 2 dimension errorbars
        l = l(:);
        if length(l) == 2*lengthY
            l = reshape(l,lengthY,2);
        end
        u = varargin{4};
        % check for 2 dimension errorbars
        u = u(:);
        if length(u) == 2*lengthY
            u = reshape(u,lengthY,2);
        end
        
        if ~all(size(u)==size(l))
            error('l, u have to be the same size!')
        end
        
end % switch number of inputs


u = abs(u);
l = abs(l);

if ischar(x) || ischar(y) || ischar(u) || ischar(l)
    error('Arguments must be numeric.')
end

if ~isequal(size(x),size(y))
    error('The sizes of X and Y must be the same.');
end

if isequal([1 2].*size(x),size(l)) && isequal([1 2].*size(x),size(u))
    xyBars = 1;
elseif isequal(size(x),size(l)) && isequal(size(x),size(u))
    xyBars = 0;
else
    error('The sizes of L and U must be equal to or twice the size of X, Y')
end

%=======================


% Plot graph and bars
hold_state = ishold;
hold on;


%find color of current plot
dataH = get(axesH,'Children');
myLineH = dataH(1);
% support also bar plots
if strcmp(get(myLineH,'Type'),'hggroup')
    latestColor = get(myLineH,'EdgeColor'); %new children are added on top!
else
    latestColor = get(myLineH,'Color'); %new children are added on top!
end

tee=0;
if ~strcmp('log',get(axesH,'XScale'))
    tee = (max(x(:))-min(x(:)))/100;  % make tee .02 x-distance for error bars
    tee = min(tee,0.3*nanmedian(diff(unique(x(:))))); % or at most 0.3*deltaX
    xl = x - tee;
    xr = x + tee;
end
if strcmp('log',get(axesH,'XScale'))
    tee = (max(log(x(:)))-min(log(x(:))))/100;  % make tee .02 x-distance for error bars
    tee = min(tee,0.3*nanmedian(diff(unique(log(x(:)))))); % or at most 0.3*deltaX
    
    xl = x *exp(tee);
    xr = x *exp(-tee);
end

if xyBars
    if ~strcmp('log',get(axesH,'YScale'))
        tee = (max(y(:))-min(y(:)))/100;  % make tee .02 y-distance for error bars
        tee = min(tee,0.3*nanmedian(diff(unique(y(:))))); % or at most 0.3*deltaY
        
        yl = y - tee;
        yr = y + tee;
    end
    if strcmp('log',get(axesH,'YScale'))
        tee = (max(log(y(:)))-min(log(y(:))))/100;  % make tee .02 y-distance for error bars
        tee = min(tee,0.3*nanmedian(diff(unique(log(y(:)))))); % or at most 0.3*deltaX
        
        yl = y *exp(tee);
        yr = y *exp(-tee);
    end
end

%specify coordinates to plot error bars
if xyBars
    xtop = x + u(:,1:size(x,2));
    xbot = x - l(:,1:size(x,2));
    ytop = y + u(:,size(x,2)+1:end);
    ybot = y - l(:,size(x,2)+1:end);
else
    ytop = y + u;
    ybot = y - l;
end
n = size(y,2);

% build up nan-separated vector for bars
xb = zeros(lengthY*9,n);
xb(1:9:end,:) = x;
xb(2:9:end,:) = x;
xb(3:9:end,:) = NaN;
xb(4:9:end,:) = xl;
xb(5:9:end,:) = xr;
xb(6:9:end,:) = NaN;
xb(7:9:end,:) = xl;
xb(8:9:end,:) = xr;
xb(9:9:end,:) = NaN;

yb = zeros(lengthY*9,n);
yb(1:9:end,:) = ytop;
yb(2:9:end,:) = ybot;
yb(3:9:end,:) = NaN;
yb(4:9:end,:) = ytop;
yb(5:9:end,:) = ytop;
yb(6:9:end,:) = NaN;
yb(7:9:end,:) = ybot;
yb(8:9:end,:) = ybot;
yb(9:9:end,:) = NaN;

h = [line(xb,yb,'parent',axesH,'Color',latestColor)];

if xyBars
    
    xb(1:9:end,:) = xtop;
    xb(2:9:end,:) = xbot;
    xb(3:9:end,:) = NaN;
    xb(4:9:end,:) = xtop;
    xb(5:9:end,:) = xtop;
    xb(6:9:end,:) = NaN;
    xb(7:9:end,:) = xbot;
    xb(8:9:end,:) = xbot;
    xb(9:9:end,:) = NaN;
    
    yb(1:9:end,:) = y;
    yb(2:9:end,:) = y;
    yb(3:9:end,:) = NaN;
    yb(4:9:end,:) = yl;
    yb(5:9:end,:) = yr;
    yb(6:9:end,:) = NaN;
    yb(7:9:end,:) = yl;
    yb(8:9:end,:) = yr;
    yb(9:9:end,:) = NaN;
    
    h = [h;line(xb,yb,'parent',axesH,'Color',latestColor)];
    
end

%set the tag of all errorBar-objects to 'errorBar'
set(h,'Tag','errorBar');

% make sure errorbar doesn't produce a legend entry
for lineH = h'
    set(get(get(lineH,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off');
end


if ~hold_state, hold off; end

if nargout>0, hh = h; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = repeatEntries_FS(val,kTimes)
%REPEATENTRIES fills a matrix with k repeats the rows of the input matrix
%
% SYNOPSIS out = repeatEntries(val,kTimes)
%
% INPUT    val    : matrix (or vectors) containing the rows to repeat (works for strings, too)
%          kTimes : number of repeats of each row (scalar or vector of size(vlaues,1))
%
% OUTPUT   out    : matrix of size [sum(kTimes) size(values,2)] containing
%                   repeated entries specified with k
%
% EXAMPLES     repeatEntries([1;2;3;4],[2;3;1;1]) returns [1;1;2;2;2;3;4]
%
%              repeatEntries([1;2;3;4],2) returns [1;1;2;2;3;3;4;4]
%
% c: jonas, 2/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% note: in case we need to speed this up: adapt the code below
% nn = cellfun(@numel,points);
% a = find(nn);
% index = zeros(sum(nn),1);
% index([1;cumsum(nn(a(1:end-1)))+1])=1;
% 
% % get the indices
% ii = a(cumsum(index));

%===========
% test input
%===========

% nargin
if nargin ~= 2 || isempty(val) || isempty(kTimes)
    error('two non-empty input arguments are needed!')
end

% size
valSize = size(val);
if length(valSize)>2
    error('only 2D arrays supported for val')
end



% decide whether we have scalar k
numK = length(kTimes);
if numK == 1
    scalarK = 1;
elseif numK ~= valSize(1)
    error('vector k must have the same length as the number of rows in val or be a scalar')
else
    % check again whether we could use scalar k
    if all(kTimes(1) == kTimes)
        scalarK = 1;
        kTimes = kTimes(1);
    else
        scalarK = 0;
    end
end

% do not care about size of k: we want to make a col vector out of it - and
% this vector should only contain nonzero positive integers
kTimes = round(kTimes(:));
% if there are any negative values or zeros, remove the entry
if scalarK && kTimes < 1
    out = [];
    return
end
if ~scalarK
    badK = kTimes < 1;
    kTimes(badK) = [];
    val(badK,:) = [];
    % update valSize
    valSize = size(val);
    if any(valSize==0)
        out = [];
        return
    end
end
%kTimes = max(kTimes,ones(size(kTimes)));


%============
% fill in out
%============

% first the elegant case: scalar k
if scalarK

    % build repeat index matrix idxMat
    idxMat = meshgrid( 1:valSize(1), 1:kTimes(1) );
    idxMat = idxMat(:); % returns [1;1...2;2;... etc]

    out = val(idxMat,:);

    % second: the loop
else

    % init out, init counter
    if iscell(val)
        out = cell(sum(kTimes) , valSize(2));
    else
    out = zeros( sum(kTimes), valSize(2) );
    end
    endct = 0;

    if valSize(2) == 1

        % vector: fill directly

        % loop and fill
        for i = 1:valSize(1)
            startct = endct + 1;
            endct   = endct + kTimes(i);
            out(startct:endct,:) = val(i);
        end % for i=1:valSize(1)

    else

        % matrix: fill via index list

        idxMat = zeros(sum(kTimes),1);

        for i = 1:valSize(1)
            startct = endct + 1;
            endct   = endct + kTimes(i);
            idxMat(startct:endct) = i;
        end % for i=1:valSize(1)
        out = val(idxMat,:);

    end

    % check for strings and transform if necessary
    if ischar(val)
        out = char(out);
    end

end % if doScalar
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = isEven_FS(in)
%ISEVEN checks whether a number is even
%
% SYNOPSIS out = isEven(in)
%
% INPUT    in :  input (array) of numbers to be tested. 
% OUTPUT   out:  array of size(in) with 
%                   1 for even integers and zero
%                   0 for odd integers
%                 NaN for non-integers
%                out is a logical array as long as the input is all integers.
%
% c: jonas 5/05
% Last modified 11/24/2009 - Jonas Dorn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out = mod(in+1, 2);
% Set NaN for non-integer data, because they are neither odd or even
out((out ~= 0) & (out ~= 1)) = NaN;

% since doubles cannot be used for logical indexing, we should convert to
% logicals if possible. 
if all(isfinite(out(:)))
    out = logical(out);
end
end