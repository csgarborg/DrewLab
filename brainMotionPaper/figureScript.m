% Script to generate all figures for paper
clear;
close all;

%% figure 1e    210218_006

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
% rawData = load('211216_002.txt');
% plotThermoCoherence_FS(combinedMovementDataBrain_2c,combinedMovementDataSkull_2c,rawData)

%% figure 2d

combinedMovementDataBrain_2d = combinedMovementDataBrain_2c;
combinedMovementDataSkull_2d = combinedMovementDataSkull_2c;
plotLocomotionXCorr_FS(combinedMovementDataBrain_2d,combinedMovementDataSkull_2d)

%% figure 3c

combinedMovementDataBrain_3c = combinedMovementDataBrain_2d;
combinedMovementDataSkull_3c = combinedMovementDataSkull_2d;
% plotMovementBrainInSkull_FS(combinedMovementDataBrain_3c,combinedMovementDataSkull_3c,rawData)

%% figure 3d

combinedMovementDataBrain_3d = combinedMovementDataBrain_3c;
combinedMovementDataSkull_3d = combinedMovementDataSkull_3c;
plotEMGMovement_FS(combinedMovementDataBrain_3d,combinedMovementDataSkull_3d)

%% figure 3e

combinedMovementDataBrain_3e = combinedMovementDataBrain_3d;
combinedMovementDataSkull_3e = combinedMovementDataSkull_3d;
plotEMGXCorr_FS(combinedMovementDataBrain_3e,combinedMovementDataSkull_3e)

%% figure 3f

plotLocomotionTriggeredAvg
plotEMGTriggeredAvg









%% SUBFUNCTIONS

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

function plotMotionTracking2P2LPCA_FS(movementData,stationaryData)
targetPositionInSkull = combineBrainSkullMovement_FS(movementData,stationaryData);
motionVec = pcaMotionAnalysis_FS(targetPositionInSkull);
h(2) = figure('Color','White');
scatter(targetPositionInSkull(:,1),targetPositionInSkull(:,2),10)
hold on
drawArrow([0;0],[motionVec(1);motionVec(2)]);
hold off
axis equal square
axis([-6 6 -6 6])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title(['\fontsize{20pt}\bf{Position of Brain in Skull with PCA Vector}'])
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
title('Rostral, \mum')
ylabel('Left, \mum')
xlabel('Caudal, \mum')
text(55,55,'Bregma')
text(55,-2500,'Lambda')
text(-2500,2300,'3 \mum')
rectangle('Position',[-2650 2000 800 800])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotThermoCoherence_FS(movementData,stationaryData,rawData)
params.tapers = [10,19];   % Tapers [n, 2n - 1]
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
Thermo = rawData(:,[1 3]); % load raw thermo data
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
[Cx,~,~,~,~,fx,confCx,~,cErrx] = coherencyc_eLife2020(resp(:,2),brainX(:,2),params);
[Cy,~,~,~,~,fy,confCy,~,cErry] = coherencyc_eLife2020(resp(:,2),brainY(:,2),params);
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
title('Brain Motion Spectrum (Rostral/Caudal)')
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

[PowerX, HzX, ErrorX] = mtspectrumc(positionData(:,1),params);
[PowerY, HzY, ErrorY] = mtspectrumc(positionData(:,2),params);
[PowerT, HzT, ErrorT] = mtspectrumc(thermo,params);

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
title('Brain Motiona dn Locomotion Cross-Correlation')
ylabel('Noramlized Cross-Correlation')
xlabel('Lags (s)')
legend({'L/M Brain Motion','R/C Brain Motion'})
ylim([0 1])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMovementBrainInSkull_FS(movementData,stationaryData,rawData)
targetPositionInSkull = combineBrainSkullMovement_FS(movementData,stationaryData);
movementData.secondsPerFrame = movementData.secondsPerFrame/2;
h(6) = figure('Color','White');
subplot(4,1,1)
plot([1:size(targetPositionInSkull,1)]*movementData.secondsPerFrame,targetPositionInSkull(:,1),'r')
title(['\fontsize{20pt}\bf{Position of Brain in Skull}'])
xlabel('Time (s)')
ylabel('X Position (\mum)')
grid on
axis([0 600 0 5])
if movementData.hemisphere == 1
    text(0,5,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(0,0,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
else
    text(0,5,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(0,0,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
end
subplot(4,1,2)
plot([1:size(targetPositionInSkull,1)]*movementData.secondsPerFrame,targetPositionInSkull(:,2),'g')
xlabel('Time (s)')
ylabel('Y Position (\mum)')
grid on
axis([0 600 -4 4])
text(0,4,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(0,-4,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
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
plot(rawData(:,1),rawData(:,3),'b')
axis([0 600 -5 5])
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
subplot(2,1,1)
scatter(movementData.emgDataInterp(1:(end-best_lag),2),movementData.targetPosition(best_lag:end-1,2))
xlabel('Log EMG Power')
ylabel('Brain Shift, \mum')
text(3.5,6,'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
text(3.5,-2,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
ylim([-2 6])
subplot(2,1,2)
hh2=histogram2(movementData.emgDataInterp(1:(end-best_lag),2),movementData.targetPosition(best_lag:end-1,2), emg_bins, pixel_bins,...
    'DisplayStyle','tile','ShowEmptyBins','on');
imagesc(emg_bins, pixel_bins,log(hh2.Values'+1))
axis xy
colormap hot
xlabel('Log EMG Power')
ylabel('Brain Shift, \mum')
ylim([-2 6])
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
title('Brain Motiona and EMG Cross-Correlation')
ylabel('Noramlized Cross-Correlation')
xlabel('Lags (s)')
legend({'L/M Brain Motion','R/C Brain Motion'})
ylim([0 1])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotLocomotionTriggeredAvg
motionEventsLocationsX = [];
motionEventsLocationsY = [];
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
end
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts(motionEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts(motionEventsLocationsY,90);
meanY = -1*meanY;
cIntFillPtsY = -1*cIntFillPtsY;
timeVecX = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(meanX));
timeVecY = linspace(round(locDataCell{1,2}(1,1)-locDataCell{1,2}(1,2)),round(locDataCell{1,2}(1,3)-locDataCell{1,2}(1,2)),length(meanY));

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
hold off
text(3,-3,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(3,3,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
title(['\fontsize{20pt}\bf{Mean Motion During Locomotion Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
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
load('LTADataCell_FS.mat')
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotEMGTriggeredAvg
motionEventsLocationsX = [];
motionEventsLocationsY = [];
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
end
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts(motionEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts(motionEventsLocationsY,90);
meanY = -1*meanY;
cIntFillPtsY = -1*cIntFillPtsY;
timeVecX = linspace(round(EMGDataCell{1,2}(1,1)-EMGDataCell{1,2}(1,2)),round(EMGDataCell{1,2}(1,3)-EMGDataCell{1,2}(1,2)),length(meanX));
timeVecY = linspace(round(EMGDataCell{1,2}(1,1)-EMGDataCell{1,2}(1,2)),round(EMGDataCell{1,2}(1,3)-EMGDataCell{1,2}(1,2)),length(meanY));

h(10) = figure('Color','White');
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
title(['\fontsize{20pt}\bf{Mean Motion During EMG Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
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
load('LTADataCell_FS.mat')
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
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts(stopMotionEventsLocationsX,90);
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts(stopMotionEventsLocationsY,90);
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
end