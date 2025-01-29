close all

% brain without movement
% fileNames = {'C:\Workspace\210218_005_processed_combined.mat','C:\Workspace\210218_006_processed_combinedMean.mat','C:\Workspace\210218_007_processed_combined.mat','C:\Workspace\210218_008_processed_combined.mat'};
% ssTimes = {[0 10;10 20;20 30;45 55;55 65;80 90],[0 10;30 40;45 55;55 65;70 80;80 90],[0 10;25 35;35 45;65 75;75 85],[40 50;50 60;60 70]};

% brain with movement
% fileNames = {'C:\Workspace\210218_005_processed_combined.mat','C:\Workspace\210218_006_processed_combinedMean.mat','C:\Workspace\210218_007_processed_combined.mat','C:\Workspace\210218_008_processed_combined.mat'};
% ssTimes = {[30 40;67 77],[16 26],[12 22;49 59;83 93],[20 30;83 93]};

% beads without movement
% fileNames = {'C:\Workspace\210301_001_processed_combined.mat','C:\Workspace\210301_002_processed_combined.mat','C:\Workspace\210301_003_processed_combined.mat','C:\Workspace\210301_004_processed_combined.mat','C:\Workspace\210301_005_processed_combined.mat','C:\Workspace\210301_006_processed_combined.mat'};
% ssTimes = {[0 10;10 20;20 30;30 40;40 50;50 60;60 70;70 80;80 90],[0 10;10 20;20 30;30 40;40 50;50 60;60 70;70 80;80 90],[0 10;10 20;20 30;30 40;40 50;50 60;60 70;70 80;80 90],[0 10;35 45;45 55],[0 10;25 35;35 45],[0 10;30 40;40 50]};

% beads with movement
fileNames = {'C:\Workspace\210301_004_processed_combined.mat','C:\Workspace\210301_005_processed_combined.mat','C:\Workspace\210301_006_processed_combined.mat'};
ssTimes = {[10 20;20 30;57 67;83 93],[12 22;47 57;63 73;83 93],[10 20;50 60;83 93]};


locClipsX = [];
locClipsY = [];
velClips = [];
velClipsX = [];
velClipsY = [];
locClipsXRaw = [];
locClipsYRaw = [];
velClipsRaw = [];
velClipsXRaw = [];
velClipsYRaw = [];
for n = 1:numel(fileNames)
    for i = 1:size(ssTimes{n},1)
        load(fileNames{n});
        timeSLoc = (1:size(movementData.targetPosition,1))*movementData.secondsPerFrame;
        locClipX = movementData.targetPosition(ssTimes{n}(i,1) <= timeSLoc & timeSLoc <= ssTimes{n}(i,2),1);
        locClipY = movementData.targetPosition(ssTimes{n}(i,1) <= timeSLoc & timeSLoc <= ssTimes{n}(i,2),2);
        if n > 1 || i > 1
            if length(locClipX) > size(locClipsX,2)
                locClipX = locClipX(1:size(locClipsX,2));
                locClipY = locClipY(1:size(locClipsX,2));
            elseif length(locClipX) < size(locClipsX,2)
                locClipsX = locClipsX(:,1:length(locClipX));
                locClipsY = locClipsY(:,1:length(locClipX));
                locClipsXRaw = locClipsXRaw(:,1:length(locClipX));
                locClipsYRaw = locClipsYRaw(:,1:length(locClipX));
            end
        end
        locClipsX(end+1,:) = detrend(medfilt1(locClipX,6))';
        locClipsY(end+1,:) = detrend(medfilt1(locClipY,6))';
        locClipsXRaw(end+1,:) = medfilt1(locClipX-locClipX(1),6)';
        locClipsYRaw(end+1,:) = medfilt1(locClipY-locClipY(1),6)';
        timeSVel = (1:size(movementData.velocity,1))*movementData.secondsPerFrame;
        velClip = movementData.velocity(ssTimes{n}(i,1) <= timeSVel & timeSVel <= ssTimes{n}(i,2),1);
        velX = (diff(movementData.targetPosition(:,1))./movementData.secondsPerFrame);
        velClipX = velX(ssTimes{n}(i,1) <= timeSVel & timeSVel <= ssTimes{n}(i,2),1);
        velY = (diff(movementData.targetPosition(:,2))./movementData.secondsPerFrame);
        velClipY = velY(ssTimes{n}(i,1) <= timeSVel & timeSVel <= ssTimes{n}(i,2),1);
        if n > 1 || i > 1
            if length(velClip) > size(velClips,2)
                velClip = velClip(1:size(velClips,2));
                velClipX = velClipX(1:size(velClips,2));
                velClipY = velClipY(1:size(velClips,2));
            elseif length(velClip) < size(velClips,2)
                velClips = velClips(:,1:length(velClip));
                velClipsX = velClipsX(:,1:length(velClip));
                velClipsY = velClipsY(:,1:length(velClip));
                velClipsRaw = velClipsRaw(:,1:length(velClip));
                velClipsXRaw = velClipsXRaw(:,1:length(velClip));
                velClipsYRaw = velClipsYRaw(:,1:length(velClip));
            end
        end
        velClips(end+1,:) = detrend(velClip)';
        velClipsX(end+1,:) = detrend(velClipX)';
        velClipsY(end+1,:) = detrend(velClipY)';
        velClipsRaw(end+1,:) = velClip-velClip(1)';
        velClipsXRaw(end+1,:) = velClipX-velClipX(1)';
        velClipsYRaw(end+1,:) = velClipY-velClipY(1)';
        secondsPerFrame = movementData.secondsPerFrame;
    end
    clear movementData
end

dt = secondsPerFrame;
tapers = [2 3];
params.Fs = 1/dt;
params.tapers = tapers;
params.fpass = [0.025 50];
params.err = [2 0.05]; %use jack-knife resampling confidence intervals p = 0.05
% (could also try multitaper frequency domain bootstrapping (MFDB) for less noisy CI - but jackknife used by so many people in NVC)

[PowerX, HzX, ErrorX] = mtspectrumc(locClipsX',params);
[PowerY, HzY, ErrorY] = mtspectrumc(locClipsY',params);
[PowerVel, HzVel, ErrorVel] = mtspectrumc(velClips',params);
[PowerVelX, HzVelX, ErrorVelX] = mtspectrumc(velClipsX',params);
[PowerVelY, HzVelY, ErrorVelY] = mtspectrumc(velClipsY',params);


h(1) = figure('Color','White');
subplot(2,1,1)
semilogy(HzX,PowerX,'k')
title(['\fontsize{20pt}\bf{X Position Frequency Domain, n = ' num2str(size(PowerX,2)) ' clips}'])
xlabel('Frequency (Hz)')
ylabel('Power')
subplot(2,1,2)
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts(PowerX',90);
cIntFillPtsX = removeZeros(cIntFillPtsX);
semilogy(HzX,meanX,'k')
hold on
f = fill([HzX flip(HzX)],cIntFillPtsX,'r','Linestyle','none');
set(f,'facea',[.2]);
hold off
set(gca, 'YScale', 'log')
xlabel('Frequency (Hz)')
ylabel('Power')
title('Mean')

h(2) = figure('Color','White');
subplot(2,1,1)
semilogy(HzY,PowerY,'k')
title(['\fontsize{20pt}\bf{Y Position Frequency Domain, n = ' num2str(size(PowerX,2)) ' clips}'])
xlabel('Frequency (Hz)')
ylabel('Power')
subplot(2,1,2)
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts(PowerY',90);
cIntFillPtsY = removeZeros(cIntFillPtsY);
semilogy(HzY,meanY,'k')
hold on
f = fill([HzY flip(HzY)],cIntFillPtsY,'r','Linestyle','none');
set(f,'facea',[.2]);
hold off
set(gca, 'YScale', 'log')
xlabel('Frequency (Hz)')
ylabel('Power')
title('Mean')

h(3) = figure('Color','White');
subplot(2,1,1)
semilogy(HzVel,PowerVel,'k')
title(['\fontsize{20pt}\bf{Total Velocity Frequency Domain, n = ' num2str(size(PowerX,2)) ' clips}'])
xlabel('Frequency (Hz)')
ylabel('Power')
subplot(2,1,2)
[meanVel,cIntFillPtsVel] = getCIntMeanAndFillPts(PowerVel',90);
cIntFillPtsVel = removeZeros(cIntFillPtsVel);
semilogy(HzVel,meanVel,'k')
hold on
f = fill([HzY flip(HzY)],cIntFillPtsVel,'r','Linestyle','none');
set(f,'facea',[.2]);
hold off
set(gca, 'YScale', 'log')
xlabel('Frequency (Hz)')
ylabel('Power')
title('Mean')

h(4) = figure('Color','White');
subplot(2,1,1)
semilogy(HzVelX,PowerVelX,'k')
title(['\fontsize{20pt}\bf{X Velocity Frequency Domain, n = ' num2str(size(PowerX,2)) ' clips}'])
xlabel('Frequency (Hz)')
ylabel('Power')
subplot(2,1,2)
[meanVelX,cIntFillPtsVelX] = getCIntMeanAndFillPts(PowerVelX',90);
cIntFillPtsVelX = removeZeros(cIntFillPtsVelX);
semilogy(HzVelX,meanVelX,'k')
hold on
f = fill([HzVelX flip(HzVelX)],cIntFillPtsVelX,'r','Linestyle','none');
set(f,'facea',[.2]);
hold off
set(gca, 'YScale', 'log')
xlabel('Frequency (Hz)')
ylabel('Power')
title('Mean')

h(5) = figure('Color','White');
subplot(2,1,1)
semilogy(HzVelY,PowerVelY,'k')
title(['\fontsize{20pt}\bf{Y Velocity Frequency Domain, n = ' num2str(size(PowerX,2)) ' clips}'])
xlabel('Frequency (Hz)')
ylabel('Power')
subplot(2,1,2)
[meanVelY,cIntFillPtsVelY] = getCIntMeanAndFillPts(PowerVelY',90);
cIntFillPtsVelY = removeZeros(cIntFillPtsVelY);
semilogy(HzVelY,meanVelY,'k')
hold on
f = fill([HzVelY flip(HzVelY)],cIntFillPtsVelY,'r','Linestyle','none');
set(f,'facea',[.2]);
hold off
set(gca, 'YScale', 'log')
xlabel('Frequency (Hz)')
ylabel('Power')
title('Mean')

h(6) = figure('Color','White');
subplot(2,1,1)
plot((1:size(locClipsXRaw,2))*secondsPerFrame,locClipsXRaw,'k')
title(['\fontsize{20pt}\bf{X Position, n = ' num2str(size(PowerX,2)) ' clips}'])
xlabel('Time (s)')
ylabel('X Position (\mum)')
xlim([0 10])
subplot(2,1,2)
[meanX,cIntFillPtsX] = getCIntMeanAndFillPts(locClipsXRaw,90);
plot((1:length(meanX))*secondsPerFrame,meanX,'k')
hold on
f = fill([(1:length(meanX))*secondsPerFrame flip((1:length(meanX))*secondsPerFrame)],cIntFillPtsX,'r','Linestyle','none');
set(f,'facea',[.2]);
hold off
xlabel('Time (s)')
ylabel('X Position (\mum)')
title('Mean')
xlim([0 10])

h(7) = figure('Color','White');
subplot(2,1,1)
plot((1:size(locClipsYRaw,2))*secondsPerFrame,locClipsYRaw,'k')
title(['\fontsize{20pt}\bf{Y Position, n = ' num2str(size(PowerX,2)) ' clips}'])
xlabel('Time (s)')
ylabel('Y Position (\mum)')
xlim([0 10])
subplot(2,1,2)
[meanY,cIntFillPtsY] = getCIntMeanAndFillPts(locClipsYRaw,90);
plot((1:length(meanY))*secondsPerFrame,meanY,'k')
hold on
f = fill([(1:length(meanY))*secondsPerFrame flip((1:length(meanY))*secondsPerFrame)],cIntFillPtsY,'r','Linestyle','none');
set(f,'facea',[.2]);
hold off
xlabel('Time (s)')
ylabel('Y Position (\mum)')
title('Mean')
xlim([0 10])

clear

function out = removeZeros(in)

negativePoints = in <= 0;
inMax = in;
inMax(negativePoints) = max(in);
minIn = min(inMax);
out = in;
out(negativePoints) = minIn;
end