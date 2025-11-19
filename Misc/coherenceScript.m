close all
params.tapers = [1,1];   % Tapers [n, 2n - 1]
params.pad = 1;
params.Fs = 20;
params.fpass = [0,10];   % Pass band [0, nyquist]
params.trialave = 1;
params.err = [2,0.05];
% you want to make sure 'LH_restData' and 'RH_restData' have a mean of 0. You can use multiple trials just make sure the matrix is oriented properly

% input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
ballData = load('H:\21-12-16_MouseExp\211216_001.txt');
load('H:\21-12-16_MouseExp\211216_001_processed_Layer1_combined.mat');
Thermo = ballData(:,[1 3]); % load raw thermo data
resp = filterThermoDataHz(Thermo,10000); % filter and resample thermo data from 10000 to 20 Hz
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
[Cx,~,~,~,~,fx,confC_RestData,~,cErr_RestData] = coherencyc_eLife2020(resp,[tBrain' brain(:,1)],params);
[Cy,~,~,~,~,fy,confC_RestData,~,cErr_RestData] = coherencyc_eLife2020(resp,[tBrain' brain(:,2)],params);
[PowerX,HzX,ErrorX,PowerY,HzY,ErrorY,PowerT,HzT,ErrorT] = motionSpectrumAnalysisHz(brain,resp(:,2));

figure(1)
subplot(3,1,1)
semilogx(HzX,PowerX,'k')
hold on
f = fill([HzX flip(HzX)],[ErrorX(1,:) flip(ErrorX(2,:))],'r','Linestyle','none');
set(f,'facea',[.2]);
hold off
subplot(3,1,2)
semilogx(HzT,PowerT,'k')
% hold on
% f = fill([HzT flip(HzT)],[ErrorT(1,:) flip(ErrorT(2,:))],'r','Linestyle','none');
% set(f,'facea',[.2]);
% hold off
subplot(3,1,3)
semilogx(fx,Cx);

figure(2)
subplot(3,1,1)
semilogx(HzY,PowerY,'k')
hold on
f = fill([HzY flip(HzY)],[ErrorY(1,:) flip(ErrorY(2,:))],'r','Linestyle','none');
set(f,'facea',[.2]);
hold off
subplot(3,1,2)
semilogx(HzT,PowerT,'k')
% hold on
% f = fill([HzT flip(HzT)],[ErrorT(1,:) flip(ErrorT(2,:))],'r','Linestyle','none');
% set(f,'facea',[.2]);
% hold off
subplot(3,1,3)
semilogx(fy,Cy);


function procData = filterThermoDataHz(Thermo,sampleRate)

close all

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

plot(procData(:,1),procData(:,2))
title('Select upper or lower limit')
limitVal = ginput(1);
close
for n = 1:size(procData,1)
    if procData(n,2) > abs(limitVal(2))
        procData(n,2) = abs(limitVal(2));
    elseif procData(n,2) < -abs(limitVal(2))
        procData(n,2) = -abs(limitVal(2));
    end
end
end

function [PowerX,HzX,ErrorX,PowerY,HzY,ErrorY,PowerT,HzT,ErrorT] = motionSpectrumAnalysisHz(positionData,thermo, plotTF)
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