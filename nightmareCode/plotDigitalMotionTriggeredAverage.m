function plotDigitalMotionTriggeredAverage(tdmsPath)

if exist('tdmsPath','var')
    tdmsDataStruct = openTDMS(tdmsPath);
else
    tdmsDataStruct = openTDMS;
end

Fs = str2double(tdmsDataStruct.CameraFrameratePerSecond);

signalsStruct = tdmsDataStruct.Digital_Data;
signalNames = fieldnames(signalsStruct);

removeFields = {'Puff','Respiration_Sum','Respiration'};
signalNames = setdiff(signalNames, removeFields, 'stable');

nSignals = length(signalNames);

durationSec = 15*60;
samplesToUse = round(durationSec*Fs);

signals = cell(nSignals,1);

for i = 1:nSignals
    data = signalsStruct.(signalNames{i});
    data = data(:);
    signals{i} = data(1:min(samplesToUse,length(data)));
end

%% Build matrix
L = min(cellfun(@length,signals));
X = zeros(L,nSignals);

for i=1:nSignals
    x = signals{i}(1:L);
    x = (x-mean(x))/std(x);
    X(:,i) = x;
end

%% PCA motion component
[coeff,score] = pca(X);
motionSignal = score(:,1);

motionSignal = (motionSignal-mean(motionSignal))/std(motionSignal);

%% Detect motion events
threshold = 2;

eventIdx = find(motionSignal > threshold);

% keep only event onsets
eventIdx = eventIdx([true; diff(eventIdx) > Fs]);

%% Window around event
preSec = 2;
postSec = 2;

preSamples = round(preSec*Fs);
postSamples = round(postSec*Fs);

t = (-preSamples:postSamples)/Fs;

%% Compute averages
avgResponses = zeros(length(t),nSignals);

nEvents = 0;

for k = 1:length(eventIdx)

    idx = eventIdx(k);

    if idx-preSamples < 1 || idx+postSamples > L
        continue
    end
    
    nEvents = nEvents + 1;
    
    for i = 1:nSignals
        seg = X(idx-preSamples:idx+postSamples,i);
        avgResponses(:,i) = avgResponses(:,i) + seg;
    end
    
end

avgResponses = avgResponses / nEvents;

%% Plot

figure
tiledlayout(nSignals,1)

plotNames = strrep(signalNames,'_',' ');

for i=1:nSignals
    
    nexttile
    
    plot(t,avgResponses(:,i),'LineWidth',2)
    xline(0,'k--')
    
    ylabel(plotNames{i})
    
    if i==nSignals
        xlabel('Time around motion event (s)')
    end
    
end

sgtitle(sprintf('Motion Triggered Average (%d events)',nEvents))

end