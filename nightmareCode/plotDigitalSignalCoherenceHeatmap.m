function plotDigitalSignalCoherenceHeatmap(method,tdmsPath)

if ~exist('method','var')
    method = 'welch';
end

if exist('tdmsPath','var')
    tdmsDataStruct = openTDMS(tdmsPath);
else
    tdmsDataStruct = openTDMS;
end

% Sampling rate
Fs = str2double(tdmsDataStruct.CameraFrameratePerSecond);

signalsStruct = tdmsDataStruct.Digital_Data;
signalNames = fieldnames(signalsStruct);

removeFields = {'Puff','Respiration_Sum','Respiration'};
signalNames = setdiff(signalNames, removeFields, 'stable');

nSignals = length(signalNames);

durationSec = 15*60;
samplesToUse = round(durationSec * Fs);

signals = cell(nSignals,1);

for i = 1:nSignals
    data = signalsStruct.(signalNames{i});
    data = data(:);
    signals{i} = data(1:min(samplesToUse,length(data)));
end

cohMatrix = zeros(nSignals);

bandLow = 0.1;
bandHigh = 3;

% Chronux params
params.Fs = Fs;
params.tapers = [3 5];
params.pad = 0;
params.fpass = [0 5];

for i = 1:nSignals
    for j = 1:nSignals
        
        x = signals{i};
        y = signals{j};
        
        L = min(length(x),length(y));
        x = x(1:L);
        y = y(1:L);
        
        if strcmp(method,'welch')
            
            [Cxy,f] = mscohere(x,y,[],[],[],Fs);
            
        elseif strcmp(method,'chronux')
            
            [Cxy,~,~,~,f] = coherencyc(x,y,params);
            
        end
        
        bandIdx = f >= bandLow & f <= bandHigh;
        cohMatrix(i,j) = mean(Cxy(bandIdx));
        
    end
end

figure
imagesc(cohMatrix)
colormap(parula)

for i = 1:nSignals
    for j = 1:nSignals
        text(j,i,sprintf('%.2f',cohMatrix(i,j)), ...
            'HorizontalAlignment','center','Color','k');
    end
end

colorbar
clim([0 1])

xticks(1:nSignals)
yticks(1:nSignals)

signalNames = strrep(signalNames,'_',' ');
xticklabels(signalNames)
yticklabels(signalNames)

xtickangle(45)

axis square
title(['Mean Coherence (0.1–3 Hz) - ',method])

end