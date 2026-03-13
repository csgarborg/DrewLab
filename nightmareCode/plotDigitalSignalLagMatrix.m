function plotDigitalSignalLagMatrix(tdmsPath)

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
samplesToUse = round(durationSec * Fs);

signals = cell(nSignals,1);

for i = 1:nSignals
    
    data = signalsStruct.(signalNames{i});
    data = data(:);
    
    signals{i} = data(1:min(samplesToUse,length(data)));
    
end

lagMatrix = zeros(nSignals);

maxLagSec = 2;
maxLag = round(maxLagSec*Fs);

for i = 1:nSignals
for j = 1:nSignals
    
    x = signals{i};
    y = signals{j};
    
    L = min(length(x),length(y));
    x = x(1:L);
    y = y(1:L);
    
    [r,lags] = xcorr(x,y,maxLag,'coeff');
    
    [~,idx] = max(r);
    
    lagMatrix(i,j) = lags(idx)/Fs;
    
end
end

%% Plot

figure
imagesc(lagMatrix)

colormap(parula)
colorbar

for i = 1:nSignals
for j = 1:nSignals
    
    text(j,i,sprintf('%.2f',lagMatrix(i,j)),...
        'HorizontalAlignment','center','Color','k');
    
end
end

plotNames = strrep(signalNames,'_',' ');

xticks(1:nSignals)
yticks(1:nSignals)

xticklabels(plotNames)
yticklabels(plotNames)

xtickangle(45)

axis square

title('Signal Lag Matrix (seconds)')

end