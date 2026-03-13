function plotDigitalSignalOverlay(tdmsPath)

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

durationSec = 600;   % shorter window for clarity
samplesToUse = round(durationSec * Fs);

signals = cell(nSignals,1);

for i = 1:nSignals
    
    data = signalsStruct.(signalNames{i});
    data = data(:);
    
    signals{i} = data(1:min(samplesToUse,length(data)));
    
end

%% Build data matrix

L = min(cellfun(@length,signals));
X = zeros(L,nSignals);

for i=1:nSignals
    x = signals{i}(1:L);
    x = (x-mean(x))/std(x);
    X(:,i) = x;
end

%% PCA
[coeff,score] = pca(X);

PC1 = score(:,1);
PC1 = PC1/std(PC1);

%% Time axis
t = (0:L-1)/Fs;

%% Plot

figure
hold on

spacing = 6;

for i = 1:nSignals
    
    plot(t, X(:,i) + spacing*(i-1),'LineWidth',1.5)
    
end

% PC1 overlay at top
plot(t, PC1 + spacing*nSignals,'r','LineWidth',2)

yticks([spacing*(0:nSignals-1) spacing*nSignals])

plotNames = strrep(signalNames,'_',' ');
plotNames{end+1} = 'PC1 Motion Component';

yticklabels(plotNames)

xlabel('Time (seconds)')
ylabel('Signals')
title('Z-Scored Signal Overlay with PCA Motion Component')

end