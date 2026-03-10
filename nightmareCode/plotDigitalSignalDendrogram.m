function plotDigitalSignalDendrogram(tdmsPath)

if exist('tdmsPath','var')
    tdmsDataStruct = openTDMS(tdmsPath);
else
    tdmsDataStruct = openTDMS;
end

Fs = str2double(tdmsDataStruct.CameraFrameratePerSecond);

signalsStruct = tdmsDataStruct.Digital_Data;
signalNames = fieldnames(signalsStruct);

removeFields = {'Puff','Respiration_Sum'};
signalNames = setdiff(signalNames, removeFields, 'stable');

nSignals = length(signalNames);

durationSec = 15*60;
samplesToUse = round(durationSec * Fs);

X = zeros(samplesToUse,nSignals);

for i = 1:nSignals
    
    data = signalsStruct.(signalNames{i});
    data = data(:);
    
    X(:,i) = data(1:min(samplesToUse,length(data)));
    
end

% correlation matrix
C = corrcoef(X);

% convert to distance
D = 1 - C;

% hierarchical clustering
Z = linkage(squareform(D),'average');

figure
dendrogram(Z,'Labels',signalNames)

set(gca,'TickLabelInterpreter','none')
title('Signal Similarity Dendrogram')

end