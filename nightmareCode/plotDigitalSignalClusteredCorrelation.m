function plotDigitalSignalClusteredCorrelation(tdmsPath)

if exist('tdmsPath','var')
    tdmsDataStruct = openTDMS(tdmsPath);
else
    tdmsDataStruct = openTDMS;
end

% Sampling rate
Fs = str2double(tdmsDataStruct.CameraFrameratePerSecond);

% Extract signals
signalsStruct = tdmsDataStruct.Digital_Data;
signalNames = fieldnames(signalsStruct);

% Remove unwanted signals
removeFields = {'Puff','Respiration_Sum'};
signalNames = setdiff(signalNames, removeFields, 'stable');

nSignals = length(signalNames);

% Duration
durationSec = 15*60;
samplesToUse = round(durationSec * Fs);

X = zeros(samplesToUse,nSignals);

for i = 1:nSignals
    
    data = signalsStruct.(signalNames{i});
    data = data(:);
    
    X(:,i) = data(1:min(samplesToUse,length(data)));
    
end

%% Correlation matrix
C = corrcoef(X);

%% Convert to distance for clustering
D = 1 - C;

%% Hierarchical clustering
Z = linkage(squareform(D),'average');

order = optimalleaforder(Z,D);

%% Reorder matrix
Csorted = C(order,order);
sortedNames = signalNames(order);

%% Plot
figure
imagesc(Csorted)
colormap(parula)
colorbar
clim([0 1])

n = length(sortedNames);

for i = 1:n
    for j = 1:n
        text(j,i,sprintf('%.2f',Csorted(i,j)), ...
            'HorizontalAlignment','center','Color','k');
    end
end

sortedNames = strrep(sortedNames,'_',' ');

xticks(1:n)
yticks(1:n)

xticklabels(sortedNames)
yticklabels(sortedNames)

xtickangle(45)

axis square
title('Clustered Signal Correlation')

end