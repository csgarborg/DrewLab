function plotDigitalSignalPCA(tdmsPath)

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

signals = zeros(samplesToUse,nSignals);

for i = 1:nSignals
    data = signalsStruct.(signalNames{i});
    data = data(:);
    signals(:,i) = data(1:min(samplesToUse,length(data)));
end

%% PCA
[coeff,score,latent,~,explained] = pca(signals);

%% Variance explained

figure
plot(explained,'o-')
xlabel('Principal Component')
ylabel('Variance Explained (%)')
title('PCA Variance Explained')

%% Component weights

figure
bar(coeff(:,1))
xticks(1:nSignals)
xticklabels(signalNames)
set(gca,'TickLabelInterpreter','none')
ylabel('Weight')
title('PC1 Loadings')

%% PC time series

figure
plot(score(:,1))
xlabel('Time (samples)')
ylabel('PC1 Activity')
title('Principal Component 1 Timecourse')

end