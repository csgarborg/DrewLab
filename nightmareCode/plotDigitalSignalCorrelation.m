function plotDigitalSignalCorrelation(tdmsPath)

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

signals = cell(nSignals,1);

for i = 1:nSignals
    data = signalsStruct.(signalNames{i});
    data = data(:);
    signals{i} = data(1:min(samplesToUse,length(data)));
end

%% Build matrix for correlation
X = zeros(samplesToUse,nSignals);

for i = 1:nSignals
    X(:,i) = signals{i}(1:samplesToUse);
end

%% Correlation matrix
corrMatrix = corrcoef(X);

figure
imagesc(corrMatrix)
for i = 1:nSignals
    for j = 1:nSignals
        text(j,i,sprintf('%.2f',corrMatrix(i,j)), ...
            'HorizontalAlignment','center', ...
            'Color','k');
    end
end
axis square
colorbar
clim([-1 1])

ax = gca;
ax.XTick = 1:nSignals;
ax.YTick = 1:nSignals;
ax.XTickLabel = signalNames;
ax.YTickLabel = signalNames;
ax.TickLabelInterpreter = 'none';

title('Signal Correlation Matrix')

%% Cross correlation checkerboard

figure
tiledlayout(nSignals,nSignals,'TileSpacing','compact')

for i = 1:nSignals
    for j = 1:nSignals
        
        nexttile
        
        [r,lags] = xcorr(signals{i},signals{j},Fs*5,'coeff'); % ±5 seconds
        lags = lags/Fs;
        
        plot(lags,r,'w')
        xlim([-5 5])
        ylim([0 1])
        
        if i == nSignals
            xlabel('Lag (s)')
        end
        
        if j == 1
            ylabel(signalNames{i},'Interpreter','none')
        end
        
        if i == 1
            title(signalNames{j},'Interpreter','none')
        end
        
    end
end

sgtitle('Pairwise Cross Correlation')

end