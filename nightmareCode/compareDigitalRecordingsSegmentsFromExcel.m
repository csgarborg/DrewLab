function compareDigitalRecordingsSegmentsFromExcel(excelPath)

%% Load Excel (no headers assumed)

T = readtable(excelPath,'Sheet',2,'ReadVariableNames',false);

nRec = height(T);

fprintf('Loaded %d recordings\n',nRec)

for r = 1:nRec
    
    fprintf('Processing recording %d/%d\n',r,nRec)
    
    tdmsPath = T{r,1}{1};
    
    %% Extract segments
    
    segments = table2array(T(r,2:end));
    segments = segments(~isnan(segments));
    
    if isempty(segments)
        error('No segments provided for row %d', r)
    end
    
    if mod(length(segments),2) ~= 0
        error('Row %d has unmatched start/stop times.', r)
    end
    
    segments = reshape(segments,2,[])';
    
    %% Compute signals from segments
    
    [coeff,score,latent,explained,corrMatrix,signalNames,nSignals] = ...
        computeDigitalSignalsFromSegments_SF(tdmsPath, segments);
    
    %% Labels
    
    [~,name,~] = fileparts(tdmsPath);
    dates{r} = name;
    
    %% Store
    
    results(r).corrMatrix = corrMatrix;
    results(r).coeff = coeff;
    results(r).score = score;
    results(r).latent = latent;
    results(r).explained = explained;
    
end

%% ==============================
%% (UNCHANGED PLOTTING SECTION)
%% ==============================

%% PC1–PC5 variance explained overlay

figure
hold on

for r = 1:nRec
    plot(1:5, results(r).explained(1:5), '-o','LineWidth',2)
end

xlabel('Principal Component')
ylabel('Variance Explained (%)')
title('PCA Variance Explained Across Recordings')

xticks(1:5)
legend(dates,'Location','best')
grid on
ylim([0 100])

%% Correlation matrix comparison

figure
tiledlayout(1,nRec)

for r = 1:nRec
    nexttile
    imagesc(results(r).corrMatrix)
    clim([0 1])
    axis square
    title(dates{r})
end

colorbar
sgtitle('Signal Correlation Comparison')

%% PC1 strength

pc1Var = arrayfun(@(x) x.explained(1), results);

figure
bar(pc1Var)

set(gca,'XTick',1:nRec,'XTickLabel',dates)
ylabel('Variance Explained (%)')
title('PC1 Component Strength')
ylim([0 100])
grid on

%% PC1 contributions

figure
hold on

for r = 1:nRec
    plot(results(r).coeff(:,1), '-o','LineWidth',2)
end

xticks(1:nSignals)
xticklabels(strrep(signalNames,'_',' '))

xlabel('Signal')
ylabel('PC1 Loading')

title('PC1 Signal Contributions Across Recordings')

legend(dates,'Location','best')
grid on
ylim([-1 1])

%% PC2 strength

pc2Var = arrayfun(@(x) x.explained(2), results);

figure
bar(pc2Var)

set(gca,'XTick',1:nRec,'XTickLabel',dates)
ylabel('Variance Explained (%)')
title('PC2 Component Strength')
ylim([0 100])
grid on

%% PC2 contributions

figure
hold on

for r = 1:nRec
    plot(results(r).coeff(:,2), '-o','LineWidth',2)
end

xticks(1:nSignals)
xticklabels(strrep(signalNames,'_',' '))

xlabel('Signal')
ylabel('PC2 Loading')

title('PC2 Signal Contributions Across Recordings')

legend(dates,'Location','best')
grid on
ylim([-1 1])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [coeff,score,latent,explained,corrMatrix,signalNames,nSignals] = ...
    computeDigitalSignalsFromSegments_SF(tdmsPath, segments)

tdmsDataStruct = openTDMS(tdmsPath);

Fs = str2double(tdmsDataStruct.CameraFrameratePerSecond);

signalsStruct = tdmsDataStruct.Digital_Data;
signalNames = fieldnames(signalsStruct);

removeFields = {'Puff','Respiration_Sum','Respiration'};
signalNames = setdiff(signalNames, removeFields, 'stable');

nSignals = length(signalNames);

%% Concatenate segments

concatSignals = cell(nSignals,1);

for i = 1:nSignals
    concatSignals{i} = [];
end

for s = 1:size(segments,1)
    
    startIdx = max(1, round(segments(s,1)*Fs));
    stopIdx  = round(segments(s,2)*Fs);
    
    for i = 1:nSignals
        
        data = signalsStruct.(signalNames{i});
        data = data(:);
        
        stopIdxSafe = min(stopIdx,length(data));
        
        seg = data(startIdx:stopIdxSafe);
        
        concatSignals{i} = [concatSignals{i}; seg];
        
    end
end

%% Build matrix

L = min(cellfun(@length,concatSignals));
X = zeros(L,nSignals);

for i = 1:nSignals
    
    X(:,i) = concatSignals{i}(1:L);
    
end

%% PCA
[coeff,score,latent,~,explained] = pca(X);

%% Correlation
corrMatrix = corrcoef(X);

end