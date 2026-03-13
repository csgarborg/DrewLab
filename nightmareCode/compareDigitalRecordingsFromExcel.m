function compareDigitalRecordingsFromExcel(excelPath)

%% Load recording paths

T = readtable(excelPath,'ReadVariableNames',false);
tdmsPaths = T{:,1};

nRec = length(tdmsPaths);

fprintf('Loaded %d recordings\n',nRec)

%% Run analyses

for r = 1:nRec
    
    fprintf('Processing recording %d/%d\n',r,nRec)
    
    tdmsPath = tdmsPaths{r};
    
    % PCA
    [coeff,score,latent,explained] = plotDigitalSignalPCA(tdmsPath,false);
    
    % Correlation
    [corrMatrix,signals,nSignals,Fs,signalNames] = plotDigitalSignalCorrelation(tdmsPath,false);

    % Recording label
    [~,name,~] = fileparts(tdmsPath);
    dates{r} = name;
    
    results(r).corrMatrix = corrMatrix;
    results(r).coeff = coeff;
    results(r).score = score;
    results(r).latent = latent;
    results(r).explained = explained;
    
end

%% PC1–PC5 variance explained overlay

figure
hold on

for r = 1:nRec
    
    pcVar = results(r).explained(1:5);
    
    plot(1:5, pcVar, '-o','LineWidth',2)
    
end

xlabel('Principal Component')
ylabel('Variance Explained (%)')
title('PCA Variance Explained Across Recordings')

xticks(1:5)

legend(dates,'Location','best')

grid on

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

%% PC1 component strength

pc1Var = zeros(nRec,1);

for r = 1:nRec
    pc1Var(r) = results(r).explained(1);
end

figure

bar(pc1Var)

set(gca,'XTick',1:nRec,'XTickLabel',dates)

ylabel('Variance Explained (%)')
title('PC1 Component Strength')
ylim([0 100])

grid on

%% PC1 signal contribution comparison

figure
hold on

for r = 1:nRec
    
    pc1Load = results(r).coeff(:,1);
    
    plot(1:length(pc1Load), pc1Load, '-o','LineWidth',2)
    
end

xlabel('Signal')
ylabel('PC1 Loading')

xticks(1:nSignals)
xticklabels(signalNames)

title('PC1 Signal Contributions Across Recordings')

legend(dates,'Location','best')

grid on

%% PC2 component strength

pc2Var = zeros(nRec,1);

for r = 1:nRec
    pc2Var(r) = results(r).explained(2);
end

figure

bar(pc2Var)

set(gca,'XTick',1:nRec,'XTickLabel',dates)

ylabel('Variance Explained (%)')
title('PC2 Component Strength')
ylim([0 100])

grid on

%% PC2 signal contribution comparison

figure
hold on

for r = 1:nRec
    
    pc2Load = results(r).coeff(:,2);
    
    plot(1:length(pc2Load), pc2Load, '-o','LineWidth',2)
    
end

xlabel('Signal')
ylabel('PC2 Loading')

xticks(1:nSignals)
xticklabels(signalNames)

title('PC2 Signal Contributions Across Recordings')

legend(dates,'Location','best')

grid on

end