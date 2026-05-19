function compareREMvsAwake(awakeExcelPath, remExcelPath, mefRemExcelPath, tdmsTF, saveResultsStructTF, roiPath, saveResultsStructPath)

%% Load data
if exist('saveResultsStructPath','var') && ~saveResultsStructTF
    load(saveResultsStructPath);
    awake = resultsStruct.awake;
    rem = resultsStruct.rem;
    mefRem = resultsStruct.mefRem;
    remWake = resultsStruct.remWake;
    mefRemWake = resultsStruct.mefRemWake;
else
    if ~exist("tdmsTF","var") || tdmsTF
        awake  = processExcelSegmentsTDMS_SF(awakeExcelPath);
        rem    = processExcelSegmentsTDMS_SF(remExcelPath);
        mefRem = processExcelSegmentsTDMS_SF(mefRemExcelPath);
    else
        awake  = processExcelSegmentsROISelect_SF(awakeExcelPath,roiPath,false);
        rem    = processExcelSegmentsROISelect_SF(remExcelPath,roiPath,false);
        mefRem = processExcelSegmentsROISelect_SF(mefRemExcelPath,roiPath,false);
        remWake = processExcelSegmentsROISelect_SF(remExcelPath,roiPath,true);
        mefRemWake = processExcelSegmentsROISelect_SF(mefRemExcelPath,roiPath,true);

        % show ROI locations once
        plotFirstFrameWithROIs(remExcelPath, roiPath);
    end
end
if saveResultsStructTF
    resultsStruct.awake = awake;
    resultsStruct.rem = rem;
    resultsStruct.mefRem = mefRem;
    resultsStruct.remWake = remWake;
    resultsStruct.mefRemWake = mefRemWake;
    save(saveResultsStructPath, 'resultsStruct');
end

% awake = processMatFiles_SF(awakeExcelPath);
% rem   = processMatFiles_SF(remExcelPath);

nA = length(awake);
nR = length(rem);
nM = length(mefRem);

signalNames = awake(1).signalNames;
nSignals = awake(1).nSignals;

%% ==============================
%% 1) PCA variance (with CI)
%% ==============================

figure
tiledlayout(1,2)

% --- Raw overlay ---
nexttile
hold on

hA = plot(nan,nan,'-o','Color',[0.3 0.5 1],'DisplayName','Awake');
hR = plot(nan,nan,'--o','Color',[1 0.4 0.4],'DisplayName','REM');
hM = plot(nan,nan,':o','Color',[0.2 0.7 0.2],'DisplayName','MefREM');

for r = 1:nA
    plot(1:nSignals, awake(r).explained(1:nSignals),'-o','Color',[0.3 0.5 1])
end
for r = 1:nR
    plot(1:nSignals, rem(r).explained(1:nSignals),'--o','Color',[1 0.4 0.4])
end
for r = 1:nM
    plot(1:nSignals, mefRem(r).explained(1:nSignals),':o','Color',[0.2 0.7 0.2])
end

legend([hA hR hM])

title('Individual Recordings')
xlabel('PC')
ylabel('% Variance')
grid on
ylim([0 100])
xlim([1 nSignals])

% --- Mean + 95% CI ---
nexttile
hold on

A = reshape([awake.explained],[],nA)';
R = reshape([rem.explained],[],nR)';
M = reshape([mefRem.explained],[],nM)';

meanA = mean(A(:,1:nSignals));
meanR = mean(R(:,1:nSignals));
meanM = mean(M(:,1:nSignals));

ciA = 1.96*std(A(:,1:nSignals))/sqrt(nA);
ciR = 1.96*std(R(:,1:nSignals))/sqrt(nR);
ciM = 1.96*std(M(:,1:nSignals))/sqrt(nM);

x = 1:nSignals;

% Awake CI
fill([x fliplr(x)], ...
     [meanA-ciA fliplr(meanA+ciA)], ...
     [0.3 0.5 1], ...
     'FaceAlpha',0.15,'EdgeColor','none');

plot(x,meanA,'-o','Color',[0.3 0.5 1],'LineWidth',2)

% REM CI
fill([x fliplr(x)], ...
     [meanR-ciR fliplr(meanR+ciR)], ...
     [1 0.4 0.4], ...
     'FaceAlpha',0.15,'EdgeColor','none');

plot(x,meanR,'--o','Color',[1 0.4 0.4],'LineWidth',2)

% MefREM CI
fill([x fliplr(x)], ...
     [meanM-ciM fliplr(meanM+ciM)], ...
     [0.2 0.7 0.2], ...
     'FaceAlpha',0.15,'EdgeColor','none');

plot(x,meanM,':o','Color',[0.2 0.7 0.2],'LineWidth',2)

legend({'Awake CI','Awake','REM CI','REM','MefREM CI','MefREM'})

title('Mean ± 95% CI')
xlabel('PC')
ylabel('% Variance')
grid on
ylim([0 100])
xlim([1 nSignals])

sgtitle('PCA Variance Comparison')

%% ==============================
%% 2) PC1 + PC2 strength (mean ± SD)
%% ==============================

pc1A = arrayfun(@(x)x.explained(1),awake);
pc1R = arrayfun(@(x)x.explained(1),rem);
pc1M = arrayfun(@(x)x.explained(1),mefRem);

pc2A = arrayfun(@(x)x.explained(2),awake);
pc2R = arrayfun(@(x)x.explained(2),rem);
pc2M = arrayfun(@(x)x.explained(2),mefRem);

figure

subplot(1,2,1)
bar([mean(pc1A), mean(pc1R), mean(pc1M)])
hold on
errorbar([1 2 3], ...
    [mean(pc1A) mean(pc1R) mean(pc1M)], ...
    [std(pc1A) std(pc1R) std(pc1M)], ...
    '.w','LineWidth',1.5)
set(gca,'XTickLabel',{'Awake','REM','MefREM'})
title('PC1 Strength')
ylabel('% Variance')
ylim([0 100])

subplot(1,2,2)
bar([mean(pc2A),mean(pc2R),mean(pc2M)])
hold on
errorbar([1 2 3],[mean(pc2A) mean(pc2R) mean(pc2M)], ...
    [std(pc2A) std(pc2R) std(pc2M)],'.w','LineWidth',1.5)
set(gca,'XTickLabel',{'Awake','REM','MefREM'})
title('PC2 Strength')
ylabel('% Variance')
ylim([0 100])

%% ==============================
%% 3) Mean correlation (mean ± SD)
%% ==============================

meanCorrA = arrayfun(@(x)mean(x.corrMatrix(~eye(size(x.corrMatrix)))),awake);
meanCorrR = arrayfun(@(x)mean(x.corrMatrix(~eye(size(x.corrMatrix)))),rem);
meanCorrM = arrayfun(@(x)mean(x.corrMatrix(~eye(size(x.corrMatrix)))),mefRem);

figure
bar([mean(meanCorrA),mean(meanCorrR),mean(meanCorrM)])
hold on
errorbar([1 2 3],[mean(meanCorrA) mean(meanCorrR) mean(meanCorrM)], ...
    [std(meanCorrA) std(meanCorrR) std(meanCorrR)],'.w','LineWidth',1.5)
set(gca,'XTickLabel',{'Awake','REM','MefREM'})
ylabel('Mean Correlation')
title('Average Correlation ± SD')
ylim([0 1])

%% ==============================
%% 4) Correlation matrices (with labels)
%% ==============================

corrA = cat(3,awake.corrMatrix);
corrR = cat(3,rem.corrMatrix);
corrM = cat(3,mefRem.corrMatrix);

meanA = mean(corrA,3);
meanR = mean(corrR,3);
meanM = mean(corrM,3);

figure
tiledlayout(1,3)

nexttile
imagesc(meanA)
title('Awake')
axis square
clim([0 1])
colorbar
xticks(1:nSignals)
yticks(1:nSignals)
xticklabels(strrep(signalNames,'_',' '))
yticklabels(strrep(signalNames,'_',' '))
xtickangle(45)

nexttile
imagesc(meanR)
title('REM')
axis square
clim([0 1])
colorbar
xticks(1:nSignals)
yticks(1:nSignals)
xticklabels(strrep(signalNames,'_',' '))
yticklabels(strrep(signalNames,'_',' '))
xtickangle(45)

nexttile
imagesc(meanM)
title('MefREM')
axis square
clim([0 1])
colorbar
xticks(1:nSignals)
yticks(1:nSignals)
xticklabels(strrep(signalNames,'_',' '))
yticklabels(strrep(signalNames,'_',' '))
xtickangle(45)

sgtitle('Mean Correlation Matrices')

%% ==============================
%% 5) PC1 loadings (with CI)
%% ==============================

figure
tiledlayout(1,2)

% raw
nexttile
hold on

hA = plot(nan,nan,'-o','Color',[0.3 0.5 1],'DisplayName','Awake');
hR = plot(nan,nan,'--o','Color',[1 0.4 0.4],'DisplayName','REM');
hM = plot(nan,nan,':o','Color',[0.2 0.7 0.2],'DisplayName','MefREM');

for r=1:nA, plot(awake(r).coeff(:,1),'-o','Color',[0.3 0.5 1]); end
for r=1:nR, plot(rem(r).coeff(:,1),'--o','Color',[1 0.4 0.4]); end
for r=1:nM, plot(mefRem(r).coeff(:,1),':o','Color',[0.2 0.7 0.2]); end

title('Individual')
xticks(1:nSignals)
xticklabels(strrep(signalNames,'_',' '))
xtickangle(45)

legend([hA hR hM])
ylim([-1 1])

% mean + CI
nexttile
hold on

A = reshape([awake.coeff],nSignals,[],nA);
R = reshape([rem.coeff],nSignals,[],nR);
M = reshape([mefRem.coeff],nSignals,[],nM);

pc1A = squeeze(A(:,1,:))';
pc1R = squeeze(R(:,1,:))';
pc1M = squeeze(M(:,1,:))';

meanA = mean(pc1A);
meanR = mean(pc1R);
meanM = mean(pc1M);

ciA = 1.96*std(pc1A)/sqrt(nA);
ciR = 1.96*std(pc1R)/sqrt(nR);
ciM = 1.96*std(pc1M)/sqrt(nM);

x = 1:nSignals;

% Awake CI
fill([x fliplr(x)], ...
     [meanA-ciA fliplr(meanA+ciA)], ...
     [0.3 0.5 1], ...
     'FaceAlpha',0.15,'EdgeColor','none');

plot(x,meanA,'-o','Color',[0.3 0.5 1],'LineWidth',2)

% REM CI
fill([x fliplr(x)], ...
     [meanR-ciR fliplr(meanR+ciR)], ...
     [1 0.4 0.4], ...
     'FaceAlpha',0.15,'EdgeColor','none');

plot(x,meanR,'--o','Color',[1 0.4 0.4],'LineWidth',2)

% MefREM CI
fill([x fliplr(x)], ...
     [meanM-ciM fliplr(meanM+ciM)], ...
     [0.2 0.7 0.2], ...
     'FaceAlpha',0.15,'EdgeColor','none');

plot(x,meanM,':o','Color',[0.2 0.7 0.2],'LineWidth',2)

legend({'Awake CI','Awake','REM CI','REM','MefREM CI','MefREM'})

title('Mean ± 95% CI')
xticks(1:nSignals)
xticklabels(strrep(signalNames,'_',' '))
xtickangle(45)

sgtitle('PC1 Loadings')
ylim([-1 1])

%% ==============================
%% 6) PC2 loadings (same structure)
%% ==============================

figure
tiledlayout(1,2)

nexttile
hold on

hA = plot(nan,nan,'-o','Color',[0.3 0.5 1],'DisplayName','Awake');
hR = plot(nan,nan,'--o','Color',[1 0.4 0.4],'DisplayName','REM');
hM = plot(nan,nan,':o','Color',[0.2 0.7 0.2],'DisplayName','REM');

for r=1:nA, plot(awake(r).coeff(:,2),'-o','Color',[0.3 0.5 1]); end
for r=1:nR, plot(rem(r).coeff(:,2),'--o','Color',[1 0.4 0.4]); end
for r=1:nM, plot(mefRem(r).coeff(:,2),':o','Color',[0.2 0.7 0.2]); end

title('Individual')
legend([hA hR hM])
ylim([-1 1])

xticks(1:nSignals)
xticklabels(strrep(signalNames,'_',' '))
xtickangle(45)

nexttile
hold on

pc2A = squeeze(A(:,2,:))';
pc2R = squeeze(R(:,2,:))';
pc2M = squeeze(M(:,2,:))';

meanA = mean(pc2A);
meanR = mean(pc2R);
meanM = mean(pc2M);

ciA = 1.96*std(pc2A)/sqrt(nA);
ciR = 1.96*std(pc2R)/sqrt(nR);
ciM = 1.96*std(pc2M)/sqrt(nM);

x = 1:nSignals;

% Awake CI
fill([x fliplr(x)], ...
     [meanA-ciA fliplr(meanA+ciA)], ...
     [0.3 0.5 1], ...
     'FaceAlpha',0.15,'EdgeColor','none');

plot(x,meanA,'-o','Color',[0.3 0.5 1],'LineWidth',2)

% REM CI
fill([x fliplr(x)], ...
     [meanR-ciR fliplr(meanR+ciR)], ...
     [1 0.4 0.4], ...
     'FaceAlpha',0.15,'EdgeColor','none');

plot(x,meanR,'--o','Color',[1 0.4 0.4],'LineWidth',2)

% MefREM CI
fill([x fliplr(x)], ...
     [meanM-ciM fliplr(meanM+ciM)], ...
     [0.2 0.7 0.2], ...
     'FaceAlpha',0.15,'EdgeColor','none');

plot(x,meanM,':o','Color',[0.2 0.7 0.2],'LineWidth',2)

legend({'Awake CI','Awake','REM CI','REM','MefREM CI','MefREM'})

title('Mean ± 95% CI')
xticks(1:nSignals)
xticklabels(strrep(signalNames,'_',' '))
xtickangle(45)
ylim([-1 1])

sgtitle('PC2 Loadings')

%% ==============================
%% 7) PC1 time series comparison
%% ==============================

% figure
% hold on
% 
% for r = 1:nA
%     plot(awake(r).score(:,1),'Color',[0.3 0.5 1])
% end
% 
% for r = 1:nR
%     plot(rem(r).score(:,1),'Color',[1 0.4 0.4])
% end
% 
% title('PC1 Time Series (Overlay)')
% xlabel('Samples')
% ylabel('PC1 Activity')

figure
hold on

% Normalize length
minLen = min([arrayfun(@(x)size(x.score,1),awake), ...
              arrayfun(@(x)size(x.score,1),rem)]);

A = zeros(nA,minLen);
R = zeros(nR,minLen);

for r = 1:nA
    A(r,:) = awake(r).score(1:minLen,1);
end

for r = 1:nR
    R(r,:) = rem(r).score(1:minLen,1);
end

meanA = mean(A);
meanR = mean(R);

ciA = 1.96*std(A)/sqrt(nA);
ciR = 1.96*std(R)/sqrt(nR);

x = 1:minLen;

% Awake
fill([x fliplr(x)], ...
     [meanA-ciA fliplr(meanA+ciA)], ...
     [0.3 0.5 1],'FaceAlpha',0.15,'EdgeColor','none')
plot(x,meanA,'Color',[0.3 0.5 1],'LineWidth',2)

% REM
fill([x fliplr(x)], ...
     [meanR-ciR fliplr(meanR+ciR)], ...
     [1 0.4 0.4],'FaceAlpha',0.15,'EdgeColor','none')
plot(x,meanR,'Color',[1 0.4 0.4],'LineWidth',2)

legend({'Awake CI','Awake','REM CI','REM'})
title('PC1 Time Series (Mean ± 95% CI)')
xlabel('Samples')
ylabel('PC1 Activity')

%% ==============================
%% 8) Global vs Local motion index
%% ==============================

% simple definition:
% global = PC1 variance
% local  = mean of PC2–PC5

globalA = pc1A;
localA  = mean(A(:,2:5,:),2);

globalR = pc1R;
localR  = mean(R(:,2:5,:),2);

globalM = pc1M;
localM  = mean(M(:,2:5,:),2);

GL_A = zeros(nA,1);
GL_R = zeros(nR,1);
GL_M = zeros(nM,1);

for r = 1:nA
    e = awake(r).explained;
    GL_A(r) = e(1) / mean(e(2:min(5,end)));
end

for r = 1:nR
    e = rem(r).explained;
    GL_R(r) = e(1) / mean(e(2:min(5,end)));
end

for r = 1:nM
    e = mefRem(r).explained;
    GL_M(r) = e(1) / mean(e(2:min(5,end)));
end

figure
bar([mean(GL_A), mean(GL_R), mean(GL_M)])
set(gca,'XTickLabel',{'Awake','REM','MefREM'})
ylabel('Global / Local Ratio')
title('Global vs Local Motion Index')

%% ==============================
%% 9) Frequency analysis (ROI + PC1)
%% ==============================

figure
tiledlayout(2,2)

% --- PC1 PSD ---
nexttile
hold on

for r = 1:nA
    x = awake(r).score(:,1);
    [Pxx,F] = pwelch(x,[],[],[],1); % normalized Fs if unknown
    plot(F,Pxx,'Color',[0.3 0.5 1])
end

for r = 1:nR
    x = rem(r).score(:,1);
    [Pxx,F] = pwelch(x,[],[],[],1);
    plot(F,Pxx,'Color',[1 0.4 0.4])
end

for r = 1:nM
    x = mefRem(r).score(:,1);
    [Pxx,F] = pwelch(x,[],[],[],1);
    plot(F,Pxx,'Color',[0.2 0.7 0.2])
end

title('PC1 Power Spectrum')
xlabel('Frequency')
ylabel('Power')
grid on

% --- Dominant frequency ---
nexttile
domA = zeros(nA,1);
domR = zeros(nR,1);
domM = zeros(nM,1);

for r = 1:nA
    [Pxx,F] = pwelch(awake(r).score(:,1),[],[],[],1);
    [~,idx] = max(Pxx);
    domA(r) = F(idx);
end

for r = 1:nR
    [Pxx,F] = pwelch(rem(r).score(:,1),[],[],[],1);
    [~,idx] = max(Pxx);
    domR(r) = F(idx);
end

for r = 1:nM
    [Pxx,F] = pwelch(mefRem(r).score(:,1),[],[],[],1);
    [~,idx] = max(Pxx);
    domM(r) = F(idx);
end

bar([mean(domA), mean(domR), mean(domM)])
hold on
errorbar([1 2 3], ...
    [mean(domA), mean(domR), mean(domM)], ...
    [std(domA), std(domR), std(domM)], ...
    '.w','LineWidth',1.5)

set(gca,'XTickLabel',{'Awake','REM','MefREM'})
ylabel('Dominant Frequency')
title('PC1 Dominant Frequency')

% --- Spectral entropy ---
nexttile
entA = zeros(nA,1);
entR = zeros(nR,1);
entM = zeros(nM,1);

for r = 1:nA
    [Pxx,~] = pwelch(awake(r).score(:,1),[],[],[],1);
    p = Pxx ./ sum(Pxx);
    entA(r) = -sum(p .* log2(p + eps));
end

for r = 1:nR
    [Pxx,~] = pwelch(rem(r).score(:,1),[],[],[],1);
    p = Pxx ./ sum(Pxx);
    entR(r) = -sum(p .* log2(p + eps));
end

for r = 1:nM
    [Pxx,~] = pwelch(mefRem(r).score(:,1),[],[],[],1);
    p = Pxx ./ sum(Pxx);
    entM(r) = -sum(p .* log2(p + eps));
end

bar([mean(entA), mean(entR), mean(entM)])
hold on
errorbar([1 2 3], ...
    [mean(entA), mean(entR), mean(entM)], ...
    [std(entA), std(entR), std(entM)], ...
    '.w','LineWidth',1.5)

set(gca,'XTickLabel',{'Awake','REM','MefREM'})
ylabel('Spectral Entropy')
title('PC1 Spectral Entropy')

% --- High / Low frequency ratio ---
nexttile
hlA = zeros(nA,1);
hlR = zeros(nR,1);
hlM = zeros(nM,1);

for r = 1:nA
    [Pxx,F] = pwelch(awake(r).score(:,1),[],[],[],1);
    low = sum(Pxx(F <= 0.1));
    high = sum(Pxx(F > 0.1));
    hlA(r) = high / (low + eps);
end

for r = 1:nR
    [Pxx,F] = pwelch(rem(r).score(:,1),[],[],[],1);
    low = sum(Pxx(F <= 0.1));
    high = sum(Pxx(F > 0.1));
    hlR(r) = high / (low + eps);
end

for r = 1:nM
    [Pxx,F] = pwelch(mefRem(r).score(:,1),[],[],[],1);
    low = sum(Pxx(F <= 0.1));
    high = sum(Pxx(F > 0.1));
    hlM(r) = high / (low + eps);
end

bar([mean(hlA), mean(hlR), mean(hlM)])
hold on
errorbar([1 2 3], ...
    [mean(hlA), mean(hlR), mean(hlM)], ...
    [std(hlA), std(hlR), std(hlM)], ...
    '.w','LineWidth',1.5)

set(gca,'XTickLabel',{'Awake','REM','MefREM'})
ylabel('High / Low Ratio')
title('High vs Low Frequency Power')

sgtitle('Frequency Analysis')

%% ==============================
%% 10) Participation ratio + effective dimensionality
%% ==============================

PR_A = zeros(nA,1);
PR_R = zeros(nR,1);
PR_M = zeros(nM,1);

Dim90_A = zeros(nA,1);
Dim90_R = zeros(nR,1);
Dim90_M = zeros(nM,1);

for r = 1:nA
    e = awake(r).explained(:);
    PR_A(r) = (sum(e)^2) / sum(e.^2);

    c = cumsum(e) / sum(e);
    Dim90_A(r) = find(c >= 0.90,1,'first');
end

for r = 1:nR
    e = rem(r).explained(:);
    PR_R(r) = (sum(e)^2) / sum(e.^2);

    c = cumsum(e) / sum(e);
    Dim90_R(r) = find(c >= 0.90,1,'first');
end

for r = 1:nM
    e = mefRem(r).explained(:);
    PR_M(r) = (sum(e)^2) / sum(e.^2);

    c = cumsum(e) / sum(e);
    Dim90_M(r) = find(c >= 0.90,1,'first');
end

figure

subplot(1,2,1)
bar([mean(PR_A), mean(PR_R), mean(PR_M)])
hold on
errorbar([1 2 3], ...
    [mean(PR_A), mean(PR_R), mean(PR_M)], ...
    [std(PR_A), std(PR_R), std(PR_M)], ...
    '.w','LineWidth',1.5)

set(gca,'XTickLabel',{'Awake','REM','MefREM'})
ylabel('Participation Ratio')
title('Effective Dimensionality')

subplot(1,2,2)
bar([mean(Dim90_A), mean(Dim90_R), mean(Dim90_M)])
hold on
errorbar([1 2 3], ...
    [mean(Dim90_A), mean(Dim90_R), mean(Dim90_M)], ...
    [std(Dim90_A), std(Dim90_R), std(Dim90_M)], ...
    '.w','LineWidth',1.5)

set(gca,'XTickLabel',{'Awake','REM','MefREM'})
ylabel('# PCs for 90% Variance')
title('Dimensionality to 90% Variance')

%% ==============================
%% 11) Variability / Fano factor
%% ==============================

F_A = zeros(nA,1);
F_R = zeros(nR,1);
F_M = zeros(nM,1);

for r = 1:nA
    x = abs(awake(r).score(:,1));
    F_A(r) = var(x) / (mean(x) + eps);
end

for r = 1:nR
    x = abs(rem(r).score(:,1));
    F_R(r) = var(x) / (mean(x) + eps);
end

for r = 1:nM
    x = abs(mefRem(r).score(:,1));
    F_M(r) = var(x) / (mean(x) + eps);
end

figure
bar([mean(F_A), mean(F_R), mean(F_M)])
hold on
errorbar([1 2 3], ...
    [mean(F_A), mean(F_R), mean(F_M)], ...
    [std(F_A), std(F_R), std(F_M)], ...
    '.w','LineWidth',1.5)

set(gca,'XTickLabel',{'Awake','REM','MefREM'})
ylabel('Fano Factor')
title('Burstiness / Variability')

%% ==============================
%% 12) Simple twitch event detection
%% ==============================

thrMult = 2; % threshold = mean + 2*std

rateA = zeros(nA,1);
rateR = zeros(nR,1);
rateM = zeros(nM,1);

for r = 1:nA
    x = abs(awake(r).score(:,1));
    thr = mean(x) + thrMult*std(x);
    events = x > thr;
    rateA(r) = sum(diff([0; events]) == 1);
end

for r = 1:nR
    x = abs(rem(r).score(:,1));
    thr = mean(x) + thrMult*std(x);
    events = x > thr;
    rateR(r) = sum(diff([0; events]) == 1);
end

for r = 1:nM
    x = abs(mefRem(r).score(:,1));
    thr = mean(x) + thrMult*std(x);
    events = x > thr;
    rateM(r) = sum(diff([0; events]) == 1);
end

figure
bar([mean(rateA), mean(rateR), mean(rateM)])
hold on
errorbar([1 2 3], ...
    [mean(rateA), mean(rateR), mean(rateM)], ...
    [std(rateA), std(rateR), std(rateM)], ...
    '.w','LineWidth',1.5)

set(gca,'XTickLabel',{'Awake','REM','MefREM'})
ylabel('Event Count')
title('Twitch Event Rate')

%% ==============================
%% 13) REM vs MefREM Event Lengths
%% Every segment as its own point
%% Using plotSpread + mean ± SD
%% ==============================

% Assumes:
% rem(r).segmentLengths
% mefRem(r).segmentLengths
%
% Each individual segment length will be plotted
% as its own point (NO per-recording averaging)

allLenR = [];
allLenM = [];

for r = 1:nR
    if isfield(rem(r),'segmentLengths') && ~isempty(rem(r).segmentLengths)
        allLenR = [allLenR; rem(r).segmentLengths(:)];
    end
end

for r = 1:nM
    if isfield(mefRem(r),'segmentLengths') && ~isempty(mefRem(r).segmentLengths)
        allLenM = [allLenM; mefRem(r).segmentLengths(:)];
    end
end

% Remove NaNs just in case
allLenR = allLenR(~isnan(allLenR));
allLenM = allLenM(~isnan(allLenM));

figure
hold on

plotSpread({allLenR, allLenM}, ...
    'xNames', {'REM','MefREM'}, ...
    'showMM', 5, ... % mean ± SD
    'distributionColors', {[1 0.4 0.4], [0.2 0.7 0.2]}, ...
    'distributionMarkers', {'o','o'}, ...
    'spreadWidth', 0.6)

ylabel('Event Length (s)')
title('REM vs MefREM Event Lengths')
grid on

%% ==============================
%% Percent Time Moving by ROI
%% Using plotSpread
%%
%% Assumes:
%% awake(r).percentAbove3
%% awake(r).percentAbove2
%% awake(r).percentAbove1point5
%%
%% These can now be MATRICES:
%% rows = multiple values/events for one ROI
%% cols = ROIs (16 columns)
%%
%% Example:
%% size(percentAbove3) = [N x 16]
%%
%% For plotting:
%% all values for a given ROI are pooled together
%% across rows AND across recordings
%%
%% Creates:
%% Figure 1 -> z > 3
%% Figure 2 -> z > 2
%% Figure 3 -> z > 1.5
%%
%% Each figure:
%% 16 ROI subplots
%% each subplot compares:
%% Awake vs REM vs MefREM
%% ==============================

roiLabels = signalNames;
nROI = length(roiLabels);

cA = [0.3 0.5 1];
cR = [1 0.4 0.4];
cM = [0.2 0.7 0.2];

%% ==============================
%% Figure 1 : z > 3
%% ==============================

figure
tiledlayout(4,4,'TileSpacing','compact')

for i = 1:nROI
    
    nexttile
    hold on
    
    valsA = [];
    valsR = [];
    valsM = [];
    
    % pool across all awake recordings
    for r = 1:nA
        if ~isempty(awake(r).percentAbove3)
            valsA = [valsA; awake(r).percentAbove3(:,i)];
        end
    end
    
    % pool across all REM recordings
    for r = 1:nR
        if ~isempty(rem(r).percentAbove3)
            valsR = [valsR; rem(r).percentAbove3(:,i)];
        end
    end
    
    % pool across all MefREM recordings
    for r = 1:nM
        if ~isempty(mefRem(r).percentAbove3)
            valsM = [valsM; mefRem(r).percentAbove3(:,i)];
        end
    end
    
    valsA = valsA(~isnan(valsA));
    valsR = valsR(~isnan(valsR));
    valsM = valsM(~isnan(valsM));
    
    plotSpread({valsA, valsR, valsM}, ...
        'xNames', {'Awake','REM','MefREM'}, ...
        'showMM', 5, ...
        'distributionColors', {cA, cR, cM}, ...
        'distributionMarkers', {'o','o','o'}, ...
        'spreadWidth', 0.5)
    
    title(strrep(roiLabels{i},'_',' '))
    grid on
    
    if i == 1
        ylabel('% Time Moving')
    end
end

sgtitle('Percent Time Moving by ROI (z > 3)')


%% ==============================
%% Figure 2 : z > 2
%% ==============================

figure
tiledlayout(4,4,'TileSpacing','compact')

for i = 1:nROI
    
    nexttile
    hold on
    
    valsA = [];
    valsR = [];
    valsM = [];
    
    for r = 1:nA
        if ~isempty(awake(r).percentAbove2)
            valsA = [valsA; awake(r).percentAbove2(:,i)];
        end
    end
    
    for r = 1:nR
        if ~isempty(rem(r).percentAbove2)
            valsR = [valsR; rem(r).percentAbove2(:,i)];
        end
    end
    
    for r = 1:nM
        if ~isempty(mefRem(r).percentAbove2)
            valsM = [valsM; mefRem(r).percentAbove2(:,i)];
        end
    end
    
    valsA = valsA(~isnan(valsA));
    valsR = valsR(~isnan(valsR));
    valsM = valsM(~isnan(valsM));
    
    plotSpread({valsA, valsR, valsM}, ...
        'xNames', {'Awake','REM','MefREM'}, ...
        'showMM', 5, ...
        'distributionColors', {cA, cR, cM}, ...
        'distributionMarkers', {'o','o','o'}, ...
        'spreadWidth', 0.5)
    
    title(strrep(roiLabels{i},'_',' '))
    grid on
    
    if i == 1
        ylabel('% Time Moving')
    end
end

sgtitle('Percent Time Moving by ROI (z > 2)')


%% ==============================
%% Figure 3 : z > 1.5
%% ==============================

figure
tiledlayout(4,4,'TileSpacing','compact')

for i = 1:nROI
    
    nexttile
    hold on
    
    valsA = [];
    valsR = [];
    valsM = [];
    
    for r = 1:nA
        if ~isempty(awake(r).percentAbove1point5)
            valsA = [valsA; awake(r).percentAbove1point5(:,i)];
        end
    end
    
    for r = 1:nR
        if ~isempty(rem(r).percentAbove1point5)
            valsR = [valsR; rem(r).percentAbove1point5(:,i)];
        end
    end
    
    for r = 1:nM
        if ~isempty(mefRem(r).percentAbove1point5)
            valsM = [valsM; mefRem(r).percentAbove1point5(:,i)];
        end
    end
    
    valsA = valsA(~isnan(valsA));
    valsR = valsR(~isnan(valsR));
    valsM = valsM(~isnan(valsM));
    
    plotSpread({valsA, valsR, valsM}, ...
        'xNames', {'Awake','REM','MefREM'}, ...
        'showMM', 5, ...
        'distributionColors', {cA, cR, cM}, ...
        'distributionMarkers', {'o','o','o'}, ...
        'spreadWidth', 0.5)
    
    title(strrep(roiLabels{i},'_',' '))
    grid on
    
    if i == 1
        ylabel('% Time Moving')
    end
end

sgtitle('Percent Time Moving by ROI (z > 1.5)')


%% ==============================
%% Mean ROI Motion Across Event
%%
%% Assumes:
%% awake(r).meanROIMotion
%% rem(r).meanROIMotion
%% mefRem(r).meanROIMotion
%%
%% Each is a cell array:
%% results(r).meanROIMotion{k}
%%
%% where each cell contains:
%% column 1 = percent of event elapsed (0 → 100)
%% column 2 = mean ROI z-score at that time
%%
%% This plots every event trace separately
%% with semi-transparent lines
%% ==============================

figure
tiledlayout(1,3,'TileSpacing','compact')

%% =====================================
%% AWAKE
%% =====================================

nexttile
hold on

for r = 1:nA
    
    if isfield(awake(r),'meanROIMotion') && ~isempty(awake(r).meanROIMotion)
        
        for k = 1:length(awake(r).meanROIMotion)
            
            thisData = awake(r).meanROIMotion{k};
            
            if isempty(thisData) || size(thisData,2) < 2
                continue
            end
            
            x = thisData(:,1);   % percent through event
            y = thisData(:,2);   % mean ROI z-score
            
            plot(x, y, ...
                'Color', [0.3 0.5 1 0.1], ...
                'LineWidth', 1.5)
        end
    end
end

title('Awake')
xlabel('% Event Elapsed')
ylabel('Mean ROI Z-score')
grid on
xlim([0 100])
ylim([-2 12])



%% =====================================
%% REM
%% =====================================

nexttile
hold on

for r = 1:nR
    
    if isfield(rem(r),'meanROIMotion') && ~isempty(rem(r).meanROIMotion)
        
        for k = 1:length(rem(r).meanROIMotion)
            
            thisData = rem(r).meanROIMotion{k};
            
            if isempty(thisData) || size(thisData,2) < 2
                continue
            end
            
            x = thisData(:,1);
            y = thisData(:,2);
            
            plot(x, y, ...
                'Color', [1 0.4 0.4 0.1], ...
                'LineWidth', 1.5)
        end
    end
end

title('REM')
xlabel('% Event Elapsed')
ylabel('Mean ROI Z-score')
grid on
xlim([0 100])
ylim([-2 12])


%% =====================================
%% MefREM
%% =====================================

nexttile
hold on

for r = 1:nM
    
    if isfield(mefRem(r),'meanROIMotion') && ~isempty(mefRem(r).meanROIMotion)
        
        for k = 1:length(mefRem(r).meanROIMotion)
            
            thisData = mefRem(r).meanROIMotion{k};
            
            if isempty(thisData) || size(thisData,2) < 2
                continue
            end
            
            x = thisData(:,1);
            y = thisData(:,2);
            
            plot(x, y, ...
                'Color', [0.2 0.7 0.2 0.1], ...
                'LineWidth', 1.5)
        end
    end
end

title('MefREM')
xlabel('% Event Elapsed')
ylabel('Mean ROI Z-score')
grid on
xlim([0 100])
ylim([-2 12])


sgtitle('Mean ROI Motion Across Events')

%% ==============================
%% Mean ROI Motion Across Event
%% Binned at 1%, 5%, and 10%
%%
%% Assumes:
%% awake(r).meanROIMotion{k}
%% rem(r).meanROIMotion{k}
%% mefRem(r).meanROIMotion{k}
%%
%% Each cell:
%% col 1 = percent through event (0–100)
%% col 2 = mean ROI z-score
%%
%% Creates:
%% Figure 1 -> 1% bins
%% Figure 2 -> 5% bins
%% Figure 3 -> 10% bins
%% ==============================

binSizes = [1 5 10];

for b = 1:length(binSizes)

    binSize = binSizes(b);

    figure
    tiledlayout(1,3,'TileSpacing','compact')

    %% =====================================
    %% AWAKE
    %% =====================================

    nexttile
    hold on

    for r = 1:nA

        if isfield(awake(r),'meanROIMotion') && ~isempty(awake(r).meanROIMotion)

            for k = 1:length(awake(r).meanROIMotion)

                thisData = awake(r).meanROIMotion{k};

                if isempty(thisData) || size(thisData,2) < 2
                    continue
                end

                x = thisData(:,1);
                y = thisData(:,2);

                % ---- binning ----
                binEdges = 0:binSize:100;
                binCenters = binEdges(1:end-1) + binSize/2;
                yBinned = nan(size(binCenters));

                for i = 1:length(binCenters)
                    idx = x >= binEdges(i) & x < binEdges(i+1);
                    if any(idx)
                        yBinned(i) = mean(y(idx),'omitnan');
                    end
                end

                plot(binCenters, yBinned, ...
                    'Color', [0.3 0.5 1 0.2], ...
                    'LineWidth', 1.5)
            end
        end
    end

    title('Awake')
    xlabel('% Event Elapsed')
    ylabel('Mean ROI Z-score')
    grid on
    xlim([0 100])
    ylim([-1 2])


    %% =====================================
    %% REM
    %% =====================================

    nexttile
    hold on

    for r = 1:nR

        if isfield(rem(r),'meanROIMotion') && ~isempty(rem(r).meanROIMotion)

            for k = 1:length(rem(r).meanROIMotion)

                thisData = rem(r).meanROIMotion{k};

                if isempty(thisData) || size(thisData,2) < 2
                    continue
                end

                x = thisData(:,1);
                y = thisData(:,2);

                % ---- binning ----
                binEdges = 0:binSize:100;
                binCenters = binEdges(1:end-1) + binSize/2;
                yBinned = nan(size(binCenters));

                for i = 1:length(binCenters)
                    idx = x >= binEdges(i) & x < binEdges(i+1);
                    if any(idx)
                        yBinned(i) = mean(y(idx),'omitnan');
                    end
                end

                plot(binCenters, yBinned, ...
                    'Color', [1 0.4 0.4 0.2], ...
                    'LineWidth', 1.5)
            end
        end
    end

    title('REM')
    xlabel('% Event Elapsed')
    ylabel('Mean ROI Z-score')
    grid on
    xlim([0 100])
    ylim([-1 2])


    %% =====================================
    %% MefREM
    %% =====================================

    nexttile
    hold on

    for r = 1:nM

        if isfield(mefRem(r),'meanROIMotion') && ~isempty(mefRem(r).meanROIMotion)

            for k = 1:length(mefRem(r).meanROIMotion)

                thisData = mefRem(r).meanROIMotion{k};

                if isempty(thisData) || size(thisData,2) < 2
                    continue
                end

                x = thisData(:,1);
                y = thisData(:,2);

                % ---- binning ----
                binEdges = 0:binSize:100;
                binCenters = binEdges(1:end-1) + binSize/2;
                yBinned = nan(size(binCenters));

                for i = 1:length(binCenters)
                    idx = x >= binEdges(i) & x < binEdges(i+1);
                    if any(idx)
                        yBinned(i) = mean(y(idx),'omitnan');
                    end
                end

                plot(binCenters, yBinned, ...
                    'Color', [0.2 0.7 0.2 0.2], ...
                    'LineWidth', 1.5)
            end
        end
    end

    title('MefREM')
    xlabel('% Event Elapsed')
    ylabel('Mean ROI Z-score')
    grid on
    xlim([0 100])
    ylim([-1 2])

    sgtitle(sprintf('Mean ROI Motion Across Events (%d%% Bins)', binSize))

end

    %% ==============================
%% EMG Power Across Event
%%
%% Same format as meanROIMotion plots
%%
%% Assumes:
%% awake(r).emgPower{k}
%% rem(r).emgPower{k}
%% mefRem(r).emgPower{k}
%%
%% Each cell contains:
%% column 1 = percent of event elapsed (0 → 100)
%% column 2 = EMG power at that time
%%
%% Creates 4 figures:
%% Figure 1 -> raw traces
%% Figure 2 -> 1% bins
%% Figure 3 -> 5% bins
%% Figure 4 -> 10% bins
%% ==============================

%% =========================================
%% FIGURE 1 : RAW TRACES
%% =========================================

figure
tiledlayout(1,3,'TileSpacing','compact')

%% ---------- AWAKE ----------
nexttile
hold on

for r = 1:nA
    if isfield(awake(r),'emgPower') && ~isempty(awake(r).emgPower)

        for k = 1:length(awake(r).emgPower)

            thisData = awake(r).emgPower{k};

            if isempty(thisData) || size(thisData,2) < 2
                continue
            end

            x = thisData(:,1);
            y = thisData(:,2);

            plot(x,y, ...
                'Color',[0.3 0.5 1 0.25], ...
                'LineWidth',1.5)
        end
    end
end

title('Awake')
xlabel('% Event Elapsed')
ylabel('EMG Power')
grid on
xlim([0 100])
ylim([-6 0])


%% ---------- REM ----------
nexttile
hold on

for r = 1:nR
    if isfield(rem(r),'emgPower') && ~isempty(rem(r).emgPower)

        for k = 1:length(rem(r).emgPower)

            thisData = rem(r).emgPower{k};

            if isempty(thisData) || size(thisData,2) < 2
                continue
            end

            x = thisData(:,1);
            y = thisData(:,2);

            plot(x,y, ...
                'Color',[1 0.4 0.4 0.25], ...
                'LineWidth',1.5)
        end
    end
end

title('REM')
xlabel('% Event Elapsed')
ylabel('EMG Power')
grid on
xlim([0 100])
ylim([-6 0])


%% ---------- MefREM ----------
nexttile
hold on

for r = 1:nM
    if isfield(mefRem(r),'emgPower') && ~isempty(mefRem(r).emgPower)

        for k = 1:length(mefRem(r).emgPower)

            thisData = mefRem(r).emgPower{k};

            if isempty(thisData) || size(thisData,2) < 2
                continue
            end

            x = thisData(:,1);
            y = thisData(:,2);

            plot(x,y, ...
                'Color',[0.2 0.7 0.2 0.25], ...
                'LineWidth',1.5)
        end
    end
end

title('MefREM')
xlabel('% Event Elapsed')
ylabel('EMG Power')
grid on
xlim([0 100])
ylim([-6 0])

sgtitle('EMG Power Across Events (Raw)')


%% =========================================
%% FIGURES 2–4 : BINNED (1%, 5%, 10%)
%% =========================================

binSizes = [1 5 10];

for b = 1:length(binSizes)

    binSize = binSizes(b);

    figure
    tiledlayout(1,3,'TileSpacing','compact')

    %% ---------- AWAKE ----------
    nexttile
    hold on

    for r = 1:nA
        if isfield(awake(r),'emgPower') && ~isempty(awake(r).emgPower)

            for k = 1:length(awake(r).emgPower)

                thisData = awake(r).emgPower{k};

                if isempty(thisData) || size(thisData,2) < 2
                    continue
                end

                x = thisData(:,1);
                y = thisData(:,2);

                binEdges = 0:binSize:100;
                binCenters = binEdges(1:end-1) + binSize/2;
                yBinned = nan(size(binCenters));

                for i = 1:length(binCenters)
                    idx = x >= binEdges(i) & x < binEdges(i+1);
                    if any(idx)
                        yBinned(i) = mean(y(idx),'omitnan');
                    end
                end

                plot(binCenters,yBinned, ...
                    'Color',[0.3 0.5 1 0.25], ...
                    'LineWidth',1.5)
            end
        end
    end

    title('Awake')
    xlabel('% Event Elapsed')
    ylabel('EMG Power')
    grid on
    xlim([0 100])
    ylim([-6 0])


    %% ---------- REM ----------
    nexttile
    hold on

    for r = 1:nR
        if isfield(rem(r),'emgPower') && ~isempty(rem(r).emgPower)

            for k = 1:length(rem(r).emgPower)

                thisData = rem(r).emgPower{k};

                if isempty(thisData) || size(thisData,2) < 2
                    continue
                end

                x = thisData(:,1);
                y = thisData(:,2);

                binEdges = 0:binSize:100;
                binCenters = binEdges(1:end-1) + binSize/2;
                yBinned = nan(size(binCenters));

                for i = 1:length(binCenters)
                    idx = x >= binEdges(i) & x < binEdges(i+1);
                    if any(idx)
                        yBinned(i) = mean(y(idx),'omitnan');
                    end
                end

                plot(binCenters,yBinned, ...
                    'Color',[1 0.4 0.4 0.25], ...
                    'LineWidth',1.5)
            end
        end
    end

    title('REM')
    xlabel('% Event Elapsed')
    ylabel('EMG Power')
    grid on
    xlim([0 100])
    ylim([-6 0])


    %% ---------- MefREM ----------
    nexttile
    hold on

    for r = 1:nM
        if isfield(mefRem(r),'emgPower') && ~isempty(mefRem(r).emgPower)

            for k = 1:length(mefRem(r).emgPower)

                thisData = mefRem(r).emgPower{k};

                if isempty(thisData) || size(thisData,2) < 2
                    continue
                end

                x = thisData(:,1);
                y = thisData(:,2);

                binEdges = 0:binSize:100;
                binCenters = binEdges(1:end-1) + binSize/2;
                yBinned = nan(size(binCenters));

                for i = 1:length(binCenters)
                    idx = x >= binEdges(i) & x < binEdges(i+1);
                    if any(idx)
                        yBinned(i) = mean(y(idx),'omitnan');
                    end
                end

                plot(binCenters,yBinned, ...
                    'Color',[0.2 0.7 0.2 0.25], ...
                    'LineWidth',1.5)
            end
        end
    end

    title('MefREM')
    xlabel('% Event Elapsed')
    ylabel('EMG Power')
    grid on
    xlim([0 100])
    ylim([-6 0])

    sgtitle(sprintf('EMG Power Across Events (%d%% Bins)', binSize))

end

%% ==============================
%% EMG Power START / END Across Events
%%
%% Layout:
%% [ Awake      REM      MefREM
%%   AwakeMean  REMMean  MefREMMean ]
%%
%% Top row    = individual event traces
%% Bottom row = mean ± 95% CI
%%
%% X-axis = Time (seconds)
%% ==============================

varsToPlot   = {'emgPowerStart','emgPowerEnd'};
titlesToPlot = {'EMG Power Start','EMG Power End'};

cA = [0.3 0.5 1];
cR = [1 0.4 0.4];
cM = [0.2 0.7 0.2];

alphaRaw = 0.25;
alphaCI  = 0.18;

for v = 1:length(varsToPlot)

    varName   = varsToPlot{v};
    mainTitle = titlesToPlot{v};

    figure
    tiledlayout(2,3,'TileSpacing','compact')

    groupStructs = {awake, rem, mefRem};
    groupNs      = [nA, nR, nM];
    groupTitles  = {'Awake','REM','MefREM'};
    groupColors  = {cA, cR, cM};

    for g = 1:3

        dataStruct = groupStructs{g};
        nGroup     = groupNs(g);
        thisTitle  = groupTitles{g};
        thisColor  = groupColors{g};

        allX = {};
        allY = {};

        %% =====================================
        %% TOP ROW: Individual traces
        %% =====================================

        nexttile(g)
        hold on

        for r = 1:nGroup

            if isfield(dataStruct(r),varName) && ~isempty(dataStruct(r).(varName))

                for k = 1:length(dataStruct(r).(varName))

                    thisData = dataStruct(r).(varName){k};

                    if isempty(thisData) || size(thisData,2) < 2
                        continue
                    end

                    x = thisData(:,1);
                    y = thisData(:,2);

                    % force column vectors
                    x = x(:);
                    y = y(:);

                    allX{end+1} = x; %#ok<AGROW>
                    allY{end+1} = y; %#ok<AGROW>

                    plot(x, y, ...
                        'Color', [thisColor alphaRaw], ...
                        'LineWidth', 1.5)
                end
            end
        end

        title(thisTitle)
        xlabel('Time (s)')
        ylabel(mainTitle)
        grid on
        ylim([-6 0])


        %% =====================================
        %% BOTTOM ROW: Mean ± 95% CI
        %% =====================================

        nexttile(g+3)
        hold on

        if ~isempty(allX)

            % ---------------------------------
            % IMPORTANT FIX:
            % interpolate onto common time axis
            % instead of truncating by minLen
            % ---------------------------------

            minStart = max(cellfun(@(x) min(x), allX));
            maxEnd   = min(cellfun(@(x) max(x), allX));

            % fallback protection
            if maxEnd > minStart

                nInterp = 300;
                xCommon = linspace(minStart, maxEnd, nInterp);

                Yinterp = nan(length(allY), nInterp);

                for i = 1:length(allY)

                    x = allX{i};
                    y = allY{i};

                    % remove duplicate x values if present
                    [xUnique, ia] = unique(x);
                    yUnique = y(ia);

                    if length(xUnique) < 2
                        continue
                    end

                    Yinterp(i,:) = interp1( ...
                        xUnique, ...
                        yUnique, ...
                        xCommon, ...
                        'linear', ...
                        nan);
                end

                yMean = mean(Yinterp,1,'omitnan');
                yStd  = std(Yinterp,0,1,'omitnan');
                nPts  = sum(~isnan(Yinterp),1);

                yCI = 1.96 * yStd ./ sqrt(nPts);

                fill([xCommon fliplr(xCommon)], ...
                     [yMean-yCI fliplr(yMean+yCI)], ...
                     thisColor, ...
                     'FaceAlpha', alphaCI, ...
                     'EdgeColor', 'none');

                plot(xCommon, yMean, ...
                    'Color', thisColor, ...
                    'LineWidth', 2.5)

                xlim([min(xCommon) max(xCommon)])

            end
        end

        title([thisTitle ' Mean ± 95% CI'])
        xlabel('Time (s)')
        ylabel(mainTitle)
        grid on
        ylim([-6 0])

    end

    sgtitle(sprintf('%s Across Events', mainTitle))

end

%% ============================================================
%% REM STOP vs MefREM STOP ("Wake Events")
%%
%% Creates versions of:
%% 1) Percent time moving
%% 2) Mean ROI motion across events
%% 3) EMG Power across events
%%
%% Uses:
%% remStop
%% mefRemStop
%%
%% Same formatting as before:
%% - no Awake plots
%% - only REM Stop + MefREM Stop
%% - add "Wake Events" to all titles
%% ============================================================

nR = length(remWake);
nM = length(mefRemWake);

roiLabels = strrep(remWake(1).signalNames,'_',' ');

%% ============================================================
%% 1) PERCENT TIME MOVING
%%
%% Assumes:
%% .percentAbove3
%% .percentAbove2
%% .percentAbove1point5
%%
%% Can be Nx16 matrices
%% ============================================================

varsToPlot = {'percentAbove3','percentAbove2','percentAbove1point5'};
titlesToPlot = { ...
    '% Time Moving (z > 3)', ...
    '% Time Moving (z > 2)', ...
    '% Time Moving (z > 1.5)'};

for v = 1:length(varsToPlot)

    varName = varsToPlot{v};

    figure
    tiledlayout(4,4,'TileSpacing','compact')

    for roi = 1:length(roiLabels)

        nexttile
        hold on

        %% REM STOP
        allR = [];

        for r = 1:nR
            if isfield(remWake(r),varName)
                thisVal = remWake(r).(varName);

                if ~isempty(thisVal)
                    if isvector(thisVal)
                        allR = [allR; thisVal(:,roi)];
                    else
                        allR = [allR; thisVal(:,roi)];
                    end
                end
            end
        end

        %% MefREM STOP
        allM = [];

        for r = 1:nM
            if isfield(mefRemWake(r),varName)
                thisVal = mefRemWake(r).(varName);

                if ~isempty(thisVal)
                    if isvector(thisVal)
                        allM = [allM; thisVal(:,roi)];
                    else
                        allM = [allM; thisVal(:,roi)];
                    end
                end
            end
        end

        plotSpread({allR, allM}, ...
            'xNames', {'REM Stop','MefREM Stop'}, ...
            'showMM', 5)

        title(roiLabels{roi})
        ylabel('% Time Moving')
        grid on

        ax = gca;
        h = findobj(ax,'Type','Line');

        for i = 1:length(h)
            if isequal(get(h(i),'Color'), [1 0 0])
                set(h(i),'Color','w','LineWidth',2)
            end
        end
    end

    sgtitle([titlesToPlot{v} ' — Wake Events'])

end


%% ============================================================
%% 2) MEAN ROI MOTION ACROSS WAKE EVENTS
%%
%% Uses:
%% .meanROIMotion
%%
%% cell array
%% col 1 = percent time
%% col 2 = mean ROI zscore
%%
%% Includes:
%% - raw traces
%% - mean ± 95% CI
%% - binned figures at 1%, 5%, 10%
%% ============================================================

%% -------------------------------
%% MAIN FIGURE (raw + mean CI)
%% -------------------------------

figure
tiledlayout(2,2,'TileSpacing','compact')

groups = {remWake, mefRemWake};
groupNames = {'REM Wake','MefREM Wake'};
groupColors = { ...
    [1 0.4 0.4], ...
    [0.2 0.7 0.2]};

for g = 1:2

    dataStruct = groups{g};
    thisTitle = groupNames{g};
    thisColor = groupColors{g};

    %% TOP: individual traces
    nexttile(g)
    hold on

    allX = {};
    allY = {};

    for r = 1:length(dataStruct)
        if isfield(dataStruct(r),'meanROIMotion') && ...
                ~isempty(dataStruct(r).meanROIMotion)

            for k = 1:length(dataStruct(r).meanROIMotion)

                thisData = dataStruct(r).meanROIMotion{k};
                if isempty(thisData)
                    continue
                end

                x = thisData(:,1);
                y = thisData(:,2);

                allX{end+1} = x;
                allY{end+1} = y;

                plot(x,y,...
                    'Color',[thisColor 0.25],...
                    'LineWidth',1.5)
            end
        end
    end

    title(thisTitle)
    xlabel('% Event')
    ylabel('Mean ROI Z-score')
    grid on
    ylim([-1 8])

    %% BOTTOM: mean ± 95% CI
    nexttile(g+2)
    hold on

    if ~isempty(allX)

        minStart = max(cellfun(@(x) min(x), allX));
        maxEnd   = min(cellfun(@(x) max(x), allX));

        xCommon = linspace(minStart,maxEnd,300);
        Yinterp = nan(length(allY),length(xCommon));

        for i = 1:length(allY)

            x = allX{i};
            y = allY{i};

            [xUnique, ia] = unique(x);
            yUnique = y(ia);

            if length(xUnique) < 2
                continue
            end

            Yinterp(i,:) = interp1( ...
                xUnique,yUnique,xCommon,'linear',nan);
        end

        yMean = mean(Yinterp,1,'omitnan');
        yStd  = std(Yinterp,0,1,'omitnan');
        nPts  = sum(~isnan(Yinterp),1);
        yCI   = 1.96 * yStd ./ sqrt(nPts);

        fill([xCommon fliplr(xCommon)], ...
             [yMean-yCI fliplr(yMean+yCI)], ...
             thisColor, ...
             'FaceAlpha',0.18, ...
             'EdgeColor','none');

        plot(xCommon,yMean,...
            'Color',thisColor,...
            'LineWidth',2.5)
    end

    title([thisTitle ' Mean ± 95% CI'])
    xlabel('% Event')
    ylabel('Mean ROI Z-score')
    grid on
    ylim([-1 1.5])

end

sgtitle('Mean ROI Motion Across Wake Events')


%% ============================================================
%% BINNED FIGURES (1%, 5%, 10%)
%%
%% Layout:
%% [REM Wake traces      MefREM Wake traces
%%  REM Wake mean+CI     MefREM Wake mean+CI]
%% ============================================================

binSizes = [1 5 10];

for b = 1:length(binSizes)

    binSize = binSizes(b);

    figure
    tiledlayout(2,2,'TileSpacing','compact')

    for g = 1:2

        dataStruct = groups{g};
        thisTitle = groupNames{g};
        thisColor = groupColors{g};

        allBinned = [];

        %% =====================================
        %% TOP ROW — INDIVIDUAL BINNED TRACES
        %% =====================================

        nexttile(g)
        hold on

        for r = 1:length(dataStruct)

            if isfield(dataStruct(r),'meanROIMotion') && ...
                    ~isempty(dataStruct(r).meanROIMotion)

                for k = 1:length(dataStruct(r).meanROIMotion)

                    thisData = dataStruct(r).meanROIMotion{k};

                    if isempty(thisData)
                        continue
                    end

                    x = thisData(:,1);
                    y = thisData(:,2);

                    edges = 0:binSize:100;
                    centers = edges(1:end-1) + binSize/2;

                    yBin = nan(size(centers));

                    for i = 1:length(centers)
                        idx = x >= edges(i) & x < edges(i+1);

                        if any(idx)
                            yBin(i) = mean(y(idx),'omitnan');
                        end
                    end

                    allBinned(end+1,:) = yBin; %#ok<AGROW>

                    plot(centers,yBin,...
                        'Color',[thisColor 0.22],...
                        'LineWidth',1.2)
                end
            end
        end

        title([thisTitle ' (' num2str(binSize) '% bins)'])
        xlabel('% Event')
        ylabel('Mean ROI Z-score')
        grid on
        ylim([-1 3.5])


        %% =====================================
        %% BOTTOM ROW — MEAN ± 95% CI
        %% =====================================

        nexttile(g+2)
        hold on

        if ~isempty(allBinned)

            yMean = mean(allBinned,1,'omitnan');
            yStd  = std(allBinned,0,1,'omitnan');
            nPts  = sum(~isnan(allBinned),1);
            yCI   = 1.96 * yStd ./ sqrt(nPts);

            fill([centers fliplr(centers)], ...
                 [yMean-yCI fliplr(yMean+yCI)], ...
                 thisColor, ...
                 'FaceAlpha',0.18, ...
                 'EdgeColor','none');

            plot(centers,yMean,...
                'Color',thisColor,...
                'LineWidth',3)
        end

        title([thisTitle ' Mean ± 95% CI'])
        xlabel('% Event')
        ylabel('Mean ROI Z-score')
        grid on
        ylim([-1 1.5])

    end

    sgtitle(['Mean ROI Motion Across Wake Events — ' ...
        num2str(binSize) '% Binning'])

end

%% ============================================================
%% EMG POWER ACROSS WAKE EVENTS — BINNED FIGURES
%%
%% Same format as Mean ROI Motion:
%% - Main figure (raw traces + mean ± 95% CI)
%% - 1% bins
%% - 5% bins
%% - 10% bins
%%
%% Uses:
%% .emgPower
%%
%% Layout:
%% [REM Wake traces      MefREM Wake traces
%%  REM Wake mean+CI     MefREM Wake mean+CI]
%% ============================================================

%% -------------------------------
%% MAIN FIGURE (raw + mean CI)
%% -------------------------------

figure
tiledlayout(2,2,'TileSpacing','compact')

for g = 1:2

    dataStruct = groups{g};
    thisTitle = groupNames{g};
    thisColor = groupColors{g};

    %% TOP ROW — individual traces
    nexttile(g)
    hold on

    allX = {};
    allY = {};

    for r = 1:length(dataStruct)

        if isfield(dataStruct(r),'emgPower') && ...
                ~isempty(dataStruct(r).emgPower)

            for k = 1:length(dataStruct(r).emgPower)

                thisData = dataStruct(r).emgPower{k};

                if isempty(thisData)
                    continue
                end

                x = thisData(:,1);
                y = thisData(:,2);

                allX{end+1} = x;
                allY{end+1} = y;

                plot(x,y,...
                    'Color',[thisColor 0.25],...
                    'LineWidth',1.5)
            end
        end
    end

    title(thisTitle)
    xlabel('% Event')
    ylabel('EMG Power')
    grid on
    ylim([-6 0])

    %% BOTTOM ROW — mean ± 95% CI
    nexttile(g+2)
    hold on

    if ~isempty(allX)

        minStart = max(cellfun(@(x) min(x), allX));
        maxEnd   = min(cellfun(@(x) max(x), allX));

        xCommon = linspace(minStart,maxEnd,300);
        Yinterp = nan(length(allY),length(xCommon));

        for i = 1:length(allY)

            x = allX{i};
            y = allY{i};

            [xUnique, ia] = unique(x);
            yUnique = y(ia);

            if length(xUnique) < 2
                continue
            end

            Yinterp(i,:) = interp1( ...
                xUnique,yUnique,xCommon,'linear',nan);
        end

        yMean = mean(Yinterp,1,'omitnan');
        yStd  = std(Yinterp,0,1,'omitnan');
        nPts  = sum(~isnan(Yinterp),1);
        yCI   = 1.96 * yStd ./ sqrt(nPts);

        fill([xCommon fliplr(xCommon)], ...
             [yMean-yCI fliplr(yMean+yCI)], ...
             thisColor, ...
             'FaceAlpha',0.18, ...
             'EdgeColor','none');

        plot(xCommon,yMean,...
            'Color',thisColor,...
            'LineWidth',2.5)
    end

    title([thisTitle ' Mean ± 95% CI'])
    xlabel('% Event')
    ylabel('EMG Power')
    grid on
    ylim([-6 0])

end

sgtitle('EMG Power Across Wake Events')


%% -------------------------------
%% BINNED FIGURES (1%, 5%, 10%)
%% -------------------------------

binSizes = [1 5 10];

for b = 1:length(binSizes)

    binSize = binSizes(b);

    figure
    tiledlayout(2,2,'TileSpacing','compact')

    for g = 1:2

        dataStruct = groups{g};
        thisTitle = groupNames{g};
        thisColor = groupColors{g};

        allBinned = [];

        %% TOP ROW — INDIVIDUAL BINNED TRACES

        nexttile(g)
        hold on

        for r = 1:length(dataStruct)

            if isfield(dataStruct(r),'emgPower') && ...
                    ~isempty(dataStruct(r).emgPower)

                for k = 1:length(dataStruct(r).emgPower)

                    thisData = dataStruct(r).emgPower{k};

                    if isempty(thisData)
                        continue
                    end

                    x = thisData(:,1);
                    y = thisData(:,2);

                    edges = 0:binSize:100;
                    centers = edges(1:end-1) + binSize/2;

                    yBin = nan(size(centers));

                    for i = 1:length(centers)
                        idx = x >= edges(i) & x < edges(i+1);

                        if any(idx)
                            yBin(i) = mean(y(idx),'omitnan');
                        end
                    end

                    allBinned(end+1,:) = yBin; %#ok<AGROW>

                    plot(centers,yBin,...
                        'Color',[thisColor 0.22],...
                        'LineWidth',1.2)
                end
            end
        end

        title([thisTitle ' (' num2str(binSize) '% bins)'])
        xlabel('% Event')
        ylabel('EMG Power')
        grid on
        ylim([-6 0])

        %% BOTTOM ROW — MEAN ± 95% CI

        nexttile(g+2)
        hold on

        if ~isempty(allBinned)

            yMean = mean(allBinned,1,'omitnan');
            yStd  = std(allBinned,0,1,'omitnan');
            nPts  = sum(~isnan(allBinned),1);
            yCI   = 1.96 * yStd ./ sqrt(nPts);

            fill([centers fliplr(centers)], ...
                 [yMean-yCI fliplr(yMean+yCI)], ...
                 thisColor, ...
                 'FaceAlpha',0.18, ...
                 'EdgeColor','none');

            plot(centers,yMean,...
                'Color',thisColor,...
                'LineWidth',3)
        end

        title([thisTitle ' Mean ± 95% CI'])
        xlabel('% Event')
        ylabel('EMG Power')
        grid on
        ylim([-6 0])

    end

    sgtitle(['EMG Power Across Wake Events — ' ...
        num2str(binSize) '% Binning'])

end

%% ============================================================
%% RESPIRATION RATE (respFreqCentroid)
%%
%% Uses:
%% results(r).respFreqCentroid
%%
%% Format matches previous meanROIMotion / emgPower plots:
%%
%% Figure 1:
%% Full event traces (0–100% event progression)
%%
%% Figure 2:
%% 1% binning
%%
%% Figure 3:
%% 5% binning
%%
%% Figure 4:
%% 10% binning
%%
%% Includes:
%% Awake / REM / MefREM
%%
%% Each cell contains:
%% column 1 = % event progression (0–100)
%% column 2 = respiration frequency centroid
%% ============================================================

groups = {awake, rem, mefRem};
groupNames = {'Awake','REM','MefREM'};
groupColors = { ...
    [0.3 0.5 1], ...
    [1 0.4 0.4], ...
    [0.2 0.7 0.2]};

%% ============================================================
%% FIGURE 1 — FULL EVENT TRACES + MEAN ± 95% CI
%% ============================================================

figure
tiledlayout(2,3,'TileSpacing','compact')

for g = 1:3

    dataStruct = groups{g};
    thisTitle = groupNames{g};
    thisColor = groupColors{g};

    allX = {};
    allY = {};

    %% --------------------------------
    %% TOP ROW: Individual traces
    %% --------------------------------

    nexttile(g)
    hold on

    for r = 1:length(dataStruct)

        if isfield(dataStruct(r),'respFreqCentroid') && ...
                ~isempty(dataStruct(r).respFreqCentroid)

            for k = 1:length(dataStruct(r).respFreqCentroid)

                thisData = dataStruct(r).respFreqCentroid{k};

                if isempty(thisData) || size(thisData,2) < 2
                    continue
                end

                x = thisData(:,1);
                y = thisData(:,2);

                x = x(:);
                y = y(:);

                allX{end+1} = x; %#ok<AGROW>
                allY{end+1} = y; %#ok<AGROW>

                plot(x, y, ...
                    'Color', [thisColor 0.25], ...
                    'LineWidth', 1.5)
            end
        end
    end

    title(thisTitle)
    xlabel('% Event')
    ylabel('Respiration Rate')
    grid on
    xlim([0 100])
    ylim([0 3.5])

    %% --------------------------------
    %% BOTTOM ROW: Mean ± 95% CI
    %% --------------------------------

    nexttile(g+3)
    hold on

    if ~isempty(allX)

        minStart = max(cellfun(@(x) min(x), allX));
        maxEnd   = min(cellfun(@(x) max(x), allX));

        if maxEnd > minStart

            xCommon = linspace(minStart, maxEnd, 300);
            Yinterp = nan(length(allY), length(xCommon));

            for i = 1:length(allY)

                x = allX{i};
                y = allY{i};

                [xUnique, ia] = unique(x);
                yUnique = y(ia);

                if length(xUnique) < 2
                    continue
                end

                Yinterp(i,:) = interp1( ...
                    xUnique, ...
                    yUnique, ...
                    xCommon, ...
                    'linear', ...
                    nan);
            end

            yMean = mean(Yinterp,1,'omitnan');
            yStd  = std(Yinterp,0,1,'omitnan');
            nPts  = sum(~isnan(Yinterp),1);

            yCI = 1.96 * yStd ./ sqrt(nPts);

            fill([xCommon fliplr(xCommon)], ...
                 [yMean-yCI fliplr(yMean+yCI)], ...
                 thisColor, ...
                 'FaceAlpha', 0.18, ...
                 'EdgeColor', 'none');

            plot(xCommon, yMean, ...
                'Color', thisColor, ...
                'LineWidth', 2.5)

            xlim([0 100])
            ylim([0 3.5])

        end
    end

    title([thisTitle ' Mean ± 95% CI'])
    xlabel('% Event')
    ylabel('Respiration Rate')
    grid on

end

sgtitle('Respiration Rate Across Events')


%% ============================================================
%% FIGURES 2–4 — BINNED RESPIRATION RATE
%% (updated: traces on top, mean ± 95% CI below)
%% ============================================================

binSizes = [1 5 10];

for b = 1:length(binSizes)

    binSize = binSizes(b);

    figure
    tiledlayout(2,3,'TileSpacing','compact')

    for g = 1:3

        dataStruct = groups{g};
        thisTitle = groupNames{g};
        thisColor = groupColors{g};

        allBinned = [];

        %% --------------------------------
        %% TOP ROW: individual binned traces
        %% --------------------------------

        nexttile(g)
        hold on

        for r = 1:length(dataStruct)

            if isfield(dataStruct(r),'respFreqCentroid') && ...
                    ~isempty(dataStruct(r).respFreqCentroid)

                for k = 1:length(dataStruct(r).respFreqCentroid)

                    thisData = dataStruct(r).respFreqCentroid{k};

                    if isempty(thisData) || size(thisData,2) < 2
                        continue
                    end

                    x = thisData(:,1);
                    y = thisData(:,2);

                    edges = 0:binSize:100;
                    centers = edges(1:end-1) + binSize/2;

                    yBinned = nan(size(centers));

                    for i = 1:length(centers)

                        idx = x >= edges(i) & x < edges(i+1);

                        if any(idx)
                            yBinned(i) = mean(y(idx),'omitnan');
                        end
                    end

                    allBinned = [allBinned; yBinned]; %#ok<AGROW>

                    plot(centers, yBinned, ...
                        'Color', [thisColor 0.25], ...
                        'LineWidth', 1.2)
                end
            end
        end

        title(thisTitle)
        xlabel('% Event')
        ylabel('Respiration Rate')
        grid on
        xlim([0 100])
        ylim([0 3.5])

        %% --------------------------------
        %% BOTTOM ROW: mean ± 95% CI
        %% --------------------------------

        nexttile(g+3)
        hold on

        if ~isempty(allBinned)

            yMean = mean(allBinned,1,'omitnan');
            yStd  = std(allBinned,0,1,'omitnan');
            nPts  = sum(~isnan(allBinned),1);

            yCI = 1.96 * yStd ./ sqrt(nPts);

            fill([centers fliplr(centers)], ...
                 [yMean-yCI fliplr(yMean+yCI)], ...
                 thisColor, ...
                 'FaceAlpha',0.18, ...
                 'EdgeColor','none');

            plot(centers, yMean, ...
                'Color', thisColor, ...
                'LineWidth', 2.5)
        end

        title([thisTitle ' Mean ± 95% CI'])
        xlabel('% Event')
        ylabel('Respiration Rate')
        grid on
        xlim([0 100])
        ylim([0 3.5])

    end

    sgtitle(sprintf('Respiration Rate — %d%% Binning', binSize))

end

%% ============================================================
%% WAKE EVENTS — RESPIRATION RATE (respFreqCentroid)
%%
%% Uses:
%% remStop(r).respFreqCentroid
%% mefRemStop(r).respFreqCentroid
%%
%% Same format as before:
%% Figure 1 -> full event traces (0–100%)
%% Figure 2 -> 1% bins
%% Figure 3 -> 5% bins
%% Figure 4 -> 10% bins
%%
%% Includes:
%% REM Stop / MefREM Stop
%%
%% Titles include "Wake Events"
%% ============================================================

groups = {remWake, mefRemWake};
groupNames = {'REM Wake','MefREM Wake'};
groupColors = { ...
    [1 0.4 0.4], ...
    [0.2 0.7 0.2]};

%% ============================================================
%% FIGURE 1 — FULL EVENT TRACES + MEAN ± 95% CI
%% ============================================================

figure
tiledlayout(2,2,'TileSpacing','compact')

for g = 1:2

    dataStruct = groups{g};
    thisTitle  = groupNames{g};
    thisColor  = groupColors{g};

    allX = {};
    allY = {};

    %% --------------------------------
    %% TOP ROW: Individual traces
    %% --------------------------------

    nexttile(g)
    hold on

    for r = 1:length(dataStruct)

        if isfield(dataStruct(r),'respFreqCentroid') && ...
                ~isempty(dataStruct(r).respFreqCentroid)

            for k = 1:length(dataStruct(r).respFreqCentroid)

                thisData = dataStruct(r).respFreqCentroid{k};

                if isempty(thisData) || size(thisData,2) < 2
                    continue
                end

                x = thisData(:,1);
                y = thisData(:,2);

                x = x(:);
                y = y(:);

                allX{end+1} = x; %#ok<AGROW>
                allY{end+1} = y; %#ok<AGROW>

                plot(x, y, ...
                    'Color', [thisColor 0.25], ...
                    'LineWidth', 1.5)
            end
        end
    end

    title(thisTitle)
    xlabel('% Event')
    ylabel('Respiration Rate')
    grid on
    xlim([0 100])
    ylim([0 3.5])

    %% --------------------------------
    %% BOTTOM ROW: Mean ± 95% CI
    %% --------------------------------

    nexttile(g+2)
    hold on

    if ~isempty(allX)

        minStart = max(cellfun(@(x) min(x), allX));
        maxEnd   = min(cellfun(@(x) max(x), allX));

        if maxEnd > minStart

            xCommon = linspace(minStart, maxEnd, 300);
            Yinterp = nan(length(allY), length(xCommon));

            for i = 1:length(allY)

                x = allX{i};
                y = allY{i};

                [xUnique, ia] = unique(x);
                yUnique = y(ia);

                if length(xUnique) < 2
                    continue
                end

                Yinterp(i,:) = interp1( ...
                    xUnique, ...
                    yUnique, ...
                    xCommon, ...
                    'linear', ...
                    nan);
            end

            yMean = mean(Yinterp,1,'omitnan');
            yStd  = std(Yinterp,0,1,'omitnan');
            nPts  = sum(~isnan(Yinterp),1);

            yCI = 1.96 * yStd ./ sqrt(nPts);

            fill([xCommon fliplr(xCommon)], ...
                 [yMean-yCI fliplr(yMean+yCI)], ...
                 thisColor, ...
                 'FaceAlpha', 0.18, ...
                 'EdgeColor', 'none');

            plot(xCommon, yMean, ...
                'Color', thisColor, ...
                'LineWidth', 2.5)

            xlim([0 100])
            ylim([0 3.5])

        end
    end

    title([thisTitle ' Mean ± 95% CI'])
    xlabel('% Event')
    ylabel('Respiration Rate')
    grid on

end

sgtitle('Respiration Rate Across Wake Events')


%% ============================================================
%% FIGURES 2–4 — WAKE EVENTS BINNED RESPIRATION RATE
%% (updated: traces on top, mean ± 95% CI below)
%% ============================================================

binSizes = [1 5 10];

for b = 1:length(binSizes)

    binSize = binSizes(b);

    figure
    tiledlayout(2,2,'TileSpacing','compact')

    for g = 1:2

        dataStruct = groups{g};
        thisTitle  = groupNames{g};
        thisColor  = groupColors{g};

        allBinned = [];

        %% --------------------------------
        %% TOP ROW: individual binned traces
        %% --------------------------------

        nexttile(g)
        hold on

        for r = 1:length(dataStruct)

            if isfield(dataStruct(r),'respFreqCentroid') && ...
                    ~isempty(dataStruct(r).respFreqCentroid)

                for k = 1:length(dataStruct(r).respFreqCentroid)

                    thisData = dataStruct(r).respFreqCentroid{k};

                    if isempty(thisData) || size(thisData,2) < 2
                        continue
                    end

                    x = thisData(:,1);
                    y = thisData(:,2);

                    edges = 0:binSize:100;
                    centers = edges(1:end-1) + binSize/2;

                    yBinned = nan(size(centers));

                    for i = 1:length(centers)

                        idx = x >= edges(i) & x < edges(i+1);

                        if any(idx)
                            yBinned(i) = mean(y(idx),'omitnan');
                        end
                    end

                    allBinned = [allBinned; yBinned]; %#ok<AGROW>

                    plot(centers, yBinned, ...
                        'Color', [thisColor 0.25], ...
                        'LineWidth', 1.2)
                end
            end
        end

        title(thisTitle)
        xlabel('% Event')
        ylabel('Respiration Rate')
        grid on
        xlim([0 100])
        ylim([0 3.5])

        %% --------------------------------
        %% BOTTOM ROW: mean ± 95% CI
        %% --------------------------------

        nexttile(g+2)
        hold on

        if ~isempty(allBinned)

            yMean = mean(allBinned,1,'omitnan');
            yStd  = std(allBinned,0,1,'omitnan');
            nPts  = sum(~isnan(allBinned),1);

            yCI = 1.96 * yStd ./ sqrt(nPts);

            fill([centers fliplr(centers)], ...
                 [yMean-yCI fliplr(yMean+yCI)], ...
                 thisColor, ...
                 'FaceAlpha',0.18, ...
                 'EdgeColor','none');

            plot(centers, yMean, ...
                'Color', thisColor, ...
                'LineWidth', 2.5)
        end

        title([thisTitle ' Mean ± 95% CI'])
        xlabel('% Event')
        ylabel('Respiration Rate')
        grid on
        xlim([0 100])
        ylim([0 3.5])

    end

    sgtitle(sprintf('Respiration Rate Across Wake Events — %d%% Binning', binSize))

end

%% ============================================================
%% MEAN TEMPERATURE COMPARISON + REM HISTOGRAM
%%
%% Top:
%% plotSpread comparison of Awake / REM / MefREM
%%
%% Bottom:
%% Histogram of REM event temperatures
%% ============================================================

figure
tiledlayout(2,1,'TileSpacing','compact')

groups = {awake, rem, mefRem};
groupNames = {'Awake','REM','MefREM'};

plotData = cell(1,3);

for g = 1:3

    dataStruct = groups{g};
    tempVals = [];

    for r = 1:length(dataStruct)

        if isfield(dataStruct(r),'meanTemp') && ...
                ~isempty(dataStruct(r).meanTemp)

            thisTemp = dataStruct(r).meanTemp;

            if iscell(thisTemp)

                for k = 1:length(thisTemp)

                    if isempty(thisTemp{k})
                        continue
                    end

                    vals = thisTemp{k};
                    vals = vals(:);
                    vals = vals(isfinite(vals));

                    tempVals = [tempVals; vals]; %#ok<AGROW>
                end

            else
                vals = thisTemp(:);
                vals = vals(isfinite(vals));
                tempVals = [tempVals; vals]; %#ok<AGROW>
            end
        end
    end

    plotData{g} = tempVals;

end

%% ============================================================
%% TOP: plotSpread
%% ============================================================

nexttile
hold on

plotSpread(plotData, ...
    'xNames', groupNames, ...
    'showMM', 5, ...
    'distributionMarkers', 'o', ...
    'distributionColors', {'b','r','g'});

ylabel('Mean Temperature')
title('Mean Temperature Across Sleep States')
grid on
box on

% make mean + SD white
h = findobj(gca,'Type','Line');

for i = 1:length(h)
    if strcmp(get(h(i),'Marker'),'+') || ...
       strcmp(get(h(i),'Marker'),'none')
        set(h(i),'Color','w','LineWidth',1.5)
    end
end

%% ============================================================
%% BOTTOM: REM histogram
%% ============================================================

nexttile
hold on

remTemps = [plotData{2}; plotData{3}]; % REM groups

histogram(remTemps, ...
    'BinMethod','auto', ...
    'Normalization','count', ...
    'BinEdges',72:.5:88)

xlabel('Mean Temperature')
ylabel('Number of REM Events')
title('REM Event Temperature Distribution')
grid on
box on

%% ============================================================
%% LEFT–RIGHT ROI CORRELATION COMPARISON
%%
%% Uses:
%% results(r).signalNames
%% results(r).zScore
%%
%% signalNames:
%% cell array of ROI names
%%
%% Example:
%% Left_Whisker
%% Right_Whisker
%%
%% zScore:
%% rows = time
%% cols = ROI signals
%%
%% Goal:
%% Compare left-right correlation across:
%% Awake / REM / MefREM
%%
%% Every left-right pair from every recording
%% becomes its own point on the plot
%%
%% Higher correlation = more symmetric movement
%% Lower correlation = more asymmetric movement
%% ============================================================

figure
hold on

groups = {awake, rem, mefRem};
groupNames = {'Awake','REM','MefREM'};

plotData = cell(1,3);

for g = 1:3

    dataStruct = groups{g};
    corrVals = [];

    for r = 1:length(dataStruct)

        if ~isfield(dataStruct(r),'signalNames') || ...
           ~isfield(dataStruct(r),'zScore') || ...
           isempty(dataStruct(r).signalNames) || ...
           isempty(dataStruct(r).zScore)

            continue
        end

        signalNames = dataStruct(r).signalNames;
        zData = dataStruct(r).zScore;

        if isempty(signalNames) || isempty(zData)
            continue
        end

        % ensure signalNames is cell
        if ~iscell(signalNames)
            continue
        end

        for i = 1:length(signalNames)

            thisName = signalNames{i};

            if startsWith(thisName,'Left_')

                baseName = erase(thisName,'Left_');
                rightName = ['Right_' baseName];

                rightIdx = find(strcmp(signalNames,rightName),1);

                if isempty(rightIdx)
                    continue
                end

                if i > size(zData,2) || rightIdx > size(zData,2)
                    continue
                end

                L = zData(:,i);
                R = zData(:,rightIdx);

                validIdx = isfinite(L) & isfinite(R);

                if sum(validIdx) < 10
                    continue
                end

                rVal = corr(L(validIdx),R(validIdx));

                if isfinite(rVal)
                    corrVals = [corrVals; rVal]; %#ok<AGROW>
                end
            end
        end
    end

    plotData{g} = corrVals;

end

%% ============================================================
%% PLOT USING plotSpread
%% ============================================================

plotSpread(plotData, ...
    'xNames', groupNames, ...
    'showMM', 5, ... % mean ± SD
    'distributionMarkers', 'o', ...
    'distributionColors', {'b','r','g'});

ylabel('Left–Right ROI Correlation (r)')
title('Left–Right Movement Symmetry Across Sleep States')
ylim([-1 1])

grid on
box on

%% Make mean ± SD white

h = findobj(gca,'Type','Line');

for i = 1:length(h)

    if strcmp(get(h(i),'Marker'),'+') || ...
       strcmp(get(h(i),'Marker'),'none')

        set(h(i),'Color','w','LineWidth',1.5)
    end
end

%% ============================================================
%% LEFT–RIGHT LAGGED CROSS-CORRELATION (IMPROVED)
%%
%% Improvements:
%% 1. detrend signals first
%% 2. only analyze active movement epochs
%%    where abs(zscore) > threshold
%%
%% This avoids:
%% - quiet period domination
%% - zero-lag baseline artifacts
%%
%% Output:
%% Figure 1: Peak cross-correlation
%% Figure 2: Lag at peak correlation
%% Figure 3: Twitch coincidence
%% ============================================================

groups = {awake, rem, mefRem};
groupNames = {'Awake','REM','MefREM'};

peakCorrData = cell(1,3);
peakLagData  = cell(1,3);
twitchCoincidenceData = cell(1,3);

maxLag = 50;          % frames
threshold = 2;        % active movement threshold

for g = 1:3

    dataStruct = groups{g};

    peakCorrVals = [];
    peakLagVals  = [];
    coincidenceVals = [];

    for r = 1:length(dataStruct)

        if ~isfield(dataStruct(r),'signalNames') || ...
           ~isfield(dataStruct(r),'zScore') || ...
           isempty(dataStruct(r).signalNames) || ...
           isempty(dataStruct(r).zScore)
            continue
        end

        signalNames = dataStruct(r).signalNames;
        zData = dataStruct(r).zScore;

        for i = 1:length(signalNames)

            thisName = signalNames{i};

            if startsWith(thisName,'Left_')

                baseName = erase(thisName,'Left_');
                rightName = ['Right_' baseName];

                rightIdx = find(strcmp(signalNames,rightName),1);

                if isempty(rightIdx)
                    continue
                end

                if i > size(zData,2) || rightIdx > size(zData,2)
                    continue
                end

                L = zData(:,i);
                R = zData(:,rightIdx);

                valid = isfinite(L) & isfinite(R);

                if sum(valid) < 20
                    continue
                end

                L = L(valid);
                R = R(valid);

                %% --------------------------------
                %% detrend first
                %% --------------------------------

                L = detrend(L);
                R = detrend(R);

                %% --------------------------------
                %% keep only active movement epochs
                %% --------------------------------

                activeIdx = abs(L) > threshold | abs(R) > threshold;

                if sum(activeIdx) < 20
                    continue
                end

                L_active = L(activeIdx);
                R_active = R(activeIdx);

                %% --------------------------------
                %% lagged cross-correlation
                %% --------------------------------

                [xc, lags] = xcorr(L_active, R_active, maxLag, 'coeff');

                [peakVal, idx] = max(xc);
                peakLag = lags(idx);

                if isfinite(peakVal)
                    peakCorrVals(end+1,1) = peakVal; %#ok<AGROW>
                    peakLagVals(end+1,1)  = peakLag; %#ok<AGROW>
                end

                %% --------------------------------
                %% twitch coincidence
                %% --------------------------------

                leftTwitch  = L > threshold;
                rightTwitch = R > threshold;

                nLeft = sum(leftTwitch);

                if nLeft > 0

                    coincidence = ...
                        sum(leftTwitch & rightTwitch) / nLeft;

                    if isfinite(coincidence)
                        coincidenceVals(end+1,1) = coincidence; %#ok<AGROW>
                    end
                end

            end
        end
    end

    peakCorrData{g} = peakCorrVals;
    peakLagData{g}  = peakLagVals;
    twitchCoincidenceData{g} = coincidenceVals;

end


%% ============================================================
%% FIGURE 1 — Peak Cross-Correlation
%% ============================================================

figure
hold on

plotSpread(peakCorrData, ...
    'xNames', groupNames, ...
    'showMM', 5, ...
    'distributionMarkers', 'o');

ylabel('Peak Cross-Correlation')
title('Left–Right Peak Lagged Cross-Correlation')
grid on
box on


%% ============================================================
%% FIGURE 2 — Lag at Peak Correlation
%% ============================================================

figure
hold on

plotSpread(peakLagData, ...
    'xNames', groupNames, ...
    'showMM', 5, ...
    'distributionMarkers', 'o');

ylabel('Lag at Peak Correlation (frames)')
title('Lag of Maximum Left–Right Correlation')
grid on
box on


%% ============================================================
%% FIGURE 3 — Twitch Coincidence
%% ============================================================

figure
hold on

plotSpread(twitchCoincidenceData, ...
    'xNames', groupNames, ...
    'showMM', 5, ...
    'distributionMarkers', 'o');

ylabel('P(Right Twitch | Left Twitch)')
title('Left–Right Twitch Coincidence')
ylim([0 1])

grid on
box on

%% ============================================================
%% LEFT vs RIGHT TWITCH ONSET DIFFERENCE
%%
%% Goal:
%% Detect discrete twitch onsets using z-score threshold crossing
%% and compare left/right timing differences directly.
%%
%% This is much better for asynchronous REM twitches than
%% simple correlation because correlation is dominated by baseline.
%%
%% Uses:
%% results(r).signalNames
%% results(r).zScore
%%
%% Assumes:
%% - zScore columns correspond to signalNames
%% - Left/right pairs are named:
%%     Left_xxx
%%     Right_xxx
%%
%% Output:
%% PlotSpread figure of absolute onset delay (frames)
%% between left/right twitches for:
%% Awake / REM / MefREM
%%
%% Smaller values = more synchronous
%% Larger values = more asymmetric
%% ============================================================

threshold = 2.5;   % z-score threshold for twitch detection
minGap    = 5;     % minimum separation between twitches (frames)

groups = {awake, rem, mefRem};
groupNames = {'Awake','REM','MefREM'};

allDelays = cell(1,3);

for g = 1:3

    dataStruct = groups{g};
    groupDelays = [];

    for r = 1:length(dataStruct)

        if ~isfield(dataStruct(r),'signalNames') || ...
           ~isfield(dataStruct(r),'zScore') || ...
           isempty(dataStruct(r).signalNames) || ...
           isempty(dataStruct(r).zScore)
            continue
        end

        names = dataStruct(r).signalNames;
        Z = dataStruct(r).zScore;

        for i = 1:length(names)

            thisName = names{i};

            if startsWith(thisName,'Left_')

                baseName = erase(thisName,'Left_');
                rightName = ['Right_' baseName];

                j = find(strcmp(names,rightName),1);

                if isempty(j)
                    continue
                end

                leftTrace  = Z(:,i);
                rightTrace = Z(:,j);

                %% ---------------------------------
                %% Detect threshold crossing onsets
                %% ---------------------------------

                leftBinary = leftTrace > threshold;
                rightBinary = rightTrace > threshold;

                leftOnsets = find(diff([0; leftBinary]) == 1);
                rightOnsets = find(diff([0; rightBinary]) == 1);

                %% enforce minimum spacing
                if ~isempty(leftOnsets)
                    leftOnsets = leftOnsets( ...
                        [true; diff(leftOnsets) > minGap]);
                end

                if ~isempty(rightOnsets)
                    rightOnsets = rightOnsets( ...
                        [true; diff(rightOnsets) > minGap]);
                end

                if isempty(leftOnsets) || isempty(rightOnsets)
                    continue
                end

                %% ---------------------------------
                %% Match nearest left/right twitches
                %% ---------------------------------

                for k = 1:length(leftOnsets)

                    d = abs(rightOnsets - leftOnsets(k));
                    nearestDelay = min(d);

                    groupDelays(end+1,1) = nearestDelay; %#ok<AGROW>
                end

            end
        end
    end

    allDelays{g} = groupDelays;

end


%% ============================================================
%% PLOT — DISTRIBUTION OF TWITCH ONSET DELAYS
%% ============================================================

figure
hold on

plotSpread(allDelays, ...
    'xNames', groupNames, ...
    'showMM', 5, ...   % mean ± SD
    'distributionMarkers', 'o')

ylabel('Left-Right Twitch Onset Delay (frames)')
title('Left vs Right Twitch Timing Asymmetry')

set(gca,'FontSize',12)
grid on


%% ============================================================
%% OPTIONAL:
%% Convert to seconds if frame rate known
%%
%% Example:
%% fps = 30;
%% delaySeconds = delayFrames / fps;
%%
%% Then label:
%% ylabel('Twitch Onset Delay (s)')
%% ============================================================

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function results = processExcelSegmentsTDMS_SF(excelPath)

T = readtable(excelPath,'Sheet',2,'ReadVariableNames',false);
nRec = height(T);

for r = 1:nRec
    
    tdmsPath = T{r,1}{1};
    
    segments = table2array(T(r,2:end));
    segments = segments(~isnan(segments));
    segments = reshape(segments,2,[])';
    
    [coeff,score,latent,explained,corrMatrix,signalNames,nSignals] = ...
        computeDigitalSignalsFromSegments_SF(tdmsPath, segments);
    
    results(r).coeff = coeff;
    results(r).explained = explained;
    results(r).corrMatrix = corrMatrix;
    results(r).score = score;
    results(r).signalNames = signalNames;
    results(r).nSignals = nSignals;
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = processExcelSegmentsROISelect_SF(excelPath, roiPath, wakeTF)

T = readtable(excelPath,'Sheet',2,'ReadVariableNames',false);
if strcmp(T{1,1}{1},'file name')
    T = T(2:end,:);
end
nRec = height(T);

for r = 1:nRec
    
    tdmsPath = T{r,1}{1};
    segments = table2array(T(r,2:end));
    segments = segments(~isnan(segments));
    segments = reshape(segments,2,[])';

    if wakeTF 
        segments = [segments segments(:,2)+30];
    end

    % ==============================
    % 1) Convert TDMS path → video path
    % ==============================
    [folder, name, ~] = fileparts(tdmsPath);
    
    % Change extension as needed (.avi, .mp4, etc.)
    videoPath = fullfile(folder, [name '.mp4']);  % <-- adjust if needed
    
    % ==============================
    % 2) Extract / Load ROI motion
    % ==============================
    roiMotion = extractROIMotionFromVideo_fast(videoPath, excelPath, roiPath, segments);
    roiMotionBaseline = extractROIMotionFromVideo_baseline(videoPath, excelPath, roiPath, segments);

    if wakeTF 
        segments = segments(:,[2 3]);
    end
    
    % Remove time field if present
    if isfield(roiMotion,'time')
        timeVec = roiMotion.time;
        roiMotion = rmfield(roiMotion,'time');
    end
    if isfield(roiMotionBaseline,'time')
        timeVecBaseline = roiMotionBaseline.time;
        roiMotionBaseline = rmfield(roiMotionBaseline,'time');
    end
    
    % Enforce consistent ordering across recordings
    signalNames = sort(fieldnames(roiMotion));
    signalNamesBaseline = sort(fieldnames(roiMotionBaseline));

    % % Remove whisker fields
    % idx = cellfun(@(s) contains(s,'whisker','IgnoreCase',true), signalNames);
    % signalNames(idx) = [];
    % idx = cellfun(@(s) contains(s,'whisker','IgnoreCase',true), signalNamesBaseline);
    % signalNamesBaseline(idx) = [];

    nSignals = numel(signalNames);
    nSignalsBaseline = numel(signalNamesBaseline);
    
    % ==============================
    % 3) Build matrix (UNCHANGED LOGIC)
    % ==============================
    L = min(structfun(@length, roiMotion));
    X = zeros(L,nSignals);
    M = min(structfun(@length, roiMotionBaseline));
    Y = zeros(M,nSignalsBaseline);
    
    for i = 1:nSignals
        X(:,i) = roiMotion.(signalNames{i})(1:L);
    end
    for i = 1:nSignalsBaseline
        Y(:,i) = roiMotionBaseline.(signalNamesBaseline{i})(1:M);
    end

    % ====== Preprocessing (IMPORTANT) ======
    X = fillmissing(X,'linear');
    X = zscore(X);   % <<< critical for ROI comparisons
    Y = fillmissing(Y,'linear');
    Y = zscore(Y);   % <<< critical for ROI comparisons
    
    % ==============================
    % 4) PCA (UNCHANGED)
    % ==============================
    [coeff,score,latent,~,explained] = pca(X);
    [coeffB,scoreB,latentB,~,explainedB] = pca(Y);
    
    % ==============================
    % 5) Correlation (UNCHANGED)
    % ==============================
    corrMatrix = corrcoef(X);
    corrMatrixBaseline = corrcoef(Y);

    % ==============================
    % 6) Segment lengths
    % ==============================
    segmentLengths = [];
    for n = 1:size(segments,1)
        segmentLengths(end+1) = segments(n,2) - segments(n,1);
    end
    
    % ==============================
    % 7) Section segments and process
    % ==============================
    % if ~isfield(roiMotion, 'time')
    %     roiMotion.time = [];
    %     for j = 1:size(segments,1)
    %         roiMotion.time = [roiMotion.time segments(j,1):0.0167:segments(j,2)];
    %     end
    % end
    segmentIdx = find(diff(timeVec)>1);
    segmentIdx = [0 segmentIdx length(timeVec)];
    segmentIdxBaseline = find(diff(timeVecBaseline)>1);
    segmentIdxBaseline = [0 segmentIdxBaseline length(timeVecBaseline)];
    emgTrigAvgLength = 20;
    [p, n, ~] = fileparts(tdmsPath);
    procDataPath = fullfile(p,[n '_ProcData.mat']);
    for n = 1:length(segmentIdx)-1
        currSegIdx = segmentIdx(n)+1:segmentIdx(n+1);
        currSegIdxBaseline = segmentIdxBaseline(n)+1:segmentIdxBaseline(n+1);
        SegX = X(currSegIdx,:);
        SegY = Y(currSegIdxBaseline,:);
        percentAbove3(n,:) = 100 * (sum(SegX > 3, 1) ./ size(SegX, 1));
        percentAbove2(n,:) = 100 * (sum(SegX > 2, 1) ./ size(SegX, 1));
        percentAbove1point5(n,:) = 100 * (sum(SegX > 1.5, 1) ./ size(SegX, 1));
        percentAbove3B(n,:) = 100 * (sum(SegY > 3, 1) ./ size(SegY, 1));
        percentAbove2B(n,:) = 100 * (sum(SegY > 2, 1) ./ size(SegY, 1));
        percentAbove1point5B(n,:) = 100 * (sum(SegY > 1.5, 1) ./ size(SegY, 1));
        meanROIMotion{n} = [linspace(0,100,size(SegX,1))' , mean(SegX,2)];
        meanROIMotionB{n} = [linspace(0,100,size(SegY,1))' , mean(SegY,2)];
        [emgPowerVec,emgPowerStartVec,emgPowerEndVec] = loadProcData(procDataPath,segments(n,1),segments(n,2),emgTrigAvgLength);
        emgPower{n} = [linspace(0,100,length(emgPowerVec))' , emgPowerVec];
        emgPowerStart{n} = [linspace(-emgTrigAvgLength/2,emgTrigAvgLength/2,length(emgPowerStartVec))' , emgPowerStartVec];
        emgPowerEnd{n} = [linspace(-emgTrigAvgLength/2,emgTrigAvgLength/2,length(emgPowerEndVec))' , emgPowerEndVec];
        [respFreqCentroidVec,meanTempVal] = respirationSpectrogramPlot_SF(tdmsPath,segments(n,1),segments(n,2));
        respFreqCentroid{n} = [linspace(0,100,length(respFreqCentroidVec))' , respFreqCentroidVec'];
        meanTemp{n} = meanTempVal;
    end
    




    % ==============================
    % 8) Store (UNCHANGED)
    % ==============================
    results(r).zScore = X;
    results(r).coeff = coeff;
    results(r).explained = explained;
    results(r).corrMatrix = corrMatrix;
    results(r).score = score;
    results(r).zScoreB = Y;
    results(r).coeffB = coeffB;
    results(r).explainedB = explainedB;
    results(r).corrMatrixB = corrMatrixBaseline;
    results(r).scoreB = scoreB;
    results(r).signalNames = signalNames;
    results(r).signalNamesB = signalNamesBaseline;
    results(r).nSignals = nSignals;
    results(r).nSignalsB = nSignalsBaseline;
    results(r).segmentLengths = segmentLengths;
    results(r).percentAbove3 = percentAbove3;
    results(r).percentAbove2 = percentAbove2;
    results(r).percentAbove1point5 = percentAbove1point5;
    results(r).meanROIMotion = meanROIMotion;
    results(r).percentAbove3B = percentAbove3B;
    results(r).percentAbove2B = percentAbove2B;
    results(r).percentAbove1point5B = percentAbove1point5B;
    results(r).meanROIMotionB = meanROIMotionB;
    results(r).emgPower = emgPower;
    results(r).emgPowerStart = emgPowerStart;
    results(r).emgPowerEnd = emgPowerEnd;
    results(r).respFreqCentroid = respFreqCentroid;
    results(r).meanTemp = meanTemp;
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function results = processExcelSegmentsROISelectStopEvents_SF(excelPath, roiPath)
% 
% T = readtable(excelPath,'Sheet',2,'ReadVariableNames',false);
% if strcmp(T{1,1}{1},'file name')
%     T = T(2:end,:);
% end
% nRec = height(T);
% 
% for r = 1:nRec
% 
%     tdmsPath = T{r,1}{1};
%     segments = table2array(T(r,2:end));
%     segments = segments(~isnan(segments));
%     segments = reshape(segments,2,[])';
% 
%     % ==============================
%     % 1) Convert TDMS path → video path
%     % ==============================
%     [folder, name, ~] = fileparts(tdmsPath);
% 
%     % Change extension as needed (.avi, .mp4, etc.)
%     videoPath = fullfile(folder, [name '.mp4']);  % <-- adjust if needed
% 
%     % ==============================
%     % 2) Extract / Load ROI motion
%     % ==============================
%     roiMotion = extractROIMotionFromVideo_fast(videoPath, excelPath, roiPath, segments, true);
%     roiMotionBaseline = extractROIMotionFromVideo_baseline(videoPath, excelPath, roiPath, segments, true);
% 
%     % Remove time field if present
%     if isfield(roiMotion,'time')
%         timeVec = roiMotion.time;
%         roiMotion = rmfield(roiMotion,'time');
%     end
%     if isfield(roiMotionBaseline,'time')
%         timeVecBaseline = roiMotionBaseline.time;
%         roiMotionBaseline = rmfield(roiMotionBaseline,'time');
%     end
% 
%     % Enforce consistent ordering across recordings
%     signalNames = sort(fieldnames(roiMotion));
% 
%     % % Remove whisker fields
%     % idx = cellfun(@(s) contains(s,'whisker','IgnoreCase',true), signalNames);
%     % signalNames(idx) = [];
% 
%     nSignals = numel(signalNames);
% 
%     % ==============================
%     % 3) Build matrix (UNCHANGED LOGIC)
%     % ==============================
%     L = min(structfun(@length, roiMotion));
%     X = zeros(L,nSignals);
% 
%     for i = 1:nSignals
%         X(:,i) = roiMotion.(signalNames{i})(1:L);
%     end
% 
%     % ====== Preprocessing (IMPORTANT) ======
%     X = fillmissing(X,'linear');
%     X = zscore(X);   % <<< critical for ROI comparisons
% 
%     % ==============================
%     % 4) PCA (UNCHANGED)
%     % ==============================
%     [coeff,score,latent,~,explained] = pca(X);
% 
%     % ==============================
%     % 5) Correlation (UNCHANGED)
%     % ==============================
%     corrMatrix = corrcoef(X);
% 
%     % ==============================
%     % 6) Segment lengths
%     % ==============================
%     segmentLengths = [];
%     for n = 1:size(segments,1)
%         segmentLengths(end+1) = segments(n,2) - segments(n,1);
%     end
% 
%     % ==============================
%     % 7) Section segments and process
%     % ==============================
%     % if ~isfield(roiMotion, 'time')
%     %     roiMotion.time = [];
%     %     for j = 1:size(segments,1)
%     %         roiMotion.time = [roiMotion.time segments(j,1):0.0167:segments(j,2)];
%     %     end
%     % end
%     segmentIdx = find(diff(timeVec)>1);
%     segmentIdx = [0 segmentIdx length(timeVec)];
%     emgTrigAvgLength = 20;
%     [p, n, ~] = fileparts(tdmsPath);
%     procDataPath = fullfile(p,[n '_ProcData.mat']);
%     for n = 1:length(segmentIdx)-1
%         currSegIdx = segmentIdx(n)+1:segmentIdx(n+1);
%         SegX = X(currSegIdx,:);
%         percentAbove3(n,:) = 100 * (sum(SegX > 3, 1) ./ size(SegX, 1));
%         percentAbove2(n,:) = 100 * (sum(SegX > 2, 1) ./ size(SegX, 1));
%         percentAbove1point5(n,:) = 100 * (sum(SegX > 1.5, 1) ./ size(SegX, 1));
%         meanROIMotion{n} = [linspace(0,100,size(SegX,1))' , mean(SegX,2)];
%         [emgPowerVec,emgPowerStartVec,emgPowerEndVec] = loadProcData(procDataPath,segments(n,1),segments(n,2),emgTrigAvgLength);
%         emgPower{n} = [linspace(0,100,length(emgPowerVec))' , emgPowerVec];
%         emgPowerStart{n} = [linspace(-emgTrigAvgLength/2,emgTrigAvgLength/2,length(emgPowerStartVec))' , emgPowerStartVec];
%         emgPowerEnd{n} = [linspace(-emgTrigAvgLength/2,emgTrigAvgLength/2,length(emgPowerEndVec))' , emgPowerEndVec];
%         [respFreqCentroidVec,meanTempVal] = respirationSpectrogramPlot_SF(tdmsPath,segments(n,1),segments(n,2));
%         respFreqCentroid{n} = [linspace(0,100,length(respFreqCentroidVec))' , respFreqCentroidVec'];
%         meanTemp{n} = meanTempVal;
%     end
% 
% 
% 
% 
% 
%     % ==============================
%     % 8) Store (UNCHANGED)
%     % ==============================
%     results(r).zScore = X;
%     results(r).coeff = coeff;
%     results(r).explained = explained;
%     results(r).corrMatrix = corrMatrix;
%     results(r).score = score;
%     results(r).signalNames = signalNames;
%     results(r).nSignals = nSignals;
%     results(r).segmentLengths = segmentLengths;
%     results(r).percentAbove3 = percentAbove3;
%     results(r).percentAbove2 = percentAbove2;
%     results(r).percentAbove1point5 = percentAbove1point5;
%     results(r).meanROIMotion = meanROIMotion;
%     results(r).emgPower = emgPower;
%     results(r).emgPowerStart = emgPowerStart;
%     results(r).emgPowerEnd = emgPowerEnd;
%     results(r).respFreqCentroid = respFreqCentroid;
%     results(r).meanTemp = meanTemp;
% 
% end
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coeff,score,latent,explained,corrMatrix,signalNames,nSignals] = ...
    computeDigitalSignalsFromSegments_SF(tdmsPath, segments)

tdmsDataStruct = openTDMS(tdmsPath);

Fs = str2double(tdmsDataStruct.CameraFrameratePerSecond);

signalsStruct = tdmsDataStruct.Digital_Data;
signalNames = fieldnames(signalsStruct);

% % Remove whisker fields
% idx = cellfun(@(s) contains(s,'whisker','IgnoreCase',true), signalNames);
% signalNames(idx) = [];

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

% ====== Preprocessing (IMPORTANT) ======
X = fillmissing(X,'linear');
X = zscore(X);   % <<< critical for ROI comparisons

%% PCA
[coeff,score,latent,~,explained] = pca(X);

%% Correlation
corrMatrix = corrcoef(X);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotFirstFrameWithROIs(excelPath, roiPath)

T = readtable(excelPath,'Sheet',2,'ReadVariableNames',false);
if strcmp(T{1,1}{1},'file name')
    T = T(2:end,:);
end

% Get first file
tdmsPath = T{1,1}{1};
[folder,name,~] = fileparts(tdmsPath);
videoPath = fullfile(folder,[name '.mp4']); % adjust extension if needed

v = VideoReader(videoPath);

% Read first frame
frame = readFrame(v);

figure
imshow(frame)
title('ROI Locations (First Frame)')
hold on

% ==============================
% Load or draw ROIs
% ==============================
if exist("roiPath","var") && ~isempty(roiPath)
    
    load(roiPath);
    roiMasks  = roiStruct.roiMasks;
    roiLabels = roiStruct.roiLabels;
    
else
    % Let user draw (same behavior as extractor)
    title('Draw ROIs → Double-click → Press Enter when done')
    
    roiMasks = {};
    roiLabels = {};
    roiCount = 0;
    
    while true
        roi = drawrectangle('Color','r');
        if isempty(roi)
            break;
        end
        
        roiCount = roiCount + 1;
        
        label = input(sprintf('Enter label for ROI %d: ', roiCount), 's');
        roiMasks{roiCount} = createMask(roi);
        roiLabels{roiCount} = label;
        
        choice = input('Add another ROI? (y/n): ','s');
        if lower(choice) ~= 'y'
            break;
        end
    end
end

% ==============================
% Overlay ROI outlines
% ==============================
for r = 1:length(roiMasks)
    
    B = bwboundaries(roiMasks{r});
    
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1.5)
    end
    
    % Label position (center of mask)
    [y,x] = find(roiMasks{r});
    % text(mean(x), mean(y), roiLabels{r}, ...
        % 'Color','y','FontSize',10,'FontWeight','bold')
end

hold off

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [emgPowerFull,emgPowerStart,emgPowerEnd] = loadProcData(procDataFileID,tStart,tEnd,emgTrigAvgLength)
%% Load ProcData
S = load(procDataFileID,'ProcData');
ProcData = S.ProcData;

%% Sampling rate
Fs = [];
if isfield(ProcData,'notes') && isfield(ProcData.notes,'dsFs')
    Fs = ProcData.notes.dsFs;
end
if isempty(Fs) || ~isfinite(Fs)
    Fs = 60;
end

%% Reference length
if isfield(ProcData,'ECoG_DS_norm')
    refLen = numel(ProcData.ECoG_DS_norm);
elseif isfield(ProcData,'ECoG_DS')
    refLen = numel(ProcData.ECoG_DS);
else
    error('Cannot determine trial duration.');
end

trialDur = refLen / Fs;
tVec = (0:refLen-1)/Fs;
tMask = tVec >= tStart & tVec <= tEnd;
tMaskStart = tVec >= (tStart-(emgTrigAvgLength/2)) & tVec <= (tStart+(emgTrigAvgLength/2));
tMaskEnd = tVec >= (tEnd-(emgTrigAvgLength/2)) & tVec <= (tEnd+(emgTrigAvgLength/2));
tPlot = tVec(tMask);

% %% Force
% if isfield(ProcData,'forceSensor_norm')
%     force = ProcData.forceSensor_norm(:);
%     forceLabel = 'Force (norm)';
% elseif isfield(ProcData,'forceSensor')
%     force = ProcData.forceSensor(:);
%     forceLabel = 'Force (V)';
% else
%     force = zeros(refLen,1);
%     forceLabel = 'Force';
% end
% force = force(tMask);
% 
% %% Binary force
% if isfield(ProcData,'binForceSensor')
%     binForce = logical(ProcData.binForceSensor(:));
% else
%     binForce = false(refLen,1);
% end
% binForce = binForce(tMask);

%% EMG power
if isfield(ProcData,'EMG') && isfield(ProcData.EMG,'emgPower_norm')
    emgPower = ProcData.EMG.emgPower_norm(:);
    emgLabel = 'EMG power (norm)';
elseif isfield(ProcData,'EMG') && isfield(ProcData.EMG,'emgPower')
    emgPower = ProcData.EMG.emgPower(:);
    emgLabel = 'EMG power';
else
    emgPower = zeros(refLen,1);
    emgLabel = 'EMG power';
end
emgPowerFull = emgPower(tMask);
emgPowerStart = emgPower(tMaskStart);
emgPowerEnd = emgPower(tMaskEnd);

% %% Raw EMG (FIXED)
% rawEMG = [];
% if isfield(ProcData,'EMG') && isfield(ProcData.EMG,'emgSignal')
%     rawEMG = ProcData.EMG.emgSignal(:);
% end
% 
% if isempty(rawEMG)
%     rawEMG = zeros(refLen,1);
% elseif numel(rawEMG) < refLen
%     rawEMG(end+1:refLen) = rawEMG(end);
% elseif numel(rawEMG) > refLen
%     rawEMG = rawEMG(1:refLen);
% end
% rawEMG = rawEMG(tMask);

% %% Load spectrogram
% [folder, base, ~] = fileparts(procDataFileID);
% prefix = regexprep(base,'_ProcData$','');
% 
% specFile = '';
% cand = dir(fullfile(folder, [prefix '*Spec*.mat']));
% if ~isempty(cand)
%     specFile = fullfile(folder, cand(1).name);
% end
% 
% Sspec = []; Fspec = []; Tspec = [];
% isNormSpec = false;
% 
% if ~isempty(specFile) && isfile(specFile)
%     L = load(specFile);
% 
%     if isfield(L,'SpecData') && isfield(L.SpecData,'ECoG')
%         SD = L.SpecData.ECoG;
%         if isfield(SD,'normS')
%             Sspec = SD.normS;
%             Fspec = SD.F;
%             Tspec = SD.T;
%             isNormSpec = true;
%         elseif all(isfield(SD,{'S','F','T'}))
%             Sspec = SD.S;
%             Fspec = SD.F;
%             Tspec = SD.T;
%         end
%     end
% end
% 
% %% --- SPECTROGRAM SAFETY FIXES ---
% 
% if ~isempty(Sspec) && ~isempty(Fspec) && ~isempty(Tspec)
% 
%     % Remove invalid freqs for log scale
%     validIdx = Fspec > 0;
%     Fspec = Fspec(validIdx);
%     Sspec = Sspec(validIdx,:);
% 
%     % Fix NaN/Inf/negative values
%     Sspec(~isfinite(Sspec)) = eps;
%     Sspec(Sspec <= 0) = eps;
% 
%     % Dimension check
%     if size(Sspec,1) ~= numel(Fspec) || size(Sspec,2) ~= numel(Tspec)
%         warning('Spectrogram dimension mismatch. Skipping.');
%         Sspec = [];
%     end
% end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [freq_centroid, meanTemp] = respirationSpectrogramPlot_SF(tdmsPath, segmentStart, segmentEnd)
% BANDPASS_PLOT filters a signal between 0.5 and 3 Hz,
% plots the power spectrum (DFT), and plots the spectrogram
% with frequency center of mass overlay.
%
tdmsDataStruct = openTDMS(tdmsPath);

Fs = str2double(tdmsDataStruct.CameraFrameratePerSecond);
signal = tdmsDataStruct.Digital_Data.Respiration_Sum;
timeVec = (0:length(signal))./Fs;
tMask = timeVec >= segmentStart & timeVec <= segmentEnd;
signal = signal(tMask);

FsTemp = str2double(tdmsDataStruct.TemperatureSamplingRate_Hz_);
signalTemp = tdmsDataStruct.Serial_Data.Temperature;
timeVecTemp = (0:length(signalTemp))./FsTemp;
tMaskTemp = timeVecTemp >= segmentStart & timeVecTemp <= segmentEnd;
meanTemp = mean(signalTemp(tMaskTemp));

% --- Design bandpass filter ---
bpFilt = designfilt('bandpassiir', ...
    'FilterOrder', 4, ...
    'HalfPowerFrequency1', 0.5, ...
    'HalfPowerFrequency2', 3, ...
    'SampleRate', Fs);

% --- Apply zero-phase filter ---
filtered_signal = filtfilt(bpFilt, signal);

% --- Compute DFT ---
N = length(filtered_signal);
X = fft(filtered_signal);
f = (0:N-1)*(Fs/N);
powerX = abs(X).^2 / N;

% figure;
% subplot(4,1,1)
% plot((1:length(signal))./Fs, signal)
% xlabel('Time (s)')
% ylabel('ROI Pixel Sum')
% title('Respiration ROI (Front Chest) Pixel Sum');
% xlim([1 15])
% 
% subplot(4,1,2)
% plot((1:length(filtered_signal))./Fs, filtered_signal)
% xlabel('Time (s)')
% ylabel('Filtered ROI Pixel Sum')
% title('Band-Passed (0.5-3Hz) Respiration ROI (Front Chest) Pixel Sum');
% xlim([1 15])

% % --- Plot Power Spectrum ---
% subplot(4,1,3)
% plot(f, powerX);
% xlim([0 10]);
% xlabel('Frequency (Hz)');
% ylabel('Power');
% title('Power Spectrum of Band-Passed Signal');
% grid on;
% 
% % --- Plot Spectrogram ---
% subplot(4,1,4)

window  = hamming(round(2*Fs));      % 2-second window
overlap = round(0.9 * length(window));
nfft    = 1024;

% Compute spectrogram explicitly (so we can use the data)
[S,F,T] = spectrogram(filtered_signal, window, overlap, nfft, Fs);

% Power spectrogram
P = abs(S).^2;

% % Plot spectrogram
% imagesc(T, F, 10*log10(P));
% axis xy
% colormap jet;
% ax = gca;
% cb = colorbar(ax, 'eastoutside');
% ax.PositionConstraint = 'innerposition';
% ylabel(cb, 'Power')
% title('Spectrogram of Band-Passed Signal');
% xlabel('Time (s)')
% ylabel('Frequency (Hz)')
% ylim([0 5])
% xlim([1 15])

% --- CENTER OF MASS CALCULATION (Frequency Centroid) ---
% Weighted mean frequency at each time slice
freq_centroid = sum(F .* P, 1) ./ sum(P, 1);

% % --- Overlay center of mass ---
% hold on;
% plot(T, freq_centroid, 'w', 'LineWidth', 2);
% plot(T, freq_centroid, 'k--', 'LineWidth', 1); % outline for contrast
% hold off;
end
