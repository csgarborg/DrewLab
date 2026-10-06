function compareREMvsAwakeNoMef(awakeExcelPath, remExcelPath, tdmsTF, saveResultsStructTF, roiPath, saveResultsStructPath, baselineROITF)

%% Load data
if exist('saveResultsStructPath','var') && ~saveResultsStructTF
    load(saveResultsStructPath);
    awake = resultsStruct.awake;
    rem = resultsStruct.rem;
    remWake = resultsStruct.remWake;
else
    if ~exist("tdmsTF","var") || tdmsTF
        awake  = processExcelSegmentsTDMS_SF(awakeExcelPath);
        rem    = processExcelSegmentsTDMS_SF(remExcelPath);
    else
        awake  = processExcelSegmentsROISelect_SF(awakeExcelPath,roiPath,false);
        rem    = processExcelSegmentsROISelect_SF(remExcelPath,roiPath,false);
        remWake = processExcelSegmentsROISelect_SF(remExcelPath,roiPath,true);

        % show ROI locations once
        plotFirstFrameWithROIs(remExcelPath, roiPath);
    end
end
if saveResultsStructTF
    resultsStruct.awake = awake;
    resultsStruct.rem = rem;
    resultsStruct.remWake = remWake;
    save(saveResultsStructPath, 'resultsStruct');
end

if baselineROITF
    % Structs to process
    structNames = {'awake','rem','remWake'};

    for s = 1:length(structNames)

        % Get struct variable
        S = eval(structNames{s});

        % Loop through struct elements
        for i = 1:numel(S)

            fn = fieldnames(S(i));

            for f = 1:length(fn)

                currentField = fn{f};

                % Check if field ends with B
                if endsWith(currentField,'B')

                    % Remove trailing B
                    baseField = extractBefore(currentField, ...
                        strlength(currentField));

                    % Replace original field with B version
                    S(i).(baseField) = S(i).(currentField);

                    % fprintf('Replaced %s with %s in %s(%d)\n', ...
                        % baseField, currentField, ...
                        % structNames{s}, i);
                end
            end
        end

        % Put modified struct back
        eval([structNames{s} ' = S;']);

    end
end

% awake = processMatFiles_SF(awakeExcelPath);
% rem   = processMatFiles_SF(remExcelPath);

nA = length(awake);
nR = length(rem);

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

for r = 1:nA
    plot(1:nSignals, awake(r).explained(1:nSignals),'-o','Color',[0.3 0.5 1])
end
for r = 1:nR
    plot(1:nSignals, rem(r).explained(1:nSignals),'--o','Color',[1 0.4 0.4])
end

legend([hA hR])

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

meanA = mean(A(:,1:nSignals));
meanR = mean(R(:,1:nSignals));

ciA = 1.96*std(A(:,1:nSignals))/sqrt(nA);
ciR = 1.96*std(R(:,1:nSignals))/sqrt(nR);

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

legend({'Awake CI','Awake','REM CI','REM'})

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

pc2A = arrayfun(@(x)x.explained(2),awake);
pc2R = arrayfun(@(x)x.explained(2),rem);

figure

subplot(1,2,1)
bar([mean(pc1A), mean(pc1R)])
hold on
errorbar([1 2], ...
    [mean(pc1A) mean(pc1R)], ...
    [std(pc1A) std(pc1R)], ...
    '.w','LineWidth',1.5)
set(gca,'XTickLabel',{'Awake','REM'})
title('PC1 Strength')
ylabel('% Variance')
ylim([0 100])

subplot(1,2,2)
bar([mean(pc2A),mean(pc2R)])
hold on
errorbar([1 2],[mean(pc2A) mean(pc2R)], ...
    [std(pc2A) std(pc2R)],'.w','LineWidth',1.5)
set(gca,'XTickLabel',{'Awake','REM'})
title('PC2 Strength')
ylabel('% Variance')
ylim([0 100])

%% ==============================
%% 3) Mean correlation (mean ± SD)
%% ==============================

meanCorrA = arrayfun(@(x)mean(x.corrMatrix(~eye(size(x.corrMatrix)))),awake);
meanCorrR = arrayfun(@(x)mean(x.corrMatrix(~eye(size(x.corrMatrix)))),rem);

figure
bar([mean(meanCorrA),mean(meanCorrR)])
hold on
errorbar([1 2],[mean(meanCorrA) mean(meanCorrR)], ...
    [std(meanCorrA) std(meanCorrR)],'.w','LineWidth',1.5)
set(gca,'XTickLabel',{'Awake','REM'})
ylabel('Mean Correlation')
title('Average Correlation ± SD')
ylim([0 1])

%% ==============================
%% 4) Correlation matrices (with labels)
%% ==============================

corrA = cat(3,awake.corrMatrix);
corrR = cat(3,rem.corrMatrix);

meanA = mean(corrA,3);
meanR = mean(corrR,3);

figure
tiledlayout(1,2)

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

sgtitle('Mean Correlation Matrices')

%% Correlation matrix stats test

n1 = length(cat(2,awake.segmentLengths));
n2 = length(cat(2,rem.segmentLengths));

[chi2, p]=JennrichTest_2Matrix(meanA, meanR, n1, n2);

% Format p-value compactly
if p < 0.001
    pstr = 'p < 0.001';
else
    pstr = sprintf('p = %.3f', p);
end

% Create title with chi-squared and p-value (LaTeX interpreter)
tstr = sprintf('$\\chi^2 = %.3f,\\; %s$', chi2, pstr);
% title(tstr, 'Interpreter', 'latex', 'FontSize', 12);

figure
imagesc(meanA-meanR)
title(['Mean Awake - Mean REM Correlation Matricies, ' tstr], 'Interpreter', 'latex')
axis square
clim([-.5 .5])
colorbar
xticks(1:nSignals)
yticks(1:nSignals)
xticklabels(strrep(signalNames,'_',' '))
yticklabels(strrep(signalNames,'_',' '))
xtickangle(45)


%% Correlation distribution comparison

mask = triu(true(nSignals),1);

awakeCorr = meanA(mask);
remCorr   = meanR(mask);

figure
hold on

awakeColor = [0 0.4470 0.7410];
remColor   = [0.8500 0.3250 0.0980];

nBins = 20;   % increase as desired

histogram(awakeCorr,...
    nBins,...
    'Normalization','pdf',...
    'FaceColor',awakeColor,...
    'FaceAlpha',0.4,...
    'EdgeColor','none',...
    'DisplayName','Awake');

histogram(remCorr,...
    nBins,...
    'Normalization','pdf',...
    'FaceColor',remColor,...
    'FaceAlpha',0.4,...
    'EdgeColor','none',...
    'DisplayName','REM');

% Normal fits
pdA = fitdist(awakeCorr,'Normal');
pdR = fitdist(remCorr,'Normal');

x = linspace( ...
    min([awakeCorr; remCorr]), ...
    max([awakeCorr; remCorr]), ...
    1000);

plot(x,pdf(pdA,x),...
    'Color',awakeColor,...
    'LineWidth',3,...
    'DisplayName',sprintf('Awake fit (\\mu=%.2f)',pdA.mu));

plot(x,pdf(pdR,x),...
    'Color',remColor,...
    'LineWidth',3,...
    'DisplayName',sprintf('REM fit (\\mu=%.2f)',pdR.mu));

xlabel('Mean Pairwise Correlation')
ylabel('Probability Density')
title('Distribution of ROI Correlations')

legend('Location','best')
box off

%% ============================================================
%% HIERARCHICAL CLUSTERING OF ROI NETWORKS
%% ============================================================

corrA = cat(3,awake.corrMatrix);
corrR = cat(3,rem.corrMatrix);

if any(corrA(:) < 0) || any(corrR(:) < 0)
    disp('Negative xcorr values detected!')
    negA = sum(corrA(:) < 0);
    totA = numel(corrA);
    pctA = 100 * negA / max(totA,1);

    negR = sum(corrR(:) < 0);
    totR = numel(corrR);
    pctR = 100 * negR / max(totR,1);

    if negA > 0 || negR > 0
        if negA > 0
            fprintf('corrA: %d negative values (%.2f%%)\n', negA, pctA);
        else
            fprintf('corrA: 0 negative values (0.00%%)\n');
        end
        if negR > 0
            fprintf('corrR: %d negative values (%.2f%%)\n', negR, pctR);
        else
            fprintf('corrR: 0 negative values (0.00%%)\n');
        end
    else
        disp('No negative xcorr values detected.')
    end
end

meanCorrA = mean(corrA,3,'omitnan');
meanCorrR = mean(corrR,3,'omitnan');

% Distance matrices
DA = 1 - meanCorrA;
DR = 1 - meanCorrR;

DA(1:size(DA,1)+1:end) = 0;
DR(1:size(DR,1)+1:end) = 0;

% linkage
ZA = linkage(squareform(DA),'average');
ZR = linkage(squareform(DR),'average');

% Optimal ordering
orderA = optimalleaforder(ZA,DA);
orderR = optimalleaforder(ZR,DR);

figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,1)

dendrogram(ZA,0,'Reorder',orderA,'Labels',strrep(signalNames(orderA),'_',' '),'Orientation','top');

title('Awake Hierarchical Clustering','Color','w')

set(gca,...
    'Color','k',...
    'XColor','w',...
    'YColor','w',...
    'LineWidth',1.5)
ax = gca;
ax.XTickLabelRotation = 45;
ylim([0 .8])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,2)

dendrogram(ZR,0,'Reorder',orderR,'Labels',strrep(signalNames(orderR),'_',' '),'Orientation','top');

title('REM Hierarchical Clustering','Color','w')

set(gca,...
    'Color','k',...
    'XColor','w',...
    'YColor','w',...
    'LineWidth',1.5)
ax = gca;
ax.XTickLabelRotation = 45;
ylim([0 .8])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3)

imagesc(meanCorrA(orderA,orderA),[0 1])
axis square
colorbar
xticks(1:nSignals)
yticks(1:nSignals)
xticklabels(strrep(signalNames(orderA),'_',' '))
yticklabels(strrep(signalNames(orderA),'_',' '))
xtickangle(45)

title('Awake Reordered Correlation Matrix','Color','w')

set(gca,...
    'Color','k',...
    'XColor','w',...
    'YColor','w')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,4)

imagesc(meanCorrR(orderR,orderR),[0 1])
axis square
colorbar
xticks(1:nSignals)
yticks(1:nSignals)
xticklabels(strrep(signalNames(orderR),'_',' '))
yticklabels(strrep(signalNames(orderR),'_',' '))
xtickangle(45)

title('REM Reordered Correlation Matrix','Color','w')

set(gca,...
    'Color','k',...
    'XColor','w',...
    'YColor','w')

%% Cophenetic correlation

cophenA = cophenet(ZA,squareform(DA));
cophenR = cophenet(ZR,squareform(DR));

figure('Color','k')

bar([1 2],[cophenA cophenR],0.6)

set(gca,...
    'Color','k',...
    'XTick',[1 2],...
    'XTickLabel',{'Awake','REM'},...
    'XColor','w',...
    'YColor','w')

ylabel('Cophenetic Correlation')

title('Hierarchical Structure','Color','w')

% %% clustered correlation matrices
% 
% figure
% 
nRec = numel(awake);
nRecREM = numel(rem);
% 
% for r=1:nRec
% 
%     subplot(2,nRec,r)
% 
%     C = awake(r).corrMatrix;
% 
%     D = 1-C;
%     D(1:size(D,1)+1:end)=0;
% 
%     Z = linkage(squareform(D),'average');
%     ord = optimalleaforder(Z,D);
% 
%     imagesc(C(ord,ord),[-1 1])
%     axis square
% 
%     title(sprintf('Awake %d',r),'Color','w')
% 
%     set(gca,'Color','k','XColor','w','YColor','w')
% 
%     subplot(2,nRec,nRec+r)
% 
%     C = rem(r).corrMatrix;
% 
%     D = 1-C;
%     D(1:size(D,1)+1:end)=0;
% 
%     Z = linkage(squareform(D),'average');
%     ord = optimalleaforder(Z,D);
% 
%     imagesc(C(ord,ord),[-1 1])
%     axis square
% 
%     title(sprintf('REM %d',r),'Color','w')
% 
%     set(gca,'Color','k','XColor','w','YColor','w')
% 
% end
% 
% colormap(parula)

%% number of large clusters

clusterCutoff = 0.55;
clusterAwake = zeros(nRec,1);
clusterREM = zeros(nRecREM,1);

for r=1:nRec

    D = 1-awake(r).corrMatrix;
    D(1:size(D,1)+1:end)=0;

    Z = linkage(squareform(D),'average');

    idx = cluster(Z,'cutoff',clusterCutoff,'Criterion','distance');

    clusterAwake(r)=numel(unique(idx));
end
for r = 1:nRecREM

    D = 1-rem(r).corrMatrix;
    D(1:size(D,1)+1:end)=0;

    Z = linkage(squareform(D),'average');

    idx = cluster(Z,'cutoff',clusterCutoff,'Criterion','distance');

    clusterREM(r)=numel(unique(idx));

end
figure

subplot(121)

bar([mean(clusterAwake) mean(clusterREM)])

hold on

errorbar([1 2],...
    [mean(clusterAwake) mean(clusterREM)],...
    [std(clusterAwake)/sqrt(nRec) std(clusterREM)/sqrt(nRecREM)],...
    'w','LineStyle','none')

set(gca,...
    'Color','k',...
    'XTick',[1 2],...
    'XTickLabel',{'Awake','REM'},...
    'XColor','w',...
    'YColor','w')

ylabel(['Number of Clusters, Cluster cutoff = ' num2str(clusterCutoff)])

subplot(122)

plotSpread({clusterAwake(:),clusterREM(:)});

set(gca,...
    'XTick',1:2,...
    'XTickLabel',{'Awake','REM'},...
    'Color','k',...
    'XColor','w',...
    'YColor','w');

ylabel(['Number of Clusters, Cluster cutoff = ' num2str(clusterCutoff)])

%% Number of clusters as a function of cutoff (0 to 1, step 0.1)
cutoffRange = 0.01:0.01:1;
nCutoffs = numel(cutoffRange);

clusterAwakeByCutoff = zeros(nRec,nCutoffs);
clusterREMByCutoff = zeros(nRecREM,nCutoffs);

for c = 1:nCutoffs
    cutoff = cutoffRange(c);
    for r = 1:nRec
        D = 1-awake(r).corrMatrix;
        D(1:size(D,1)+1:end)=0;
        Z = linkage(squareform(D),'average');
        idx = cluster(Z,'cutoff',cutoff,'Criterion','distance');
        clusterAwakeByCutoff(r,c) = numel(unique(idx));
    end
    for r = 1:nRecREM
        D = 1-rem(r).corrMatrix;
        D(1:size(D,1)+1:end)=0;
        Z = linkage(squareform(D),'average');
        idx = cluster(Z,'cutoff',cutoff,'Criterion','distance');
        clusterREMByCutoff(r,c) = numel(unique(idx));
    end
end

figure
hold on
errorbar(cutoffRange, mean(clusterAwakeByCutoff,1), ...
    std(clusterAwakeByCutoff,[],1)/sqrt(nRec), ...
    '-o','Color',awakeColor,'MarkerFaceColor',awakeColor,'LineWidth',1.5)
errorbar(cutoffRange, mean(clusterREMByCutoff,1), ...
    std(clusterREMByCutoff,[],1)/sqrt(nRecREM), ...
    '-o','Color',remColor,'MarkerFaceColor',remColor,'LineWidth',1.5)

set(gca,'Color','k','XColor','w','YColor','w')
xlabel('Cluster cutoff')
ylabel('Number of Clusters')
legend({'Awake','REM'},'TextColor','w','Location','best')
title('Cluster Count Sensitivity to Cutoff: Awake vs. REM','Color','w')
xlim([0 1])

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

for r=1:nA, plot(awake(r).coeff(:,1),'-o','Color',[0.3 0.5 1]); end
for r=1:nR, plot(rem(r).coeff(:,1),'--o','Color',[1 0.4 0.4]); end

title('Individual')
xticks(1:nSignals)
xticklabels(strrep(signalNames,'_',' '))
xtickangle(45)

legend([hA hR])
ylim([-1 1])

% mean + CI
nexttile
hold on

A = reshape([awake.coeff],nSignals,[],nA);
R = reshape([rem.coeff],nSignals,[],nR);

pc1A = squeeze(A(:,1,:))';
pc1R = squeeze(R(:,1,:))';

meanA = mean(pc1A);
meanR = mean(pc1R);

ciA = 1.96*std(pc1A)/sqrt(nA);
ciR = 1.96*std(pc1R)/sqrt(nR);

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

legend({'Awake CI','Awake','REM CI','REM'})

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

for r=1:nA, plot(awake(r).coeff(:,2),'-o','Color',[0.3 0.5 1]); end
for r=1:nR, plot(rem(r).coeff(:,2),'--o','Color',[1 0.4 0.4]); end

title('Individual')
legend([hA hR])
ylim([-1 1])

xticks(1:nSignals)
xticklabels(strrep(signalNames,'_',' '))
xtickangle(45)

nexttile
hold on

pc2A = squeeze(A(:,2,:))';
pc2R = squeeze(R(:,2,:))';

meanA = mean(pc2A);
meanR = mean(pc2R);

ciA = 1.96*std(pc2A)/sqrt(nA);
ciR = 1.96*std(pc2R)/sqrt(nR);

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

legend({'Awake CI','Awake','REM CI','REM'})

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

GL_A = zeros(nA,1);
GL_R = zeros(nR,1);

for r = 1:nA
    e = awake(r).explained;
    GL_A(r) = e(1) / mean(e(2:min(5,end)));
end

for r = 1:nR
    e = rem(r).explained;
    GL_R(r) = e(1) / mean(e(2:min(5,end)));
end

figure
bar([mean(GL_A), mean(GL_R)])
set(gca,'XTickLabel',{'Awake','REM'})
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

title('PC1 Power Spectrum')
xlabel('Frequency')
ylabel('Power')
grid on

% --- Dominant frequency ---
nexttile
domA = zeros(nA,1);
domR = zeros(nR,1);

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

bar([mean(domA), mean(domR)])
hold on
errorbar([1 2], ...
    [mean(domA), mean(domR)], ...
    [std(domA), std(domR)], ...
    '.w','LineWidth',1.5)

set(gca,'XTickLabel',{'Awake','REM'})
ylabel('Dominant Frequency')
title('PC1 Dominant Frequency')

% --- Spectral entropy ---
nexttile
entA = zeros(nA,1);
entR = zeros(nR,1);

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

bar([mean(entA), mean(entR)])
hold on
errorbar([1 2], ...
    [mean(entA), mean(entR)], ...
    [std(entA), std(entR)], ...
    '.w','LineWidth',1.5)

set(gca,'XTickLabel',{'Awake','REM'})
ylabel('Spectral Entropy')
title('PC1 Spectral Entropy')

% --- High / Low frequency ratio ---
nexttile
hlA = zeros(nA,1);
hlR = zeros(nR,1);

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

bar([mean(hlA), mean(hlR)])
hold on
errorbar([1 2], ...
    [mean(hlA), mean(hlR)], ...
    [std(hlA), std(hlR)], ...
    '.w','LineWidth',1.5)

set(gca,'XTickLabel',{'Awake','REM'})
ylabel('High / Low Ratio')
title('High vs Low Frequency Power')

sgtitle('Frequency Analysis')

%% ==============================
%% 10) Participation ratio + effective dimensionality
%% ==============================

PR_A = zeros(nA,1);
PR_R = zeros(nR,1);

Dim90_A = zeros(nA,1);
Dim90_R = zeros(nR,1);

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

figure

subplot(1,2,1)
bar([mean(PR_A), mean(PR_R)])
hold on
errorbar([1 2], ...
    [mean(PR_A), mean(PR_R)], ...
    [std(PR_A), std(PR_R)], ...
    '.w','LineWidth',1.5)

set(gca,'XTickLabel',{'Awake','REM'})
ylabel('Participation Ratio')
title('Effective Dimensionality')

subplot(1,2,2)
bar([mean(Dim90_A), mean(Dim90_R)])
hold on
errorbar([1 2], ...
    [mean(Dim90_A), mean(Dim90_R)], ...
    [std(Dim90_A), std(Dim90_R)], ...
    '.w','LineWidth',1.5)

set(gca,'XTickLabel',{'Awake','REM'})
ylabel('# PCs for 90% Variance')
title('Dimensionality to 90% Variance')

%% ==============================
%% 11) Variability / Fano factor
%% ==============================

F_A = zeros(nA,1);
F_R = zeros(nR,1);

for r = 1:nA
    x = abs(awake(r).score(:,1));
    F_A(r) = var(x) / (mean(x) + eps);
end

for r = 1:nR
    x = abs(rem(r).score(:,1));
    F_R(r) = var(x) / (mean(x) + eps);
end

figure
bar([mean(F_A), mean(F_R)])
hold on
errorbar([1 2], ...
    [mean(F_A), mean(F_R)], ...
    [std(F_A), std(F_R)], ...
    '.w','LineWidth',1.5)

set(gca,'XTickLabel',{'Awake','REM'})
ylabel('Fano Factor')
title('Burstiness / Variability')

%% ==============================
%% 12) Simple twitch event detection
%% ==============================

thrMult = 2; % threshold = mean + 2*std

rateA = zeros(nA,1);
rateR = zeros(nR,1);

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

figure
bar([mean(rateA), mean(rateR)])
hold on
errorbar([1 2], ...
    [mean(rateA), mean(rateR)], ...
    [std(rateA), std(rateR)], ...
    '.w','LineWidth',1.5)

set(gca,'XTickLabel',{'Awake','REM'})
ylabel('Event Count')
title('Twitch Event Rate')

%% ==============================
%% 13) REM Event Lengths
%% Every segment as its own point
%% Using plotSpread + mean ± SD
%% ==============================

% Assumes:
% rem(r).segmentLengths
%
% Each individual segment length will be plotted
% as its own point (NO per-recording averaging)

allLenR = [];

for r = 1:nR
    if isfield(rem(r),'segmentLengths') && ~isempty(rem(r).segmentLengths)
        allLenR = [allLenR; rem(r).segmentLengths(:)];
    end
end

% Remove NaNs just in case
allLenR = allLenR(~isnan(allLenR));

figure
hold on

plotSpread({allLenR}, ...
    'xNames', {'REM'}, ...
    'showMM', 5, ... % mean ± SD
    'distributionColors', {[1 0.4 0.4]}, ...
    'distributionMarkers', {'o'}, ...
    'spreadWidth', 0.6)

ylabel('Event Length (s)')
title('REM Event Lengths')
grid on

%% ==============================
%% 14) Event "Energy" + First-Responding ROI
%% ==============================
%
% Energy = area under the z-scored ROI curve for one event
%          (trapz over time, units = z-score * seconds)
%
% For every Awake / REM event (pooled across recordings):
%   - mean energy across all ROIs
%   - energy of the single biggest-energy ROI (and which ROI that was)
%   - which ROI first shows "significant" motion (z > energyThresh
%     sustained for energyMinDur seconds) and its onset latency
%
% Z-scoring is controlled by energyZMode:
%   'perEvent'     - every event is z-scored on its own (its own mean/SD)
%   'concatenated' - z-score computed over all events of a recording
%                    (the same z-score used for PCA etc.)
% Uses results(r).eventZSingle / eventZConcat / eventTime, which are saved
% by processExcelSegmentsROISelect_SF (re-run processing once to create them).

energyZMode        = 'perEvent';  % 'perEvent' or 'concatenated' (see above)
energyThresh       = 2;      % z-score threshold for "significant" motion
energyMinDur       = 0.1;    % s the ROI must stay above threshold
energyPositiveOnly = false;  % true -> only count area above z = 0
energyPerSecond    = false;  % true -> divide energy by event duration (z)

if strcmp(energyZMode,'perEvent')
    % An event z-scored on its own has mean 0, so its signed area is ~0 by
    % construction. Energy is therefore the area ABOVE z = 0 in this mode.
    energyPositiveOnly = true;
    zModeText = 'per event';
else
    zModeText = 'across concatenated events';
end

% Reference ROI order (every recording must match this exactly; recordings
% with the same ROIs in a different order are reordered, anything else is
% skipped with a warning)
refROINames = awake(1).signalNames;

evtA = computeEventEnergyAndOnset_SF(awake, refROINames, energyZMode, energyThresh, energyMinDur, energyPositiveOnly, energyPerSecond);
evtR = computeEventEnergyAndOnset_SF(rem,   refROINames, energyZMode, energyThresh, energyMinDur, energyPositiveOnly, energyPerSecond);

eventEnergy.awake = evtA;
eventEnergy.rem   = evtR;

evtROINames = strrep(refROINames,'_',' ');
nEvtROI     = numel(evtROINames);
nEvtA       = numel(evtA.meanEnergy);
nEvtR       = numel(evtR.meanEnergy);

cEvtA = [0.3 0.5 1];
cEvtR = [1 0.4 0.4];

if energyPerSecond
    energyLabel = 'Energy (z-score)';
else
    energyLabel = 'Energy (z-score \cdot s)';
end
if energyPositiveOnly
    energyLabel = strrep(energyLabel, 'Energy', 'Positive energy');
end

% ---- Figure 1: mean energy + max-ROI energy per event ----
figure
tiledlayout(1,2,'TileSpacing','compact')

nexttile
hold on
plotSpread({evtA.meanEnergy, evtR.meanEnergy}, ...
    'xNames', {'Awake','REM'}, ...
    'showMM', 5, ...
    'distributionColors', {cEvtA, cEvtR}, ...
    'distributionMarkers', {'o','o'}, ...
    'spreadWidth', 0.6)
ylabel(energyLabel)
title('Mean Energy Across ROIs')
grid on

nexttile
hold on
plotSpread({evtA.maxEnergy, evtR.maxEnergy}, ...
    'xNames', {'Awake','REM'}, ...
    'showMM', 5, ...
    'distributionColors', {cEvtA, cEvtR}, ...
    'distributionMarkers', {'o','o'}, ...
    'spreadWidth', 0.6)
ylabel(energyLabel)
title('Biggest-Energy ROI')
grid on

sgtitle(sprintf('Event Energy (area under ROI curve, z-scored %s)', zModeText))

% ---- Figure 2: which ROI is biggest / which ROI is first ----
pctMaxA = 100 * histcounts(evtA.maxROI, 0.5:1:nEvtROI+0.5) / max(nEvtA,1);
pctMaxR = 100 * histcounts(evtR.maxROI, 0.5:1:nEvtROI+0.5) / max(nEvtR,1);

nOnsetA = sum(any(~isnan(evtA.onsetLatency),2));
nOnsetR = sum(any(~isnan(evtR.onsetLatency),2));
pctFirstA = 100 * sum(evtA.firstWeight,1) / max(nOnsetA,1);
pctFirstR = 100 * sum(evtR.firstWeight,1) / max(nOnsetR,1);

figure
tiledlayout(2,1,'TileSpacing','compact')

nexttile
hold on
bMax = bar(1:nEvtROI, [pctMaxA(:) pctMaxR(:)], 'grouped');
bMax(1).FaceColor = cEvtA;
bMax(2).FaceColor = cEvtR;
xticks(1:nEvtROI)
xticklabels(evtROINames)
xtickangle(45)
ylabel('% of Events')
title('ROI With Biggest Energy')
legend({sprintf('Awake (n=%d)',nEvtA), sprintf('REM (n=%d)',nEvtR)})
grid on

nexttile
hold on
bFirst = bar(1:nEvtROI, [pctFirstA(:) pctFirstR(:)], 'grouped');
bFirst(1).FaceColor = cEvtA;
bFirst(2).FaceColor = cEvtR;
xticks(1:nEvtROI)
xticklabels(evtROINames)
xtickangle(45)
ylabel('% of Events')
title(sprintf('First ROI With Significant Motion (z > %g for %g s)', energyThresh, energyMinDur))
legend({sprintf('Awake (n=%d)',nOnsetA), sprintf('REM (n=%d)',nOnsetR)})
grid on

% ---- Figure 3: onset latency of every ROI in every event ----
% Black = ROI never reached threshold during that event
figure
tiledlayout(1,2,'TileSpacing','compact')

nexttile
imagesc(evtA.onsetLatency, 'AlphaData', ~isnan(evtA.onsetLatency))
set(gca,'Color','k')
xticks(1:nEvtROI)
xticklabels(evtROINames)
xtickangle(45)
xlabel('ROI')
ylabel('Awake Event #')
title('Awake')
cb = colorbar;
cb.Label.String = 'Onset Latency (s)';

nexttile
imagesc(evtR.onsetLatency, 'AlphaData', ~isnan(evtR.onsetLatency))
set(gca,'Color','k')
xticks(1:nEvtROI)
xticklabels(evtROINames)
xtickangle(45)
xlabel('ROI')
ylabel('REM Event #')
title('REM')
cb = colorbar;
cb.Label.String = 'Onset Latency (s)';

sgtitle('Time From Event Start Until Each ROI Shows Significant Motion')


%% ---- Statistics: Awake vs REM ----
% Events are the unit of analysis, but events from one recording are not
% independent, so each measure is tested three ways:
%   1) Mann-Whitney on all events        (ranksum)
%   2) Mann-Whitney on per-recording means
%   3) Linear mixed model, random intercept per recording
% Effect size = Cliff's delta (positive = Awake larger, range -1..1)

rng(0)            % reproducible permutation tests
nPermEvt = 10000;

nRecA = numel(unique(evtA.recIdx));
nRecR = numel(unique(evtR.recIdx));

% fprintf('\n===== Event energy statistics: Awake vs REM =====\n');
% fprintf('Awake: %d events from %d recordings | REM: %d events from %d recordings\n', ...
%     nEvtA, nRecA, nEvtR, nRecR);
% fprintf('Z-scoring: %s | onset threshold: z > %g for %g s\n\n', zModeText, energyThresh, energyMinDur);

evtMeasures = {
    'Mean energy across ROIs',  evtA.meanEnergy,                       evtR.meanEnergy
    'Biggest-ROI energy',       evtA.maxEnergy,                        evtR.maxEnergy
    'Event duration (s)',       evtA.duration,                         evtR.duration
    'ROIs recruited per event', sum(~isnan(evtA.onsetLatency),2),      sum(~isnan(evtR.onsetLatency),2)};

% fprintf('%-26s %9s %9s | %10s %10s %10s | %8s\n', ...
%     'Measure','Awake med','REM med','p events','p per-rec','p mixed','Cliff d');
for mm = 1:size(evtMeasures,1)
    vA = evtMeasures{mm,2};
    vR = evtMeasures{mm,3};

    pEvt = ranksum(vA, vR);
    pRec = ranksum(splitapply(@mean, vA, findgroups(evtA.recIdx)), ...
                   splitapply(@mean, vR, findgroups(evtR.recIdx)));
    pLME = lmeStatePValue_SF(vA, evtA.recIdx, vR, evtR.recIdx);
    dEff = cliffsDelta_SF(vA, vR);

    % fprintf('%-26s %9.3g %9.3g | %10.3g %10.3g %10.3g | %8.2f\n', ...
    %     evtMeasures{mm,1}, median(vA), median(vR), pEvt, pRec, pLME, dEff);

    eventEnergy.stats.(sprintf('measure%d',mm)) = struct('name',evtMeasures{mm,1}, ...
        'medianAwake',median(vA),'medianREM',median(vR), ...
        'pEvents',pEvt,'pPerRecording',pRec,'pMixedModel',pLME,'cliffsDelta',dEff);
end

% Fraction of events in which ROIs reached threshold at all
nYesA = sum(any(~isnan(evtA.onsetLatency),2));
nYesR = sum(any(~isnan(evtR.onsetLatency),2));
[~, pAny] = fishertest([nYesA nEvtA-nYesA; nYesR nEvtR-nYesR]);
% fprintf('\nEvents with >=1 ROI above threshold: Awake %d/%d (%.0f%%), REM %d/%d (%.0f%%), Fisher p = %.3g\n', ...
%     nYesA, nEvtA, 100*nYesA/nEvtA, nYesR, nEvtR, 100*nYesR/nEvtR, pAny);

% Does the identity of the biggest / first ROI differ between states?
% (permutation chi-square on event labels; handles small counts)
oneHotA = accumarray([(1:nEvtA)' evtA.maxROI], 1, [nEvtA nEvtROI]);
oneHotR = accumarray([(1:nEvtR)' evtR.maxROI], 1, [nEvtR nEvtROI]);
[chiMax, pMaxROI]     = permChi2Counts_SF(oneHotA, oneHotR, nPermEvt);
[chiFirst, pFirstROI] = permChi2Counts_SF(evtA.firstWeight, evtR.firstWeight, nPermEvt);
fprintf('Biggest-energy ROI distribution differs Awake vs REM: chi2 = %.1f, permutation p = %.3g\n', chiMax, pMaxROI);
% fprintf('First-responding ROI distribution differs Awake vs REM: chi2 = %.1f, permutation p = %.3g\n', chiFirst, pFirstROI);

% Energy depends on how long the event is
for gg = 1:2
    if gg == 1, ev = evtA; nm = 'Awake'; else, ev = evtR; nm = 'REM'; end
    if numel(ev.duration) >= 3
        [rho, pRho] = corr(ev.duration, ev.meanEnergy, 'Type','Spearman');
        % fprintf('%s: Spearman duration vs mean energy rho = %.2f (p = %.3g)\n', nm, rho, pRho);
    end
end

% Per-ROI energy, Awake vs REM (Benjamini-Hochberg FDR across ROIs)
pROI    = nan(nEvtROI,1);
dROI    = nan(nEvtROI,1);
for j = 1:nEvtROI
    pROI(j) = ranksum(evtA.energy(:,j), evtR.energy(:,j));
    dROI(j) = cliffsDelta_SF(evtA.energy(:,j), evtR.energy(:,j));
end
qROI = bhFDR_SF(pROI);

% fprintf('\nPer-ROI energy (BH-FDR across %d ROIs)\n', nEvtROI);
% fprintf('%-14s %9s %9s %10s %10s %8s\n','ROI','Awake med','REM med','p','q (FDR)','Cliff d');
% for j = 1:nEvtROI
%     fprintf('%-14s %9.3g %9.3g %10.3g %10.3g %8.2f\n', evtROINames{j}, ...
%         median(evtA.energy(:,j)), median(evtR.energy(:,j)), pROI(j), qROI(j), dROI(j));
% end
eventEnergy.stats.perROI = table(evtROINames(:), pROI, qROI, dROI, ...
    'VariableNames', {'ROI','p','qFDR','cliffsDelta'});
eventEnergy.stats.pAnyOnset  = pAny;
eventEnergy.stats.pMaxROI    = pMaxROI;
eventEnergy.stats.pFirstROI  = pFirstROI;

% ---- Figure 4: per-ROI energy (mean +/- SEM), * = FDR q < 0.05 ----
mA   = mean(evtA.energy,1);
mR   = mean(evtR.energy,1);
seA  = std(evtA.energy,0,1) / sqrt(nEvtA);
seR  = std(evtR.energy,0,1) / sqrt(nEvtR);

figure
hold on
bROI = bar(1:nEvtROI, [mA(:) mR(:)], 'grouped');
bROI(1).FaceColor = cEvtA;
bROI(2).FaceColor = cEvtR;
% grouped-bar x offsets for the error bars
xA = bROI(1).XEndPoints;
xR = bROI(2).XEndPoints;
errorbar(xA, mA, seA, '.', 'Color', [0.8 0.8 0.8], 'LineWidth', 1.2, 'HandleVisibility','off')
errorbar(xR, mR, seR, '.', 'Color', [0.8 0.8 0.8], 'LineWidth', 1.2, 'HandleVisibility','off')

yTop = max([mA+seA; mR+seR], [], 1);
yPad = 0.05 * (max(yTop) - min([0 min(mA-seA) min(mR-seR)]));
for j = 1:nEvtROI
    if qROI(j) < 0.05
        text(j, yTop(j)+yPad, '*', 'HorizontalAlignment','center', 'FontSize',16)
    end
end

xticks(1:nEvtROI)
xticklabels(evtROINames)
xtickangle(45)
ylabel(energyLabel)
title('Per-ROI Event Energy (mean \pm SEM; * = FDR q < 0.05)')
legend({sprintf('Awake (n=%d)',nEvtA), sprintf('REM (n=%d)',nEvtR)})
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
%% Awake vs REM
%% ==============================

roiLabels = signalNames;
nROI = length(roiLabels);

cA = [0.3 0.5 1];
cR = [1 0.4 0.4];

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
    
    
    valsA = valsA(~isnan(valsA));
    valsR = valsR(~isnan(valsR));
    
    plotSpread({valsA, valsR}, ...
        'xNames', {'Awake','REM'}, ...
        'showMM', 5, ...
        'distributionColors', {cA, cR}, ...
        'distributionMarkers', {'o','o'}, ...
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
    
    
    valsA = valsA(~isnan(valsA));
    valsR = valsR(~isnan(valsR));
    
    plotSpread({valsA, valsR}, ...
        'xNames', {'Awake','REM'}, ...
        'showMM', 5, ...
        'distributionColors', {cA, cR}, ...
        'distributionMarkers', {'o','o'}, ...
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
    
    
    valsA = valsA(~isnan(valsA));
    valsR = valsR(~isnan(valsR));
    
    plotSpread({valsA, valsR}, ...
        'xNames', {'Awake','REM'}, ...
        'showMM', 5, ...
        'distributionColors', {cA, cR}, ...
        'distributionMarkers', {'o','o'}, ...
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
tiledlayout(1,2,'TileSpacing','compact')

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


sgtitle('Mean ROI Motion Across Events')

%% ==============================
%% Mean ROI Motion Across Event
%% Binned at 1%, 5%, and 10%
%%
%% Assumes:
%% awake(r).meanROIMotion{k}
%% rem(r).meanROIMotion{k}
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
    tiledlayout(1,2,'TileSpacing','compact')

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
tiledlayout(1,2,'TileSpacing','compact')

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

sgtitle('EMG Power Across Events (Raw)')


%% =========================================
%% FIGURES 2–4 : BINNED (1%, 5%, 10%)
%% =========================================

binSizes = [1 5 10];

for b = 1:length(binSizes)

    binSize = binSizes(b);

    figure
    tiledlayout(1,2,'TileSpacing','compact')

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

    sgtitle(sprintf('EMG Power Across Events (%d%% Bins)', binSize))

end

%% ==============================
%% EMG Power START / END Across Events
%%
%% Layout:
%% [ Awake      REM
%%   AwakeMean  REMMean ]
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

alphaRaw = 0.25;
alphaCI  = 0.18;

for v = 1:length(varsToPlot)

    varName   = varsToPlot{v};
    mainTitle = titlesToPlot{v};

    figure
    tiledlayout(2,2,'TileSpacing','compact')

    groupStructs = {awake, rem};
    groupNs      = [nA, nR];
    groupTitles  = {'Awake','REM'};
    groupColors  = {cA, cR};

    for g = 1:2

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

        nexttile(g+2)
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
%% PERCENT OF EMG POWER ABOVE THRESHOLD
%%
%% Threshold definition:
%% threshold = mean(lowest 10% of EMG values) + 0.5
%%
%% Calculates:
%% percent of EMG points above threshold
%%
%% Uses:
%% rem(r).emgPower
%%
%% Each emgPower cell:
%% col 1 = time
%% col 2 = EMG power
%%
%% Output:
%% plotSpread plot for REM
%% each event plotted as its own point
%% ============================================================

groups = {rem};
groupNames = {'REM'};

allPercents = cell(1,length(groups));

for g = 1:length(groups)

    dataStruct = groups{g};

    percentVals = [];

    for r = 1:length(dataStruct)

        if isfield(dataStruct(r),'emgPower') && ...
                ~isempty(dataStruct(r).emgPower)

            for k = 1:length(dataStruct(r).emgPower)

                thisData = dataStruct(r).emgPower{k};

                if isempty(thisData) || size(thisData,2) < 2
                    continue
                end

                emg = thisData(:,2);

                emg = emg(~isnan(emg));

                if isempty(emg)
                    continue
                end

                %% --------------------------------
                %% Calculate threshold
                %% --------------------------------

                sortedVals = sort(emg);

                nLow = max(round(length(sortedVals)*0.10),1);

                lowVals = sortedVals(1:nLow);

                threshold = mean(lowVals) + 1;

                %% --------------------------------
                %% Percent above threshold
                %% --------------------------------

                percentAbove = ...
                    sum(emg > threshold) / length(emg) * 100;

                percentVals(end+1,1) = percentAbove; %#ok<AGROW>

            end
        end
    end

    allPercents{g} = percentVals;

end


%% ============================================================
%% PLOT
%% ============================================================

figure
hold on

plotSpread(allPercents, ...
    'xNames',groupNames, ...
    'showMM',5, ...
    'distributionMarkers','o')

ylabel('% EMG Power Above Threshold')
title('Percent of EMG Activity Above Baseline Threshold')

grid on
set(gca,'FontSize',12)

%% Optional:
%% make mean/error bars white

h = findobj(gca,'Color','r');
set(h,'Color','w')


%% ============================================================
%% REM STOP ("Wake Events")
%%
%% Creates versions of:
%% 1) Percent time moving
%% 2) Mean ROI motion across events
%% 3) EMG Power across events
%%
%% Uses:
%% remStop
%%
%% Same formatting as before:
%% - no Awake plots
%% - only REM Stop
%% - add "Wake Events" to all titles
%% ============================================================

nR = length(remWake);

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
                    allR = [allR; thisVal(:,roi)];
                end
            end
        end

        plotSpread({allR}, ...
            'xNames', {'REM Stop'}, ...
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
tiledlayout(2,1,'TileSpacing','compact')

groups = {remWake};
groupNames = {'REM Wake'};
groupColors = { ...
    [1 0.4 0.4]};

for g = 1:length(groups)

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
    nexttile(g+1)
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
%% [REM Wake traces
%%  REM Wake mean+CI]
%% ============================================================

binSizes = [1 5 10];

for b = 1:length(binSizes)

    binSize = binSizes(b);

    figure
    tiledlayout(2,1,'TileSpacing','compact')

    for g = 1:length(groups)

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

        nexttile(g+1)
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
%% [REM Wake traces
%%  REM Wake mean+CI]
%% ============================================================

%% -------------------------------
%% MAIN FIGURE (raw + mean CI)
%% -------------------------------

figure
tiledlayout(2,1,'TileSpacing','compact')

for g = 1:length(groups)

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
    nexttile(g+1)
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
    tiledlayout(2,1,'TileSpacing','compact')

    for g = 1:length(groups)

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

        nexttile(g+1)
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
%% Awake / REM
%%
%% Each cell contains:
%% column 1 = % event progression (0–100)
%% column 2 = respiration frequency centroid
%% ============================================================

groups = {awake, rem};
groupNames = {'Awake','REM'};
groupColors = { ...
    [0.3 0.5 1], ...
    [1 0.4 0.4]};

%% ============================================================
%% FIGURE 1 — FULL EVENT TRACES + MEAN ± 95% CI
%% ============================================================

figure
tiledlayout(2,2,'TileSpacing','compact')

for g = 1:2

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

sgtitle('Respiration Rate Across Events')


%% ============================================================
%% FIGURES 2–4 — BINNED RESPIRATION RATE
%% (updated: traces on top, mean ± 95% CI below)
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

    sgtitle(sprintf('Respiration Rate — %d%% Binning', binSize))

end

%% ============================================================
%% WAKE EVENTS — RESPIRATION RATE (respFreqCentroid)
%%
%% Uses:
%% remStop(r).respFreqCentroid
%%
%% Same format as before:
%% Figure 1 -> full event traces (0–100%)
%% Figure 2 -> 1% bins
%% Figure 3 -> 5% bins
%% Figure 4 -> 10% bins
%%
%% Includes:
%% REM Stop
%%
%% Titles include "Wake Events"
%% ============================================================

groups = {remWake};
groupNames = {'REM Wake'};
groupColors = { ...
    [1 0.4 0.4]};

%% ============================================================
%% FIGURE 1 — FULL EVENT TRACES + MEAN ± 95% CI
%% ============================================================

figure
tiledlayout(2,1,'TileSpacing','compact')

for g = 1:length(groups)

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

    nexttile(g+1)
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
    tiledlayout(2,1,'TileSpacing','compact')

    for g = 1:length(groups)

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

        nexttile(g+1)
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
%% plotSpread comparison of Awake / REM
%%
%% Bottom:
%% Histogram of REM event temperatures
%% ============================================================

figure
tiledlayout(2,1,'TileSpacing','compact')

groups = {awake, rem};
groupNames = {'Awake','REM'};

plotData = cell(1,2);

for g = 1:2

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
    'distributionColors', {'b','r'});

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

remTemps = plotData{2}; % REM group

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
%% Awake / REM
%%
%% Every left-right pair from every recording
%% becomes its own point on the plot
%%
%% Higher correlation = more symmetric movement
%% Lower correlation = more asymmetric movement
%% ============================================================

figure
hold on

groups = {awake, rem};
groupNames = {'Awake','REM'};

plotData = cell(1,2);

for g = 1:2

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
    'distributionColors', {'b','r'});

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

groups = {awake, rem};
groupNames = {'Awake','REM'};

peakCorrData = cell(1,2);
peakLagData  = cell(1,2);
twitchCoincidenceData = cell(1,2);

maxLag = 50;          % frames
threshold = 2;        % active movement threshold

for g = 1:2

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
%% Awake / REM
%%
%% Smaller values = more synchronous
%% Larger values = more asymmetric
%% ============================================================

threshold = 2.5;   % z-score threshold for twitch detection
minGap    = 5;     % minimum separation between twitches (frames)

groups = {awake, rem};
groupNames = {'Awake','REM'};

allDelays = cell(1,2);

for g = 1:2

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
function evt = computeEventEnergyAndOnset_SF(S, refNames, zMode, thr, minDur, posOnly, perSec)
% Per-event energy (area under the z-scored ROI curve) and ROI onset
% latencies, pooled across all recordings in S. Every event is processed on
% its own. ROI columns are forced into the order given by refNames.
%
% zMode: 'perEvent'     - each event z-scored on its own
%        'concatenated' - z-score taken over all events of the recording
%
% Output (one row per event):
%   energy         [nEvents x nROI]  area under z-score curve (trapz, z*s)
%   meanEnergy     [nEvents x 1]     mean energy across ROIs
%   maxEnergy      [nEvents x 1]     energy of the biggest-energy ROI
%   maxROI         [nEvents x 1]     index of the biggest-energy ROI
%   onsetLatency   [nEvents x nROI]  s from event start until z > thr is
%                                    sustained for minDur (NaN = never)
%   firstWeight    [nEvents x nROI]  1 for the first ROI to respond
%                                    (split evenly if tied; 0 if none)
%   duration       [nEvents x 1]     event length (s)
%   recIdx, evtIdx [nEvents x 1]     recording / event number

switch zMode
    case 'perEvent'
        zField = 'eventZSingle';
    case 'concatenated'
        zField = 'eventZConcat';
    otherwise
        error('zMode must be ''perEvent'' or ''concatenated''.');
end
if ~isfield(S, zField) || ~isfield(S, 'eventTime')
    error(['Results are missing per-event data (%s). Re-run with ' ...
        'saveResultsStructTF = true so the processing is recomputed and saved.'], zField);
end

evt = struct('energy',[],'meanEnergy',[],'maxEnergy',[],'maxROI',[], ...
    'onsetLatency',[],'firstWeight',[],'duration',[],'recIdx',[],'evtIdx',[], ...
    'roiNames',{refNames});

nROI = numel(refNames);

for r = 1:numel(S)
    Zcells = S(r).(zField);
    Tcells = S(r).eventTime;

    % ---- make sure ROIs are identical and in the same order ----
    [found, loc] = ismember(refNames, S(r).signalNames);
    if ~all(found) || numel(S(r).signalNames) ~= nROI
        warning('Recording %d: ROI set differs from reference; skipping.', r);
        continue
    end

    for k = 1:numel(Zcells)
        Z = Zcells{k};
        t = Tcells{k};
        n = size(Z,1);
        if n < 2 || numel(t) ~= n || t(end) <= 0
            continue
        end
        Z = Z(:, loc);

        dur = t(end);
        dt  = dur/(n-1);

        % ---- energy ----
        Zint = Z;
        if posOnly
            Zint = max(Z,0);
        end
        e = trapz(t, Zint, 1);
        if perSec
            e = e / dur;
        end
        [mx, mxIdx] = max(e);

        % ---- onset of significant motion ----
        nMin = max(1, round(minDur/dt));
        sustained = movmin(double(Z > thr), [0 nMin-1], 1, 'Endpoints', 0);
        lat = nan(1,nROI);
        for j = 1:nROI
            idx = find(sustained(:,j), 1, 'first');
            if ~isempty(idx)
                lat(j) = t(idx);
            end
        end
        w = zeros(1,nROI);
        if any(~isnan(lat))
            isFirst = lat == min(lat);
            w(isFirst) = 1 / sum(isFirst);
        end

        evt.energy(end+1,:)       = e;
        evt.meanEnergy(end+1,1)   = mean(e);
        evt.maxEnergy(end+1,1)    = mx;
        evt.maxROI(end+1,1)       = mxIdx;
        evt.onsetLatency(end+1,:) = lat;
        evt.firstWeight(end+1,:)  = w;
        evt.duration(end+1,1)     = dur;
        evt.recIdx(end+1,1)       = r;
        evt.evtIdx(end+1,1)       = k;
    end
end

end

function d = cliffsDelta_SF(a, b)
% Cliff's delta = P(a > b) - P(a < b); positive when a tends to be larger.
a = a(:); b = b(:);
D = sign(a - b');
d = sum(D(:)) / (numel(a)*numel(b));
end

function q = bhFDR_SF(p)
% Benjamini-Hochberg FDR-adjusted p-values (NaNs ignored).
p = p(:);
q = nan(size(p));
ok = ~isnan(p);
pv = p(ok);
n  = numel(pv);
[ps, idx] = sort(pv);
qs = ps .* n ./ (1:n)';
qs = min(1, flipud(cummin(flipud(qs))));
qv = zeros(n,1);
qv(idx) = qs;
q(ok) = qv;
end

function p = lmeStatePValue_SF(vA, recA, vR, recR)
% p-value for Awake vs REM from a linear mixed model with a random
% intercept per recording (needs Statistics and Machine Learning Toolbox).
p = NaN;
try
    y     = [vA(:); vR(:)];
    state = categorical([zeros(numel(vA),1); ones(numel(vR),1)], [0 1], {'Awake','REM'});
    rec   = categorical([recA(:); recR(:) + 1000]);   % keep Awake / REM recording IDs distinct
    tbl   = table(y, state, rec);
    lme   = fitlme(tbl, 'y ~ state + (1|rec)');
    p     = lme.Coefficients.pValue(2);
catch
end
end

function [chi2, p] = permChi2Counts_SF(wA, wR, nPerm)
% Permutation chi-square test of whether the distribution over ROIs differs
% between two groups of events. Each row of wA / wR is one event's
% (possibly fractional) assignment to ROIs; all-zero rows are ignored.
W    = [wA; wR];
grp  = [false(size(wA,1),1); true(size(wR,1),1)];
keep = any(W > 0, 2);
W    = W(keep,:);
grp  = grp(keep);

chi2 = chi2FromWeights_SF(W, grp);
permStat = zeros(nPerm,1);
for i = 1:nPerm
    permStat(i) = chi2FromWeights_SF(W, grp(randperm(numel(grp))));
end
p = (1 + sum(permStat >= chi2)) / (1 + nPerm);
end

function c = chi2FromWeights_SF(W, grp)
O = [sum(W(~grp,:),1); sum(W(grp,:),1)];
E = sum(O,2) * sum(O,1) / sum(O(:));
m = E > 0;
c = sum((O(m) - E(m)).^2 ./ E(m));
end

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
    XRaw = X;        % un-zscored motion, kept for per-event z-scoring
    X = zscore(X);   % <<< critical for ROI comparisons
    Y = fillmissing(Y,'linear');
    YRaw = Y;        % un-zscored motion, kept for per-event z-scoring
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
    % per-event data (fresh for every recording)
    eventZConcat  = cell(1,length(segmentIdx)-1);   % z-score over concatenated events
    eventZSingle  = cell(1,length(segmentIdx)-1);   % z-score of each event on its own
    eventTime     = cell(1,length(segmentIdx)-1);   % time within event (s)
    eventZConcatB = cell(1,length(segmentIdxBaseline)-1);
    eventZSingleB = cell(1,length(segmentIdxBaseline)-1);
    eventTimeB    = cell(1,length(segmentIdxBaseline)-1);
    for n = 1:length(segmentIdx)-1
        currSegIdx = segmentIdx(n)+1:segmentIdx(n+1);
        currSegIdxBaseline = segmentIdxBaseline(n)+1:segmentIdxBaseline(n+1);
        SegX = X(currSegIdx,:);
        SegY = Y(currSegIdxBaseline,:);
        eventZConcat{n}  = SegX;
        eventZSingle{n}  = zscore(XRaw(currSegIdx,:));
        tSeg = timeVec(currSegIdx);
        eventTime{n}     = tSeg(:) - tSeg(1);
        eventZConcatB{n} = SegY;
        eventZSingleB{n} = zscore(YRaw(currSegIdxBaseline,:));
        tSegB = timeVecBaseline(currSegIdxBaseline);
        eventTimeB{n}    = tSegB(:) - tSegB(1);
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
    results(r).eventZConcat = eventZConcat;
    results(r).eventZSingle = eventZSingle;
    results(r).eventTime = eventTime;
    results(r).eventZConcatB = eventZConcatB;
    results(r).eventZSingleB = eventZSingleB;
    results(r).eventTimeB = eventTimeB;
    
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