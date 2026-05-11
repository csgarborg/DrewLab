function compareREMvsAwake(awakeExcelPath, remExcelPath, mefRemExcelPath, tdmsTF, roiPath)

%% Load data

if ~exist("tdmsTF","var") || tdmsTF
    awake  = processExcelSegmentsTDMS_SF(awakeExcelPath);
    rem    = processExcelSegmentsTDMS_SF(remExcelPath);
    mefRem = processExcelSegmentsTDMS_SF(mefRemExcelPath);
else
    awake  = processExcelSegmentsROISelect_SF(awakeExcelPath,roiPath);
    rem    = processExcelSegmentsROISelect_SF(remExcelPath,roiPath);
    mefRem = processExcelSegmentsROISelect_SF(mefRemExcelPath,roiPath);

    % show ROI locations once
    plotFirstFrameWithROIs(remExcelPath, roiPath);
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

function results = processExcelSegmentsROISelect_SF(excelPath, roiPath)

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
    
    % Remove time field if present
    if isfield(roiMotion,'time')
        roiMotion = rmfield(roiMotion,'time');
    end
    
    % Enforce consistent ordering across recordings
    signalNames = sort(fieldnames(roiMotion));

    % % Remove whisker fields
    % idx = cellfun(@(s) contains(s,'whisker','IgnoreCase',true), signalNames);
    % signalNames(idx) = [];

    nSignals = numel(signalNames);
    
    % ==============================
    % 3) Build matrix (UNCHANGED LOGIC)
    % ==============================
    L = min(structfun(@length, roiMotion));
    X = zeros(L,nSignals);
    
    for i = 1:nSignals
        X(:,i) = roiMotion.(signalNames{i})(1:L);
    end

    % ====== Preprocessing (IMPORTANT) ======
    X = fillmissing(X,'linear');
    X = zscore(X);   % <<< critical for ROI comparisons
    
    % ==============================
    % 4) PCA (UNCHANGED)
    % ==============================
    [coeff,score,latent,~,explained] = pca(X);
    
    % ==============================
    % 5) Correlation (UNCHANGED)
    % ==============================
    corrMatrix = corrcoef(X);
    
    % ==============================
    % 6) Store (UNCHANGED)
    % ==============================
    results(r).coeff = coeff;
    results(r).explained = explained;
    results(r).corrMatrix = corrMatrix;
    results(r).score = score;
    results(r).signalNames = signalNames;
    results(r).nSignals = nSignals;
    
end

end

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