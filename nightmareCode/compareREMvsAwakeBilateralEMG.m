function results = compareREMvsAwakeBilateralEMG(excelPath)

%% Parameters

awakeDur = 90;      % seconds after REM offset
maxLagSec = 5;      % cross-correlation lag
nDensityBins = 75;

%% Read Excel

T = readtable(excelPath,'Sheet',2,'ReadVariableNames',false);

if strcmpi(T{1,1}{1},'file name')
    T = T(2:end,:);
end

nRec = height(T);

%% Storage

REMxcorrs = {};
Awakexcorrs = {};

REMzeroLag = [];
AwakezeroLag = [];

REMcoh = {};
Awakecoh = {};

remLeftPool = [];
remRightPool = [];

awakeLeftPool = [];
awakeRightPool = [];

REM_AI = [];
Awake_AI = [];

REM_AI_traces = {};
Awake_AI_traces = {};

REM_ratio = [];
Awake_ratio = [];

risingEdgesREM = [];
risingEdgesAwake = [];

awakeColor = [0 0.4470 0.7410];
remColor   = [0.8500 0.3250 0.0980];

%% Loop recordings

for r = 1:nRec
    tdmsPath = T{r,1}{1};

    segments = table2array(T(r,2:end));
    segments = segments(~isnan(segments));

    if isempty(segments)
        continue
    end

    segments = reshape(segments,2,[])';

    procDataPath = strrep(tdmsPath,'.tdms','_ProcData.mat');

    if ~exist(procDataPath,'file')
        warning('Missing ProcData file: %s',procDataPath);
        continue
    end

    S = load(procDataPath,'ProcData');
    ProcData = S.ProcData;

    Fs = ProcData.notes.dsFs;

    leftPower  = ProcData.EMG.LeftPower(:);
    rightPower = ProcData.EMG.RightPower(:);

    leftSignal  = ProcData.EMG.LeftSignal(:);
    rightSignal = ProcData.EMG.RightSignal(:);

    nSamples = length(leftPower);

    maxLagSamples = round(maxLagSec*Fs);

    %% Loop REM bouts

    for s = 1:size(segments,1)

        remStart = segments(s,1);
        remStop  = segments(s,2);

        awakeStart = remStop;
        awakeStop  = remStop + awakeDur;

        remInd = round(remStart*Fs)+1 : min(round(remStop*Fs),nSamples);

        awakeInd = round(awakeStart*Fs)+1 : ...
            min(round(awakeStop*Fs),nSamples);

        if numel(remInd) < Fs*5
            continue
        end

        if numel(awakeInd) < Fs*5
            continue
        end

        %% =========================
        %% REM CROSS-CORRELATION
        %% =========================

        L = leftPower(remInd);
        R = rightPower(remInd);

        % Mean removal
        LSub = L - mean(L);
        RSub = R - mean(R);

        [xc,lags] = xcorr(LSub,RSub,maxLagSamples,'coeff');

        REMxcorrs{end+1} = xc(:)';

        REMzeroLag(end+1,1) = xc(lags==0);

        %% =========================
        %% REM AI
        %% =========================
        
        L_z = zscore(L-(min(L)-1));
        R_z = zscore(R-(min(R)-1));

        % AI = abs(L_z-R_z)./(L_z+R_z+eps);
        AI = abs(L_z - R_z);
        % threshold = mean(AI) + 2.5*std(AI);
        threshold = 3;
        aboveThresh = AI > threshold;
        risingEdgesREM(end+1) = sum(diff([0; aboveThresh(:)]) == 1);

        figure
        subplot(5,1,1)
        plot(L)
        title('EMG Left Power')
        subplot(5,1,2)
        plot(L_z)
        title('EMG zscore Left Power')
        subplot(5,1,3)
        plot(R)
        title('EMG right Power')
        subplot(5,1,4)
        plot(R_z)
        title('EMG zscore right Power')
        subplot(5,1,5)
        plot(AI)
        title('EMG zscore power diff')
        hold on
        plot([1 length(AI)],[threshold threshold],'w')
        ylim([0 10])
        close
        
        REM_AI_traces{end+1} = AI;
        REM_AI(end+1) = std(AI);

        ratio = log10((L+eps)./(R+eps));
        REM_ratio = [REM_ratio; ratio(:)];

        %% =========================
        %% AWAKE CROSS-CORRELATION
        %% =========================

        L = leftPower(awakeInd);
        R = rightPower(awakeInd);

        % Mean removal
        LSub = L - mean(L);
        RSub = R - mean(R);

        [xc,lags] = xcorr(LSub,RSub,maxLagSamples,'coeff');

        Awakexcorrs{end+1} = xc(:)';

        AwakezeroLag(end+1,1) = xc(lags==0);


        %% =========================
        %% AWAKE AI
        %% =========================

        L_z = zscore(L-(min(L)-1));
        R_z = zscore(R-(min(R)-1));

        % AI = abs(L_z-R_z)./(L_z+R_z+eps);
        AI = abs(L_z - R_z);
        % threshold = mean(AI) + 2.5*std(AI);
        threshold = 3;
        aboveThresh = AI > threshold;
        risingEdgesAwake(end+1) = sum(diff([0; aboveThresh(:)]) == 1);

        figure
        subplot(5,1,1)
        plot(L)
        title('EMG Left Power')
        subplot(5,1,2)
        plot(L_z)
        title('EMG zscore Left Power')
        subplot(5,1,3)
        plot(R)
        title('EMG right Power')
        subplot(5,1,4)
        plot(R_z)
        title('EMG zscore right Power')
        subplot(5,1,5)
        plot(AI)
        title('EMG zscore power diff')
        hold on
        plot([1 length(AI)],[threshold threshold],'w')
        ylim([0 10])
        close


        Awake_AI_traces{end+1} = AI;
        Awake_AI(end+1) = std(AI);

        ratio = log10((L+eps)./(R+eps));
        Awake_ratio = [Awake_ratio; ratio(:)];

        %% =========================
        %% REM COHERENCE
        %% =========================

        Lsig = leftSignal(remInd);
        Rsig = rightSignal(remInd);

        [Cxy,f] = mscohere(Lsig,...
            Rsig,...
            hamming(1024),...
            512,...
            1024,...
            Fs);

        REMcoh{end+1} = Cxy(:)';

        %% =========================
        %% AWAKE COHERENCE
        %% =========================

        Lsig = leftSignal(awakeInd);
        Rsig = rightSignal(awakeInd);

        [Cxy,f] = mscohere(Lsig,...
            Rsig,...
            hamming(1024),...
            512,...
            1024,...
            Fs);

        Awakecoh{end+1} = Cxy(:)';

        %% =========================
        %% DENSITY DATA
        %% =========================

        L = leftPower(remInd);
        R = rightPower(remInd);

        % L = L./median(L);
        % R = R./median(R);
        L = zscore(L);
        R = zscore(R);

        remLeftPool  = [remLeftPool;  L(:)];
        remRightPool = [remRightPool; R(:)];

        L = leftPower(awakeInd);
        R = rightPower(awakeInd);

        % L = L./median(L);
        % R = R./median(R);
        L = zscore(L);
        R = zscore(R);

        awakeLeftPool  = [awakeLeftPool;  L(:)];
        awakeRightPool = [awakeRightPool; R(:)];

    end

end

%% Convert to matrices

REMxcorrMat   = vertcat(REMxcorrs{:});
AwakexcorrMat = vertcat(Awakexcorrs{:});

REMcohMat   = vertcat(REMcoh{:});
AwakecohMat = vertcat(Awakecoh{:});

lagSec = lags/Fs;

%% ====================================================
%% FIGURE 1 - Cross-correlation
%% ====================================================

figure;

subplot(2,1,1)
hold on

plot(lagSec,REMxcorrMat','Color', remColor)


% 95% CI
remCI = 1.96 * (std(REMxcorrMat,[],1) / sqrt(size(REMxcorrMat,1)));
fill([lagSec fliplr(lagSec)], ...
     [mean(REMxcorrMat,1)+remCI fliplr(mean(REMxcorrMat,1)-remCI)], ...
     [0.6 0.6 0.6], ...
     'FaceAlpha', 0.45, ...
     'EdgeColor', 'none');


plot(lagSec,...
    mean(REMxcorrMat,1),...
    'w',...
    'LineWidth',2)

xlabel('Lag (s)')
ylabel('Corr')
title('REM')
ylim([-.4 1])

subplot(2,1,2)
hold on

plot(lagSec,AwakexcorrMat','Color',awakeColor)

% 95% CI
awakeCI = 1.96 * (std(AwakexcorrMat,[],1) / sqrt(size(AwakexcorrMat,1)));
fill([lagSec fliplr(lagSec)], ...
     [mean(AwakexcorrMat,1)+awakeCI fliplr(mean(AwakexcorrMat,1)-awakeCI)], ...
     [0.6 0.6 0.6], ...
     'FaceAlpha', 0.45, ...
     'EdgeColor', 'none');

plot(lagSec,...
    mean(AwakexcorrMat,1),...
    'w',...
    'LineWidth',2)

xlabel('Lag (s)')
ylabel('Corr')
title('Awake')
ylim([-.4 1])

%% ====================================================
%% FIGURE 2 - Zero lag correlation
%% ====================================================

figure
hold on

plotSpread({REMzeroLag,AwakezeroLag});

means = [mean(REMzeroLag) mean(AwakezeroLag)];
sems  = [std(REMzeroLag)/sqrt(numel(REMzeroLag)) ...
    std(AwakezeroLag)/sqrt(numel(AwakezeroLag))];

errorbar(1:2,...
    means,...
    sems,...
    'k',...
    'LineStyle','none',...
    'LineWidth',2);

set(gca,'XTick',[1 2],...
    'XTickLabel',{'REM','Awake'})

ylabel('Zero-lag correlation')

%% ====================================================
%% FIGURE 3 - Coherence
%% ====================================================

figure
hold on

plot(f,...
    mean(REMcohMat,1),...
    'LineWidth',2)

plot(f,...
    mean(AwakecohMat,1),...
    'LineWidth',2)

xlabel('Frequency (Hz)')
ylabel('Coherence')

legend({'REM','Awake'})
box off

%% ====================================================
%% FIGURE 4 - Density plots
%% ====================================================

figure

subplot(1,2,1)

histogram2(remLeftPool,...
    remRightPool,...
    nDensityBins,...
    'DisplayStyle','tile',...
    'ShowEmptyBins','off');

xlabel('Left Power (zscore)')
ylabel('Right Power (zscore)')
title('REM')

colorbar

% xlim([0 2])
% ylim([0 2])

subplot(1,2,2)

histogram2(awakeLeftPool,...
    awakeRightPool,...
    nDensityBins,...
    'DisplayStyle','tile',...
    'ShowEmptyBins','off');

xlabel('Left Power (zscore)')
ylabel('Right Power (zscore)')
title('Awake')

colorbar

% xlim([0 2])
% ylim([0 2])

%% ====================================================
%% FIGURE 5 - Asymmetry index STD
%% ====================================================

figure;

plotSpread({REM_AI(:),Awake_AI(:)},...
    'distributionColors',[remColor;awakeColor]);

hold on

errorbar(1,mean(REM_AI),...
    std(REM_AI)/sqrt(length(REM_AI)),...
    'wo','LineWidth',2,'MarkerFaceColor','w')

errorbar(2,mean(Awake_AI),...
    std(Awake_AI)/sqrt(length(Awake_AI)),...
    'wo','LineWidth',2,'MarkerFaceColor','w')

xlim([0.5 2.5])

set(gca,'XTick',[1 2],...
    'XTickLabel',{'REM','Awake'})

ylabel('Standard Deviation of EMG Power Difference')
title('Left-Right EMG Asymmetry')
box off

%% ====================================================
%% FIGURE 6 - Log ratio distribution
%% ====================================================

figure

edges = -3:0.1:3;

histogram(REM_ratio,...
    edges,...
    'Normalization','probability')

hold on

histogram(Awake_ratio,...
    edges,...
    'Normalization','probability')

xlabel('log_{10}(Left / Right EMG Power)')
ylabel('Probability')

legend({'REM','Awake'})

title('Bilateral Power Ratio Distribution')

box off

%% ====================================================
%% FIGURE 7 - AI traces
%% ====================================================

figure;

%% REM
subplot(1,2,1)
hold on

for i = 1:length(REM_AI_traces)

    AI = REM_AI_traces{i};

    x = linspace(0,100,length(AI));

    plot(x,AI,'LineWidth',1);

end

xlabel('Percent Through Event')
ylabel('Asymmetry Index')
title('REM')
xlim([0 100])
% ylim([-1000 1000])

%% Awake
subplot(1,2,2)
hold on

for i = 1:length(Awake_AI_traces)

    AI = Awake_AI_traces{i};

    x = linspace(0,100,length(AI));

    plot(x,AI,'LineWidth',1);

end

xlabel('Percent Through Event')
ylabel('Asymmetry Index')
title('Awake')
xlim([0 100])

% ylim([-1000 1000])

%% ====================================================
%% FIGURE 8 - Threshold crossing histogram
%% ====================================================

figure

%% Histogram + Poisson fits
subplot(1,2,1)
hold on

allCounts = [risingEdgesAwake(:); risingEdgesREM(:)];

edges = (-0.5):(max(allCounts)+0.5);

% Histograms
histogram(risingEdgesAwake,...
    edges,...
    'Normalization','pdf',...
    'FaceAlpha',0.4, 'FaceColor',awakeColor);

histogram(risingEdgesREM,...
    edges,...
    'Normalization','pdf',...
    'FaceAlpha',0.4,'FaceColor',remColor);

% Poisson fits
lambdaAwake = mean(risingEdgesAwake);
lambdaREM   = mean(risingEdgesREM);

xFit = 0:max(allCounts);

plot(xFit,...
    poisspdf(xFit,lambdaAwake),...
    'LineWidth',3,'Color',awakeColor)

plot(xFit,...
    poisspdf(xFit,lambdaREM),...
    'LineWidth',3,'Color',remColor)

xlabel('Threshold Crossing Count')
ylabel('Probability')
title('Crossing Count Distribution')
legend('Awake','REM','Awake Poisson','REM Poisson')
box off

%% plotSpread
subplot(1,2,2)
hold on

plotSpread({risingEdgesAwake,risingEdgesREM},...
    'distributionColors',[awakeColor; remColor]);

s = findobj(gca,'Type','Scatter');

for k = 1:length(s)
    s(k).SizeData = 100;   % increase marker size
end

means = [mean(risingEdgesAwake) mean(risingEdgesREM)];
sems = [ ...
    std(risingEdgesAwake)/sqrt(numel(risingEdgesAwake)), ...
    std(risingEdgesREM)/sqrt(numel(risingEdgesREM))];

errorbar(1:2,...
    means,...
    sems,...
    'w.',...
    'LineWidth',2,...
    'CapSize',5)

plot(1:2,...
    means,...
    'wo',...
    'MarkerFaceColor','w',...
    'MarkerSize',5)

xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Awake','REM'})

ylabel('Threshold Crossing Count')
title('Per Event Counts')

box off

sgtitle('Large Bilateral Difference Events (ZScore EMG Power Diff > 3)')

%% Output structure

results.lagSec = lagSec;
results.xcorrREM = REMxcorrMat;
results.xcorrAwake = AwakexcorrMat;

results.zeroLagREM = REMzeroLag;
results.zeroLagAwake = AwakezeroLag;

results.f = f;
results.cohREM = REMcohMat;
results.cohAwake = AwakecohMat;

results.remLeftPool = remLeftPool;
results.remRightPool = remRightPool;

results.awakeLeftPool = awakeLeftPool;
results.awakeRightPool = awakeRightPool;
end
