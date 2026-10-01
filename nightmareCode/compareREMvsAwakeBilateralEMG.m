function results = compareREMvsAwakeBilateralEMG(excelPath, roiPath)

%% Parameters

awakeDur = 90;      % seconds after REM offset
maxLagSec = 5;      % cross-correlation lag
nDensityBins = 75;

% =========================
% ROI trigger parameters
% =========================
roiWin = [-5 5];              % sec around REM async rising edge
roiPostWin = [0 1];           % summary window for ROI response
roiPreWin  = [-5 -1];         % baseline window
minEventSeparationSec = 1;    % merge nearby rising edges if desired

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

% Direction of REM threshold crossings
REM_leftCrossings  = [];
REM_rightCrossings = [];

Awake_leftCrossings  = [];
Awake_rightCrossings = [];

REM_crossingsPerMin = [];
Awake_crossingsPerMin = [];

awakeColor = [0 0.4470 0.7410];
remColor   = [0.8500 0.3250 0.0980];

% =========================
% ROI storage
% =========================
allROINames = {};
roiEventData = struct();        % roiEventData.(roiName){k} = peri-event trace
roiPeakVals = struct();         % roiPeakVals.(roiName) = vector of event peak responses
roiMeanVals = struct();         % roiMeanVals.(roiName) = vector of event mean responses
roiEventCount = 0;

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

    % ============================================================
    % ROI LOAD (same style as compareREMvsAwake pipeline)
    % ============================================================
    roiMotion = [];
    roiTime = [];
    roiNames = {};

    if nargin >= 2 && ~isempty(roiPath)
        try
            [folder,name,~] = fileparts(tdmsPath);
            videoPath = fullfile(folder,[name '.mp4']);

            roiMotion = extractROIMotionFromVideo_fast(videoPath, excelPath, roiPath, segments);

            if isstruct(roiMotion) && isfield(roiMotion,'time')
                roiTime = roiMotion.time(:);
                roiMotion = rmfield(roiMotion,'time');
                roiNames = sort(fieldnames(roiMotion));

                % --------------------------------------------------------
                % Z-SCORE ROI DATA AFTER LOADING
                % --------------------------------------------------------
                for rr = 1:numel(roiNames)
                    x = roiMotion.(roiNames{rr})(:);
                    x = fillmissing(x,'linear');
                    if std(x,'omitnan') > 0
                        x = zscore(x);
                    else
                        x = x - mean(x,'omitnan');
                    end
                    roiMotion.(roiNames{rr}) = x;
                end

                % initialize global ROI fields if first time seen
                for rr = 1:numel(roiNames)
                    rn = roiNames{rr};
                    if ~ismember(rn,allROINames)
                        allROINames{end+1} = rn; %#ok<AGROW>
                    end
                    if ~isfield(roiEventData,rn)
                        roiEventData.(rn) = {};
                        roiPeakVals.(rn) = [];
                        roiMeanVals.(rn) = [];
                    end
                end
            else
                warning('ROI motion returned empty or without time for %s',tdmsPath);
            end
        catch ME
            warning('ROI extraction/loading failed for %s\n%s',tdmsPath,ME.message);
            roiMotion = [];
            roiTime = [];
            roiNames = {};
        end
    end

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

        % %% =========================
        % %% REM AI
        % %% =========================
        % 
        % L_z = zscore(L-(min(L)-1));
        % R_z = zscore(R-(min(R)-1));
        % 
        % AI = abs(L_z - R_z);
        % threshold = 3;
        % aboveThresh = AI > threshold;
        % risingEdgesREM(end+1) = sum(diff([0; aboveThresh(:)]) == 1);
        % 
        % REM_AI_traces{end+1} = AI;
        % REM_AI(end+1) = std(AI);
        % 
        % ratio = log10((L+eps)./(R+eps));
        % REM_ratio = [REM_ratio; ratio(:)];

        %% =========================
        %% REM AI
        %% =========================

        L_z = zscore(L-(min(L)-1));
        R_z = zscore(R-(min(R)-1));

        % Signed difference:
        %   positive = Left > Right
        %   negative = Right > Left
        signedDiff = L_z - R_z;

        % Absolute asymmetry index
        AI = abs(signedDiff);

        threshold = 3;
        aboveThresh = AI > threshold;

        % Find threshold crossings
        riseIdx = find(diff([0; aboveThresh(:)]) == 1);

        % Total number of large bilateral difference events
        risingEdgesREM(end+1) = numel(riseIdx);

        % Determine which side was responsible at each crossing
        leftCrossings  = sum(signedDiff(riseIdx) > threshold);
        rightCrossings = sum(signedDiff(riseIdx) < -threshold);

        REM_leftCrossings(end+1)  = leftCrossings;
        REM_rightCrossings(end+1) = rightCrossings;

        REM_AI_traces{end+1} = AI;
        REM_AI(end+1) = std(AI);

        ratio = log10((L+eps)./(R+eps));
        REM_ratio = [REM_ratio; ratio(:)];

        % Crossings normalized to REM duration
        remDurationMin = numel(remInd) / Fs / 60;
        REM_crossingsPerMin(end+1) = numel(riseIdx) / remDurationMin;

        % ============================================================
        % ROI EVENT-TRIGGERED ANALYSIS FOR REM ASYNC EVENTS
        % ============================================================
        if ~isempty(roiTime) && ~isempty(roiNames)

            riseIdx = find(diff([0; aboveThresh(:)]) == 1);

            if ~isempty(riseIdx)
                riseTimes = remStart + (riseIdx(:)-1)./Fs;

                % optional minimum event separation
                if numel(riseTimes) > 1
                    keep = [true; diff(riseTimes) >= minEventSeparationSec];
                    riseTimes = riseTimes(keep);
                end

                % keep only events whose full ROI window lies inside this REM bout
                keep = (riseTimes + roiWin(1) >= remStart) & ...
                       (riseTimes + roiWin(2) <= remStop);
                riseTimes = riseTimes(keep);

                % restrict ROI time to current REM segment
                roiSegMask = roiTime >= remStart & roiTime <= remStop;
                roiTimeSeg = roiTime(roiSegMask);

                if ~isempty(roiTimeSeg)
                    dtROI = median(diff(roiTimeSeg));
                    if ~isfinite(dtROI) || dtROI <= 0
                        dtROI = 1/60;
                    end

                    tRel = roiWin(1):dtROI:roiWin(2);

                    for ev = 1:numel(riseTimes)
                        t0 = riseTimes(ev);

                        roiEventCount = roiEventCount + 1;

                        for rr = 1:numel(roiNames)
                            rn = roiNames{rr};

                            fullTrace = roiMotion.(rn)(:);
                            roiTraceSeg = fullTrace(roiSegMask);

                            if numel(roiTraceSeg) ~= numel(roiTimeSeg)
                                m = min(numel(roiTraceSeg),numel(roiTimeSeg));
                                roiTraceSeg_use = roiTraceSeg(1:m);
                                roiTimeSeg_use = roiTimeSeg(1:m);
                            else
                                roiTraceSeg_use = roiTraceSeg;
                                roiTimeSeg_use = roiTimeSeg;
                            end

                            % interpolate onto common peri-event axis
                            x = interp1(roiTimeSeg_use - t0,...
                                        roiTraceSeg_use,...
                                        tRel,...
                                        'linear',...
                                        NaN);

                            % % baseline subtract using pre-event window
                            % preMask = tRel >= roiPreWin(1) & tRel <= roiPreWin(2);
                            % if any(preMask) && any(isfinite(x(preMask)))
                            %     x = x - mean(x(preMask),'omitnan');
                            % end

                            roiEventData.(rn){end+1,1} = x;

                            postMask = tRel >= roiPostWin(1) & tRel <= roiPostWin(2);
                            roiPeakVals.(rn)(end+1,1) = max(x(postMask),[],'omitnan');
                            roiMeanVals.(rn)(end+1,1) = mean(x(postMask),'omitnan');
                        end
                    end
                end
            end
        end

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

        % %% =========================
        % %% AWAKE AI
        % %% =========================
        % 
        % L_z = zscore(L-(min(L)-1));
        % R_z = zscore(R-(min(R)-1));
        % 
        % AI = abs(L_z - R_z);
        % threshold = 3;
        % aboveThresh = AI > threshold;
        % risingEdgesAwake(end+1) = sum(diff([0; aboveThresh(:)]) == 1);
        % 
        % Awake_AI_traces{end+1} = AI;
        % Awake_AI(end+1) = std(AI);
        %
        % ratio = log10((L+eps)./(R+eps));
        % Awake_ratio = [Awake_ratio; ratio(:)];

        %% =========================
        %% AWAKE AI
        %% =========================

        L_z = zscore(L-(min(L)-1));
        R_z = zscore(R-(min(R)-1));

        signedDiff = L_z - R_z;
        AI = abs(signedDiff);

        threshold = 3;
        aboveThresh = AI > threshold;

        riseIdx = find(diff([0; aboveThresh(:)]) == 1);

        risingEdgesAwake(end+1) = numel(riseIdx);

        leftCrossings  = sum(signedDiff(riseIdx) > threshold);
        rightCrossings = sum(signedDiff(riseIdx) < -threshold);

        Awake_leftCrossings(end+1)  = leftCrossings;
        Awake_rightCrossings(end+1) = rightCrossings;

        Awake_AI_traces{end+1} = AI;
        Awake_AI(end+1) = std(AI);

        ratio = log10((L+eps)./(R+eps));
        Awake_ratio = [Awake_ratio; ratio(:)];

        % Crossings normalized to REM duration
        awakeDurationMin = numel(awakeInd) / Fs / 60;
        Awake_crossingsPerMin(end+1) = numel(riseIdx) / awakeDurationMin;

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

        L = zscore(L);
        R = zscore(R);

        remLeftPool  = [remLeftPool;  L(:)];
        remRightPool = [remRightPool; R(:)];

        L = leftPower(awakeInd);
        R = rightPower(awakeInd);

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

remCI = 1.96 * (std(REMxcorrMat,[],1) / sqrt(size(REMxcorrMat,1)));
fill([lagSec fliplr(lagSec)], ...
     [mean(REMxcorrMat,1)+remCI fliplr(mean(REMxcorrMat,1)-remCI)], ...
     [0.6 0.6 0.6], ...
     'FaceAlpha', 0.45, ...
     'EdgeColor', 'none');

plot(lagSec, mean(REMxcorrMat,1), 'w', 'LineWidth',2)
xlabel('Lag (s)')
ylabel('Corr')
title('REM')
ylim([-.4 1])

subplot(2,1,2)
hold on
plot(lagSec,AwakexcorrMat','Color',awakeColor)

awakeCI = 1.96 * (std(AwakexcorrMat,[],1) / sqrt(size(AwakexcorrMat,1)));
fill([lagSec fliplr(lagSec)], ...
     [mean(AwakexcorrMat,1)+awakeCI fliplr(mean(AwakexcorrMat,1)-awakeCI)], ...
     [0.6 0.6 0.6], ...
     'FaceAlpha', 0.45, ...
     'EdgeColor', 'none');

plot(lagSec, mean(AwakexcorrMat,1), 'w', 'LineWidth',2)
xlabel('Lag (s)')
ylabel('Corr')
title('Awake')
ylim([-.4 1])

%% ====================================================
%% FIGURE 2 - Zero lag correlation
%% ====================================================

figure
hold on

plotSpread({REMzeroLag,AwakezeroLag},'distributionColors',[remColor;awakeColor]);

means = [mean(REMzeroLag) mean(AwakezeroLag)];
sems  = [std(REMzeroLag)/sqrt(numel(REMzeroLag)) ...
    std(AwakezeroLag)/sqrt(numel(AwakezeroLag))];

errorbar(1:2, means, sems, 'w', 'LineStyle','none', 'LineWidth',2);

set(gca,'XTick',[1 2], 'XTickLabel',{'REM','Awake'})
ylabel('Zero-lag correlation')

%% ====================================================
%% FIGURE 3 - Coherence
%% ====================================================

figure
hold on

plot(f, mean(REMcohMat,1), 'LineWidth',2)
plot(f, mean(AwakecohMat,1), 'LineWidth',2)

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
set(gca,'XTick',[1 2], 'XTickLabel',{'REM','Awake'})
ylabel('Standard Deviation of EMG Power Difference')
title('Left-Right EMG Asymmetry')
box off

%% ====================================================
%% FIGURE 6 - Log ratio distribution
%% ====================================================

figure

edges = -3:0.1:3;

histogram(REM_ratio, edges, 'Normalization','probability')
hold on
histogram(Awake_ratio, edges, 'Normalization','probability')

xlabel('log_{10}(Left / Right EMG Power)')
ylabel('Probability')
legend({'REM','Awake'})
title('Bilateral Power Ratio Distribution')
box off

%% ====================================================
%% FIGURE 7 - AI traces
%% ====================================================

figure;

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

%% ====================================================
%% FIGURE 8 - Threshold crossing histogram
%% ====================================================

figure

subplot(2,2,1)
hold on

allCounts = [risingEdgesAwake(:); risingEdgesREM(:)];
edges = (-0.5):(max(allCounts)+0.5);

histogram(risingEdgesAwake,...
    edges,...
    'Normalization','pdf',...
    'FaceAlpha',0.4, 'FaceColor',awakeColor);

histogram(risingEdgesREM,...
    edges,...
    'Normalization','pdf',...
    'FaceAlpha',0.4,'FaceColor',remColor);

lambdaAwake = mean(risingEdgesAwake);
lambdaREM   = mean(risingEdgesREM);

xFit = 0:max(allCounts);

plot(xFit, poisspdf(xFit,lambdaAwake), 'LineWidth',3,'Color',awakeColor)
plot(xFit, poisspdf(xFit,lambdaREM),   'LineWidth',3,'Color',remColor)

xlabel('Threshold Crossing Count')
ylabel('Probability')
title('Crossing Count Distribution')
legend('Awake','REM','Awake Poisson','REM Poisson')
box off

subplot(2,2,2)
hold on

plotSpread({risingEdgesAwake,risingEdgesREM},...
    'distributionColors',[awakeColor; remColor]);

s = findobj(gca,'Type','Scatter');
for k = 1:length(s)
    s(k).SizeData = 100;
end

means = [mean(risingEdgesAwake) mean(risingEdgesREM)];
sems = [ ...
    std(risingEdgesAwake)/sqrt(numel(risingEdgesAwake)), ...
    std(risingEdgesREM)/sqrt(numel(risingEdgesREM))];

errorbar(1:2, means, sems, 'w.', 'LineWidth',2, 'CapSize',5)
plot(1:2, means, 'wo', 'MarkerFaceColor','w', 'MarkerSize',5)

xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Awake','REM'})

ylabel('Threshold Crossing Count')
title('Per Event Counts')
box off

sgtitle('Large Bilateral Difference Events (ZScore EMG Power Diff > 3)')

subplot(2,2,3)
hold on

allCounts = [Awake_crossingsPerMin(:); REM_crossingsPerMin(:)];
edges = (-0.5):(max(allCounts)+0.5);

histogram(Awake_crossingsPerMin,...
    edges,...
    'Normalization','pdf',...
    'FaceAlpha',0.4, 'FaceColor',awakeColor);

histogram(REM_crossingsPerMin,...
    edges,...
    'Normalization','pdf',...
    'FaceAlpha',0.4,'FaceColor',remColor);

lambdaAwake = mean(Awake_crossingsPerMin);
lambdaREM   = mean(REM_crossingsPerMin);

xFit = 0:max(allCounts);

plot(xFit, poisspdf(xFit,lambdaAwake), 'LineWidth',3,'Color',awakeColor)
plot(xFit, poisspdf(xFit,lambdaREM),   'LineWidth',3,'Color',remColor)

xlabel('Threshold Crossing Count/min')
ylabel('Probability')
title('Crossing Count Distribution')
legend('Awake','REM','Awake Poisson','REM Poisson')
box off

subplot(2,2,4)
hold on

plotSpread({Awake_crossingsPerMin,REM_crossingsPerMin},...
    'distributionColors',[awakeColor; remColor]);

s = findobj(gca,'Type','Scatter');
for k = 1:length(s)
    s(k).SizeData = 100;
end

means = [mean(Awake_crossingsPerMin,'omitnan'),...
         mean(REM_crossingsPerMin,'omitnan')];

sems = [...
    std(Awake_crossingsPerMin,'omitnan') / ...
        sqrt(sum(isfinite(Awake_crossingsPerMin))),...
    std(REM_crossingsPerMin,'omitnan') / ...
        sqrt(sum(isfinite(REM_crossingsPerMin)))];

errorbar(1:2,...
    means,...
    sems,...
    'w.',...
    'LineWidth',2,...
    'CapSize',5);

plot(1:2,...
    means,...
    'wo',...
    'MarkerFaceColor','w',...
    'MarkerSize',5);

xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Awake','REM'})

ylabel('Threshold Crossings / min')
title('Crossing Rate')
box off

sgtitle('Large Bilateral Difference Events (ZScore EMG Power Diff > 3)')

%% ====================================================
%% FIGURE 9 - Which side causes the asymmetry?
%% ====================================================

figure
subplot(1,2,1)
hold on

plotSpread({REM_leftCrossings(:), REM_rightCrossings(:)}, ...
    'distributionColors',[remColor; remColor]);

% Mean ± SEM
means = [mean(REM_leftCrossings,'omitnan'), ...
         mean(REM_rightCrossings,'omitnan')];

sems = [ ...
    std(REM_leftCrossings,'omitnan') / sqrt(sum(isfinite(REM_leftCrossings))), ...
    std(REM_rightCrossings,'omitnan') / sqrt(sum(isfinite(REM_rightCrossings)))];

errorbar(1:2,means,sems,...
    'wo',...
    'LineStyle','none',...
    'LineWidth',2,...
    'CapSize',5);

plot(1:2,means,'wo',...
    'MarkerFaceColor','w',...
    'MarkerSize',5);

xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Left','Right'})

ylabel('Number of threshold crossings')
title('Side Responsible for REM EMG Asymmetry')
box off

subplot(1,2,2)
hold on

nTrials = min(numel(REM_leftCrossings),...
              numel(REM_rightCrossings));

for i = 1:nTrials
    plot([1 2],...
         [REM_leftCrossings(i) REM_rightCrossings(i)],...
         '-',...
         'Color',[0.6 0.6 0.6],...
         'LineWidth',1);
    % plot([1],REM_leftCrossings(i),'Color',remColor,'');
    % plot([2],REM_rightCrossings(i),'Color',remColor);
end

plotSpread({REM_leftCrossings(:),REM_rightCrossings(:)},...
    'distributionColors',[remColor; remColor]);

s = findobj(gca,'Type','Scatter');
for k = 1:length(s)
    s(k).SizeData = 90;
end

means = [mean(REM_leftCrossings,'omitnan'),...
         mean(REM_rightCrossings,'omitnan')];

sems = [...
    std(REM_leftCrossings,'omitnan') / ...
        sqrt(sum(isfinite(REM_leftCrossings))),...
    std(REM_rightCrossings,'omitnan') / ...
        sqrt(sum(isfinite(REM_rightCrossings)))];

errorbar(1:2,...
    means,...
    sems,...
    'wo',...
    'LineStyle','none',...
    'LineWidth',2,...
    'CapSize',5);

plot(1:2,...
    means,...
    'wo',...
    'MarkerFaceColor','w',...
    'MarkerSize',5);

% bar([1 2],means)

xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Left','Right'})

ylabel('Threshold Crossings')
title('REM: Side Responsible')
box off

%% ====================================================
%% ROI FIGURE 9 - Triggered averages (mean ± SEM)
%% ====================================================
if ~isempty(allROINames)

    nROI = numel(allROINames);
    nCols = ceil(sqrt(nROI));
    nRows = ceil(nROI/nCols);

    figure('Name','ROI Triggered Averages Around REM Async EMG Events');

    for i = 1:nROI
        rn = allROINames{i};
        subplot(nRows,nCols,i); hold on

        if isfield(roiEventData,rn) && ~isempty(roiEventData.(rn))
            M = vertcat(roiEventData.(rn){:});
            M = M(~all(isnan(M),2),:);

            if ~isempty(M)
                tRel = linspace(roiWin(1),roiWin(2),size(M,2));

                mu  = mean(M,1,'omitnan');
                sem = std(M,[],1,'omitnan') ./ sqrt(sum(isfinite(M),1));

                fill([tRel fliplr(tRel)],...
                     [mu+sem fliplr(mu-sem)],...
                     remColor,...
                     'FaceAlpha',0.25,...
                     'EdgeColor','none');

                plot(tRel,mu,'Color',remColor,'LineWidth',2)
                xline(0,'w-')
            end
        end

        title(strrep(rn,'_',' '))
        xlabel('Time from REM async event (s)')
        ylabel('ROI z-score')
        box off
        ylim([-.25 .75])
    end

    sgtitle('ROI Triggered Averages Around REM Async EMG Events')
end

%% ====================================================
%% ROI FIGURE 10 - All z-scored triggered traces (semi-transparent)
%% ====================================================
if ~isempty(allROINames)

    nROI = numel(allROINames);
    nCols = ceil(sqrt(nROI));
    nRows = ceil(nROI/nCols);

    figure('Name','All ROI Triggered Traces Around REM Async EMG Events');

    for i = 1:nROI
        rn = allROINames{i};
        subplot(nRows,nCols,i); hold on

        if isfield(roiEventData,rn) && ~isempty(roiEventData.(rn))
            M = vertcat(roiEventData.(rn){:});
            M = M(~all(isnan(M),2),:);

            if ~isempty(M)
                tRel = linspace(roiWin(1),roiWin(2),size(M,2));

                for k = 1:size(M,1)
                    plot(tRel, M(k,:), ...
                        'Color',[remColor 0.15], ...
                        'LineWidth',0.8);
                end

                mu = mean(M,1,'omitnan');
                plot(tRel,mu,'Color',remColor,'LineWidth',2.5)
                xline(0,'w-')
            end
        end

        title(strrep(rn,'_',' '))
        xlabel('Time from REM async event (s)')
        ylabel('ROI z-score')
        box off
        ylim([-5 40])
    end

    sgtitle('All ROI Triggered Traces Around REM Async EMG Events')
end

%% ====================================================
%% FIGURE 11 - ROI Peak Response (plotSpread)
%% ====================================================
% Mean zScore
if ~isempty(allROINames)

    figure('Name','ROI Mean Responses After REM Async Events');
    hold on

    nROI = numel(allROINames);

    data = cell(1,nROI);
    means = nan(1,nROI);
    sems  = nan(1,nROI);

    for i = 1:nROI

        rn = allROINames{i};

        if isfield(roiPeakVals,rn)

            v = roiMeanVals.(rn);
            v = v(isfinite(v));

            data{i} = v;

            if ~isempty(v)
                means(i) = mean(v);
                sems(i)  = std(v)/sqrt(numel(v));
            end

        else

            data{i} = [];

        end

    end

    plotSpread(data,...
        'distributionColors',repmat(remColor,nROI,1));

    % Increase marker size
    s = findobj(gca,'Type','Scatter');
    for k = 1:numel(s)
        s(k).SizeData = 70;
    end

    % Mean ± SEM
    errorbar(1:nROI,...
        means,...
        sems,...
        'wo',...
        'LineStyle','none',...
        'MarkerFaceColor','w',...
        'LineWidth',2,...
        'CapSize',5);

    plot(1:nROI,...
        means,...
        'w-',...
        'LineWidth',2);

    xticks(1:nROI)
    xticklabels(strrep(allROINames,'_',' '))
    xtickangle(45)

    ylabel('Mean ROI response (z-score)')
    title('Mean Facial Motion Following REM Asynchronous EMG Events')

    box off

end

% Peak zScore
if ~isempty(allROINames)

    figure('Name','ROI Max Responses After REM Async Events');
    hold on

    nROI = numel(allROINames);

    data = cell(1,nROI);
    means = nan(1,nROI);
    sems  = nan(1,nROI);

    for i = 1:nROI

        rn = allROINames{i};

        if isfield(roiPeakVals,rn)

            v = roiPeakVals.(rn);
            v = v(isfinite(v));

            data{i} = v;

            if ~isempty(v)
                means(i) = mean(v);
                sems(i)  = std(v)/sqrt(numel(v));
            end

        else

            data{i} = [];

        end

    end

    plotSpread(data,...
        'distributionColors',repmat(remColor,nROI,1));

    % Increase marker size
    s = findobj(gca,'Type','Scatter');
    for k = 1:numel(s)
        s(k).SizeData = 70;
    end

    % Mean ± SEM
    errorbar(1:nROI,...
        means,...
        sems,...
        'wo',...
        'LineStyle','none',...
        'MarkerFaceColor','w',...
        'LineWidth',2,...
        'CapSize',5);

    plot(1:nROI,...
        means,...
        'w-',...
        'LineWidth',2);

    xticks(1:nROI)
    xticklabels(strrep(allROINames,'_',' '))
    xtickangle(45)

    ylabel('Max ROI response (z-score)')
    title('Max Facial Motion Following REM Asynchronous EMG Events')

    box off

end

%% ====================================================
%% ROI FIGURE 12 - ROI summary response mean (0-1 s post event)
%% ====================================================
if ~isempty(allROINames)

    roiMeanResp = nan(numel(allROINames),1);
    roiSemResp  = nan(numel(allROINames),1);

    for i = 1:numel(allROINames)
        rn = allROINames{i};
        if isfield(roiMeanVals,rn) && ~isempty(roiMeanVals.(rn))
            v = roiMeanVals.(rn);
            roiMeanResp(i) = mean(v,'omitnan');
            roiSemResp(i)  = std(v,'omitnan')/sqrt(sum(~isnan(v)));
        end
    end

    [~,ord] = sort(roiMeanResp,'descend','MissingPlacement','last');

    figure('Name','ROI Summary Response After REM Async Events');
    hold on

    bar(roiMeanResp(ord),'FaceColor',remColor,'FaceAlpha',0.8)
    errorbar(1:numel(ord),...
        roiMeanResp(ord),...
        roiSemResp(ord),...
        'w.',...
        'LineWidth',1.5)

    xticks(1:numel(ord))
    xticklabels(strrep(allROINames(ord),'_',' '))
    xtickangle(45)

    ylabel(sprintf('Mean ROI z-score (%.1f to %.1f s)',roiPostWin(1),roiPostWin(2)))
    title('ROI response following REM async EMG events')
    box off
end

%% ====================================================
%% ROI FIGURE 12 - ROI summary response max (0-1 s post event)
%% ====================================================
if ~isempty(allROINames)

    roiMeanResp = nan(numel(allROINames),1);
    roiSemResp  = nan(numel(allROINames),1);

    for i = 1:numel(allROINames)
        rn = allROINames{i};
        if isfield(roiPeakVals,rn) && ~isempty(roiPeakVals.(rn))
            v = roiPeakVals.(rn);
            roiMeanResp(i) = mean(v,'omitnan');
            roiSemResp(i)  = std(v,'omitnan')/sqrt(sum(~isnan(v)));
        end
    end

    [~,ord] = sort(roiMeanResp,'descend','MissingPlacement','last');

    figure('Name','ROI Summary Response After REM Async Events');
    hold on

    bar(roiMeanResp(ord),'FaceColor',remColor,'FaceAlpha',0.8)
    errorbar(1:numel(ord),...
        roiMeanResp(ord),...
        roiSemResp(ord),...
        'w.',...
        'LineWidth',1.5)

    xticks(1:numel(ord))
    xticklabels(strrep(allROINames(ord),'_',' '))
    xtickangle(45)

    ylabel(sprintf('Max ROI z-score (%.1f to %.1f s)',roiPostWin(1),roiPostWin(2)))
    title('ROI response following REM async EMG events')
    box off
end

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

results.roiEventData = roiEventData;
results.roiPeakVals = roiPeakVals;
results.roiMeanVals = roiMeanVals;
results.allROINames = allROINames;
results.roiEventCount = roiEventCount;

end