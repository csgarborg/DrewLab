function results = analyzeREMKeypointDisplacement(excelPath, csvFolder, fps, options)
% ANALYZEREMKEYPOINTDISPLACEMENT
% For each REM event listed in an excel file, finds the matching DLC
% keypoint CSV, computes a 10-second-prior baseline (x,y) position for
% each keypoint, then computes the frame-by-frame displacement magnitude
% from that baseline -- across BOTH the REM event itself AND a matching
% "Awake" window (AwakeDuration seconds immediately following REM
% offset), using the SAME pre-REM baseline for both, so displacement
% magnitude is directly comparable between the two states.
%
% For each condition (REM, Awake) this computes:
%   - per-keypoint displacement traces (full resolution, no downsampling)
%   - left-vs-right symmetry/dissimilarity metrics for any "left_X"/
%     "right_X" keypoint pairs: zero-lag correlation, cross-correlation
%     lag profile, coherence spectrum, DTW distance, an asymmetry-index
%     trace, and paired peak-displacement values
%   - a cross-keypoint correlation matrix + dendrogram across ALL
%     keypoints
% ...and then produces REM-vs-Awake comparison plots for all of the above
% so you can see whether L/R symmetry (or lack of it) is a REM-specific
% phenomenon.
%
% USAGE:
%   results = analyzeREMKeypointDisplacement(excelPath, csvFolder, fps)
%   results = analyzeREMKeypointDisplacement(excelPath, csvFolder, fps, ...
%                 'Sheet', 2, 'BaselineDuration', 10, 'AwakeDuration', 30, ...
%                 'MaxLagSec', 2, 'OutputPath', 'results.mat', 'MakePlots', true)
%
% INPUTS:
%   excelPath  - path to REM excel file (Sheet 2, col 1 = tdms/recording
%                path, cols 2:end = pairs of [start end] REM segment
%                times in seconds, NaN-padded)
%   csvFolder  - folder containing DLC keypoint CSVs. CSV filenames are
%                expected to START WITH the same base filename as column 1
%                of the excel file (extension/path stripped).
%   fps        - frame rate (frames/sec) of the video/keypoint tracking.
%                REQUIRED because it is not encoded in either file.
%
% NAME-VALUE OPTIONS:
%   'Sheet'             - excel sheet to read (default 2)
%   'BaselineDuration'  - seconds prior to REM onset used as baseline (default 10)
%   'AwakeDuration'     - seconds after REM offset defining the "Awake" comparison
%                         window (default 30). Uses the SAME pre-REM baseline as
%                         the REM event it follows.
%   'MaxLagSec'         - max lag (seconds) for the L/R cross-correlation profile (default 2)
%   'CoherenceNFFT'     - FFT length for L/R coherence (default 256)
%   'OutputPath'        - if provided, saves `results` struct to this .mat path
%   'MakePlots'         - true/false, whether to generate plots (default true)
%   'LineAlpha'         - transparency of each individual event trace (default 0.15)
%   'PairwiseAnalysis'  - true/false, compute L/R pairwise metrics (default true)
%   'GlobalCorrAnalysis'- true/false, compute cross-keypoint corr matrix + dendrogram (default true)
%
% OUTPUT: results.REM and results.Awake, each with the same sub-structure:
%   .keypoints  - struct array, one per keypoint:
%       .name, .rawTraces (cell, full-res displacement magnitude per event),
%       .xNorm (cell, that event's own 0-100% x-axis), .eventInfo
%   .pairwise   - struct array, one per detected L/R keypoint pair:
%       .pairName, .leftName, .rightName
%       .zeroLagCorr, .lagsSec, .xcorrProfiles, .coherenceFreq,
%       .coherenceProfiles, .dtwNormDist, .AItraces, .AIxNorm,
%       .peakLeft, .peakRight, .eventInfo
%   .globalCorr - .meanMatrix, .keypointNames, .perEventMatrices
%
% NOTES / ASSUMPTIONS:
%   - Segment times in the excel file are assumed to be in SECONDS.
%   - Baseline = mean x,y position of the keypoint over BaselineDuration
%     seconds immediately preceding REM onset. This SAME baseline is used
%     to compute displacement during both the REM event and the
%     following Awake window, so a REM-vs-Awake displacement comparison
%     is relative to a common reference position.
%   - Within a single event, L and R keypoint traces are naturally the
%     same length (same frame range), so pairwise metrics use full
%     native resolution -- no resampling anywhere.
%   - Events with insufficient baseline data, an Awake window that runs
%     past the end of the CSV data, a missing/ambiguous CSV match, or too
%     few frames for a given metric are skipped for that metric with a
%     warning, not silently dropped.

arguments
    excelPath (1,:) char
    csvFolder (1,:) char
    fps (1,1) double {mustBePositive}
    options.Sheet (1,1) double = 2
    options.BaselineDuration (1,1) double = 10
    options.AwakeDuration (1,1) double = 30
    options.MaxLagSec (1,1) double = 2
    options.CoherenceNFFT (1,1) double = 256
    options.OutputPath (1,:) char = ''
    options.MakePlots (1,1) logical = true
    options.LineAlpha (1,1) double = 0.15
    options.PairwiseAnalysis (1,1) logical = true
    options.GlobalCorrAnalysis (1,1) logical = true
end

haveDTW = (exist('dtw', 'file') == 2);
if options.PairwiseAnalysis && ~haveDTW
    warning('dtw() not found (Signal Processing Toolbox) -- DTW distance will be skipped.');
end
haveLinkage = (exist('linkage', 'file') == 2) && (exist('dendrogram', 'file') == 2);
if options.GlobalCorrAnalysis && ~haveLinkage
    warning('linkage()/dendrogram() not found (Statistics and Machine Learning Toolbox) -- will show correlation heatmaps only, no dendrograms.');
end

remColor   = [0.85 0.33 0.10];   % orange/red, matches prior REM convention
awakeColor = [0.30 0.55 0.90];   % blue, matches prior Awake convention

%% ---- 1) Read REM excel file ----
T = readtable(excelPath, 'Sheet', options.Sheet, 'ReadVariableNames', false);
if strcmp(T{1,1}{1}, 'file name')
    T = T(2:end,:);
end
nRec = height(T);

%% ---- 2) Fixed axes for pairwise metrics (same across all events/pairs) ----
fullMaxLagSamp = max(1, round(options.MaxLagSec * fps));
lagsSecFixed   = (-fullMaxLagSamp:fullMaxLagSamp)' / fps;
nfftCoh        = options.CoherenceNFFT;
fCohFixed      = (0:floor(nfftCoh/2))' * (fps / nfftCoh);

%% ---- 3) Containers (separate for REM and Awake) ----
kpNames  = {};
pairsList = struct('pairName', {}, 'leftName', {}, 'leftIdx', {}, 'rightName', {}, 'rightIdx', {});

kpDataREM   = struct('rawTraces', {}, 'xNorm', {}, 'eventInfo', {});
kpDataAwake = struct('rawTraces', {}, 'xNorm', {}, 'eventInfo', {});

pairDataREM   = struct('zeroLagCorr', {}, 'xcorrProfiles', {}, 'coherenceProfiles', {}, ...
    'dtwNormDist', {}, 'AItraces', {}, 'AIxNorm', {}, 'peakLeft', {}, 'peakRight', {}, 'eventInfo', {});
pairDataAwake = pairDataREM;

globalCorrMatsREM   = {};
globalCorrMatsAwake = {};

%% ---- 4) Loop over recordings / REM segments ----
for r = 1:nRec
    tdmsPath = T{r,1}{1};
    if isempty(tdmsPath) || (isstring(tdmsPath) && strlength(tdmsPath) == 0)
        continue
    end

    segRaw = table2array(T(r, 2:end));
    segRaw = segRaw(~isnan(segRaw));
    if isempty(segRaw)
        continue
    end
    segments = reshape(segRaw, 2, [])';   % nSeg x 2 : [start end] in seconds

    [~, baseName, ~] = fileparts(tdmsPath);
    csvMatches = dir(fullfile(csvFolder, [baseName '*.csv']));
    if isempty(csvMatches)
        warning('No CSV found for recording "%s" (row %d). Skipping.', baseName, r);
        continue
    elseif numel(csvMatches) > 1
        warning('Multiple CSVs match "%s" (row %d): using "%s". Rename files if this is wrong.', ...
            baseName, r, csvMatches(1).name);
    end
    csvFile = fullfile(csvMatches(1).folder, csvMatches(1).name);

    try
        kp = loadDLCcsv(csvFile);
    catch ME
        warning('Failed to load CSV "%s" (row %d): %s. Skipping.', csvFile, r, ME.message);
        continue
    end

    if isempty(kpNames)
        kpNames = {kp.name};
        for k = 1:numel(kpNames)
            kpDataREM(k).rawTraces = {};   kpDataREM(k).xNorm = {};
            kpDataREM(k).eventInfo = struct('recording',{}, 'rowIndex',{}, 'segIndex',{}, ...
                'startTime',{}, 'endTime',{}, 'durationSec',{}, 'nFramesBaseline',{});
            kpDataAwake(k).rawTraces = {}; kpDataAwake(k).xNorm = {};
            kpDataAwake(k).eventInfo = struct('recording',{}, 'rowIndex',{}, 'segIndex',{}, ...
                'startTime',{}, 'endTime',{}, 'durationSec',{}, 'nFramesBaseline',{});
        end

        pairsList = findSymmetricPairs(kpNames);
        for p = 1:numel(pairsList)
            pairDataREM(p) = emptyPairEntry();
            pairDataAwake(p) = emptyPairEntry();
        end
        if isempty(pairsList)
            warning('No left_X / right_X keypoint pairs detected in "%s" -- pairwise L/R analysis will be empty.', csvFile);
        end
    elseif ~isequal(kpNames, {kp.name})
        warning('Keypoint set in "%s" differs from earlier CSVs. Skipping recording.', csvFile);
        continue
    end

    nFramesTotal = numel(kp(1).frameIdx);

    for s = 1:size(segments,1)
        remStart = segments(s,1);
        remEnd   = segments(s,2);
        if isnan(remStart) || isnan(remEnd) || remEnd <= remStart
            continue
        end

        remStartFrame = round(remStart * fps);
        remEndFrame   = round(remEnd   * fps);
        baseStartFrame = remStartFrame - round(options.BaselineDuration * fps);

        baseStartRow  = baseStartFrame + 1;
        remStartRow   = remStartFrame + 1;
        remEndRow     = remEndFrame + 1;

        if baseStartRow < 1
            warning('Row %d, segment %d: insufficient baseline (event too close to recording start). Skipping.', r, s);
            continue
        end
        if remEndRow > nFramesTotal
            warning('Row %d, segment %d: REM event extends past end of CSV data. Truncating.', r, s);
            remEndRow = nFramesTotal;
        end
        if remStartRow > remEndRow
            warning('Row %d, segment %d: no valid frames in REM window. Skipping.', r, s);
            continue
        end

        baselineRows = baseStartRow:(remStartRow-1);
        remRows      = remStartRow:remEndRow;

        eiREM.recording       = csvMatches(1).name;
        eiREM.rowIndex        = r;
        eiREM.segIndex        = s;
        eiREM.startTime       = remStart;
        eiREM.endTime         = remEnd;
        eiREM.durationSec     = remEnd - remStart;
        eiREM.nFramesBaseline = numel(baselineRows);

        [kpDataREM, eventDispREM] = computeAndStoreKeypointDisplacement(kp, kpDataREM, baselineRows, remRows, eiREM);
        [pairDataREM, globalCorrMatsREM] = accumulatePairwiseAndGlobal(pairDataREM, globalCorrMatsREM, ...
            eventDispREM, pairsList, eiREM, fullMaxLagSamp, nfftCoh, fCohFixed, fps, haveDTW, options);

        % ---- Awake window: AwakeDuration seconds following REM offset, same baseline ----
        awakeStart = remEnd;
        awakeEnd   = remEnd + options.AwakeDuration;
        awakeStartFrame = round(awakeStart * fps);
        awakeEndFrame   = round(awakeEnd   * fps);
        awakeStartRow = awakeStartFrame + 1;
        awakeEndRow   = awakeEndFrame + 1;

        if awakeStartRow > nFramesTotal
            warning('Row %d, segment %d: Awake window starts past end of CSV data. Skipping Awake for this event.', r, s);
            continue
        end
        if awakeEndRow > nFramesTotal
            warning('Row %d, segment %d: Awake window extends past end of CSV data. Truncating.', r, s);
            awakeEndRow = nFramesTotal;
        end
        if awakeStartRow > awakeEndRow
            warning('Row %d, segment %d: no valid frames in Awake window. Skipping.', r, s);
            continue
        end
        awakeRows = awakeStartRow:awakeEndRow;

        eiAwake.recording       = csvMatches(1).name;
        eiAwake.rowIndex        = r;
        eiAwake.segIndex        = s;
        eiAwake.startTime       = awakeStart;
        eiAwake.endTime         = awakeEnd;
        eiAwake.durationSec     = awakeEnd - awakeStart;
        eiAwake.nFramesBaseline = numel(baselineRows);

        [kpDataAwake, eventDispAwake] = computeAndStoreKeypointDisplacement(kp, kpDataAwake, baselineRows, awakeRows, eiAwake);
        [pairDataAwake, globalCorrMatsAwake] = accumulatePairwiseAndGlobal(pairDataAwake, globalCorrMatsAwake, ...
            eventDispAwake, pairsList, eiAwake, fullMaxLagSamp, nfftCoh, fCohFixed, fps, haveDTW, options);
    end
end

%% ---- 5) Package results ----
results.REM   = packageCondition(kpNames, kpDataREM, pairsList, pairDataREM, globalCorrMatsREM, lagsSecFixed, fCohFixed);
results.Awake = packageCondition(kpNames, kpDataAwake, pairsList, pairDataAwake, globalCorrMatsAwake, lagsSecFixed, fCohFixed);

if ~isempty(options.OutputPath)
    save(options.OutputPath, 'results');
end

%% ---- 6) Plot ----
if options.MakePlots
    plotKeypointTraces(results.REM.keypoints, options.LineAlpha, 'REM', remColor);
    plotKeypointTraces(results.Awake.keypoints, options.LineAlpha, 'Awake', awakeColor);

    if options.PairwiseAnalysis
        plotAItraces(results.REM.pairwise, options.LineAlpha, 'REM', remColor);
        plotAItraces(results.Awake.pairwise, options.LineAlpha, 'Awake', awakeColor);
        plotPairwiseComparison(results.REM.pairwise, results.Awake.pairwise, remColor, awakeColor);
    end

    if options.GlobalCorrAnalysis
        plotGlobalCorrComparison(results.REM.globalCorr, results.Awake.globalCorr, haveLinkage);
    end
end

end % analyzeREMKeypointDisplacement


%% =====================================================================
function pe = emptyPairEntry()
pe.zeroLagCorr       = [];
pe.xcorrProfiles     = [];
pe.coherenceProfiles = [];
pe.dtwNormDist       = [];
pe.AItraces          = {};
pe.AIxNorm           = {};
pe.peakLeft          = [];
pe.peakRight         = [];
pe.eventInfo         = struct('recording',{}, 'rowIndex',{}, 'segIndex',{}, ...
    'startTime',{}, 'endTime',{}, 'durationSec',{});
end % emptyPairEntry


%% =====================================================================
function cond = packageCondition(kpNames, kpData, pairsList, pairData, globalCorrMats, lagsSecFixed, fCohFixed)
cond.keypoints = struct('name', {}, 'rawTraces', {}, 'xNorm', {}, 'eventInfo', {});
for k = 1:numel(kpNames)
    cond.keypoints(k).name      = kpNames{k};
    cond.keypoints(k).rawTraces = kpData(k).rawTraces;
    cond.keypoints(k).xNorm     = kpData(k).xNorm;
    cond.keypoints(k).eventInfo = kpData(k).eventInfo;
end

cond.pairwise = struct('pairName', {}, 'leftName', {}, 'rightName', {}, ...
    'zeroLagCorr', {}, 'lagsSec', {}, 'xcorrProfiles', {}, ...
    'coherenceFreq', {}, 'coherenceProfiles', {}, 'dtwNormDist', {}, ...
    'AItraces', {}, 'AIxNorm', {}, 'peakLeft', {}, 'peakRight', {}, 'eventInfo', {});
for p = 1:numel(pairsList)
    cond.pairwise(p).pairName          = pairsList(p).pairName;
    cond.pairwise(p).leftName          = pairsList(p).leftName;
    cond.pairwise(p).rightName         = pairsList(p).rightName;
    cond.pairwise(p).zeroLagCorr       = pairData(p).zeroLagCorr;
    cond.pairwise(p).lagsSec           = lagsSecFixed;
    cond.pairwise(p).xcorrProfiles     = pairData(p).xcorrProfiles;
    cond.pairwise(p).coherenceFreq     = fCohFixed;
    cond.pairwise(p).coherenceProfiles = pairData(p).coherenceProfiles;
    cond.pairwise(p).dtwNormDist       = pairData(p).dtwNormDist;
    cond.pairwise(p).AItraces          = pairData(p).AItraces;
    cond.pairwise(p).AIxNorm           = pairData(p).AIxNorm;
    cond.pairwise(p).peakLeft          = pairData(p).peakLeft;
    cond.pairwise(p).peakRight         = pairData(p).peakRight;
    cond.pairwise(p).eventInfo         = pairData(p).eventInfo;
end

if ~isempty(globalCorrMats)
    meanCorrMat = mean(cat(3, globalCorrMats{:}), 3, 'omitnan');
else
    meanCorrMat = [];
end
cond.globalCorr.meanMatrix       = meanCorrMat;
cond.globalCorr.keypointNames    = kpNames;
cond.globalCorr.perEventMatrices = globalCorrMats;
end % packageCondition


%% =====================================================================
function [kpData, eventDispRawTemp] = computeAndStoreKeypointDisplacement(kp, kpData, baselineRows, eventRows, ei)
eventDispRawTemp = cell(1, numel(kp));
for k = 1:numel(kp)
    x0 = mean(kp(k).x(baselineRows), 'omitnan');
    y0 = mean(kp(k).y(baselineRows), 'omitnan');

    dx = kp(k).x(eventRows) - x0;
    dy = kp(k).y(eventRows) - y0;
    dispTrace = sqrt(dx.^2 + dy.^2);

    nRaw = numel(dispTrace);
    if nRaw < 2
        continue
    end
    eventDispRawTemp{k} = dispTrace(:);

    xNormThisEvent = linspace(0, 100, nRaw);
    kpData(k).rawTraces{end+1} = dispTrace(:);
    kpData(k).xNorm{end+1}     = xNormThisEvent(:);
    kpData(k).eventInfo(end+1) = ei;
end
end % computeAndStoreKeypointDisplacement


%% =====================================================================
function [pairData, globalCorrMats] = accumulatePairwiseAndGlobal(pairData, globalCorrMats, ...
    eventDispRawTemp, pairsList, ei, fullMaxLagSamp, nfftCoh, fCohFixed, fps, haveDTW, options)

if options.GlobalCorrAnalysis
    lens = cellfun(@numel, eventDispRawTemp);
    if all(lens >= 5)
        M = cat(2, eventDispRawTemp{:});
        try
            Cmat = corrcoef(M);
            globalCorrMats{end+1} = Cmat; %#ok<AGROW>
        catch ME
            warning('Row %d, segment %d: correlation matrix failed (%s). Skipping.', ei.rowIndex, ei.segIndex, ME.message);
        end
    end
end

if options.PairwiseAnalysis
    for p = 1:numel(pairsList)
        L = eventDispRawTemp{pairsList(p).leftIdx};
        R = eventDispRawTemp{pairsList(p).rightIdx};
        if isempty(L) || isempty(R) || numel(L) < 5 || numel(R) < 5
            continue
        end
        n = min(numel(L), numel(R));
        L = L(1:n); R = R(1:n);

        cc = corrcoef(L, R);
        zeroLagVal = cc(1,2);

        if n - 1 >= fullMaxLagSamp
            [xc, ~] = xcorr(L - mean(L), R - mean(R), fullMaxLagSamp, 'coeff');
            xcVec = xc(:);
        else
            xcVec = nan(2*fullMaxLagSamp+1, 1);
        end

        try
            [Cxy, ~] = mscohere(L, R, [], [], nfftCoh, fps);
            CxyVec = Cxy(:);
            if numel(CxyVec) ~= numel(fCohFixed)
                CxyVec = nan(size(fCohFixed));
            end
        catch
            CxyVec = nan(size(fCohFixed));
        end

        if haveDTW
            try
                dtwDist = dtw(L, R) / n;
            catch
                dtwDist = NaN;
            end
        else
            dtwDist = NaN;
        end

        Lz = zscoreSafe(L);
        Rz = zscoreSafe(R);
        AItrace = abs(Lz - Rz);
        xNormAI = linspace(0, 100, n)';

        pairData(p).zeroLagCorr(end+1,1)       = zeroLagVal;
        pairData(p).xcorrProfiles(:,end+1)     = xcVec;
        pairData(p).coherenceProfiles(:,end+1) = CxyVec;
        pairData(p).dtwNormDist(end+1,1)       = dtwDist;
        pairData(p).AItraces{end+1}            = AItrace;
        pairData(p).AIxNorm{end+1}             = xNormAI;
        pairData(p).peakLeft(end+1,1)          = max(L);
        pairData(p).peakRight(end+1,1)         = max(R);

        eiPair.recording   = ei.recording;
        eiPair.rowIndex    = ei.rowIndex;
        eiPair.segIndex    = ei.segIndex;
        eiPair.startTime   = ei.startTime;
        eiPair.endTime     = ei.endTime;
        eiPair.durationSec = ei.durationSec;
        pairData(p).eventInfo(end+1) = eiPair;
    end
end
end % accumulatePairwiseAndGlobal


%% =====================================================================
function kp = loadDLCcsv(csvFile)
% Loads a DeepLabCut-style csv with 3 header rows:
%   row1: scorer / row2: bodyparts (x,y,likelihood x3 per keypoint) / row3: coords
% Returns kp: struct array, one per keypoint, fields .name, .frameIdx, .x, .y, .likelihood

fid = fopen(csvFile, 'r');
if fid == -1
    error('Could not open file: %s', csvFile);
end
line1 = fgetl(fid); %#ok<NASGU>
line2 = fgetl(fid);
line3 = fgetl(fid);
fclose(fid);

if numel(strfind(line2, ',')) >= numel(strfind(line2, sprintf('\t')))
    delim = ',';
else
    delim = sprintf('\t');
end

bodypartsCells = strsplit(line2, delim);
bodypartsCells(1) = [];

if mod(numel(bodypartsCells), 3) ~= 0
    error('Unexpected header format in %s: bodyparts columns not a multiple of 3.', csvFile);
end

keypointNames = bodypartsCells(1:3:end);

data = readmatrix(csvFile, 'NumHeaderLines', 3, 'Delimiter', delim);
frameIdx = data(:,1);
coordData = data(:,2:end);

nKp = numel(keypointNames);
kp = struct('name', cell(1,nKp), 'frameIdx', cell(1,nKp), 'x', cell(1,nKp), ...
    'y', cell(1,nKp), 'likelihood', cell(1,nKp));
for k = 1:nKp
    colBase = (k-1)*3;
    kp(k).name       = keypointNames{k};
    kp(k).frameIdx   = frameIdx;
    kp(k).x          = coordData(:, colBase+1);
    kp(k).y          = coordData(:, colBase+2);
    kp(k).likelihood = coordData(:, colBase+3);
end

end % loadDLCcsv


%% =====================================================================
function pairs = findSymmetricPairs(kpNames)
% Detects "left_X" / "right_X" keypoint pairs by name.
pairs = struct('pairName', {}, 'leftName', {}, 'leftIdx', {}, 'rightName', {}, 'rightIdx', {});
for i = 1:numel(kpNames)
    nm = kpNames{i};
    if startsWith(nm, 'left_')
        suffix = nm(6:end);
        rightNm = ['right_' suffix];
        j = find(strcmp(kpNames, rightNm), 1);
        if ~isempty(j)
            pairs(end+1).pairName  = suffix; %#ok<AGROW>
            pairs(end).leftName    = nm;
            pairs(end).leftIdx     = i;
            pairs(end).rightName   = rightNm;
            pairs(end).rightIdx    = j;
        end
    end
end
end % findSymmetricPairs


%% =====================================================================
function z = zscoreSafe(x)
s = std(x);
if s > 0
    z = (x - mean(x)) / s;
else
    z = x - mean(x);
end
end % zscoreSafe


%% =====================================================================
function applyDarkAxes(ax)
bgColor   = [0.12 0.12 0.12];
textColor = [0.92 0.92 0.92];
gridColor = [0.4 0.4 0.4];
set(ancestor(ax,'figure'), 'Color', bgColor);
set(ax, 'Color', bgColor, 'XColor', textColor, 'YColor', textColor, ...
    'GridColor', gridColor, 'Box', 'on');
end % applyDarkAxes


%% =====================================================================
function plotKeypointTraces(kpResults, lineAlpha, condLabel, lineColor)
% One figure per keypoint. All events overlaid, semi-transparent, each
% event plotted at its own native sample count on its own 0-100% x-axis.
textColor = [0.92 0.92 0.92];

for k = 1:numel(kpResults)
    traces = kpResults(k).rawTraces;
    xNormAll = kpResults(k).xNorm;
    if isempty(traces)
        continue
    end

    fig = figure('Name', sprintf('%s_%s', condLabel, kpResults(k).name));
    ax = axes(fig);
    applyDarkAxes(ax);
    hold(ax, 'on');

    for e = 1:numel(traces)
        h = plot(ax, xNormAll{e}, traces{e}, 'Color', lineColor);
        h.Color(4) = lineAlpha;
    end

    xlabel(ax, sprintf('%% of %s window completed', condLabel), 'Color', textColor);
    ylabel(ax, 'Displacement from baseline (pixels)', 'Color', textColor);
    title(ax, sprintf('%s: %s (n = %d events)', condLabel, kpResults(k).name, numel(traces)), ...
        'Interpreter', 'none', 'Color', textColor);
    xlim(ax, [0 100]);
    grid(ax, 'on');
    hold(ax, 'off');
end
end % plotKeypointTraces


%% =====================================================================
function plotAItraces(pairResults, lineAlpha, condLabel, lineColor)
% One figure per L/R pair: asymmetry-index |zscore(L)-zscore(R)| traces,
% overlaid semi-transparent, full resolution, each on its own 0-100% axis.
textColor = [0.92 0.92 0.92];

for p = 1:numel(pairResults)
    traces = pairResults(p).AItraces;
    xNormAll = pairResults(p).AIxNorm;
    if isempty(traces)
        continue
    end

    fig = figure('Name', sprintf('AI_%s_%s', condLabel, pairResults(p).pairName));
    ax = axes(fig);
    applyDarkAxes(ax);
    hold(ax, 'on');

    for e = 1:numel(traces)
        h = plot(ax, xNormAll{e}, traces{e}, 'Color', lineColor);
        h.Color(4) = lineAlpha;
    end

    xlabel(ax, sprintf('%% of %s window completed', condLabel), 'Color', textColor);
    ylabel(ax, '|z(left) - z(right)|  (asymmetry index)', 'Color', textColor);
    title(ax, sprintf('%s: %s L/R asymmetry (n = %d events)', condLabel, pairResults(p).pairName, numel(traces)), ...
        'Interpreter', 'none', 'Color', textColor);
    xlim(ax, [0 100]);
    grid(ax, 'on');
    hold(ax, 'off');
end
end % plotAItraces


%% =====================================================================
function jitterScatter(ax, x0, data, color)
if isempty(data)
    return
end
jitter = (rand(size(data)) - 0.5) * 0.3;
scatter(ax, x0 + jitter, data, 25, color, 'filled', 'MarkerFaceAlpha', 0.5);
end % jitterScatter


%% =====================================================================
function plotPairwiseComparison(remPairs, awakePairs, remColor, awakeColor)
% REM-vs-Awake comparison of L/R symmetry metrics across pairs:
%   1) zero-lag correlation + DTW distance (jittered, REM vs Awake side by side)
%   2) mean cross-correlation lag profile, REM vs Awake overlay, per pair
%   3) mean coherence spectrum, REM vs Awake overlay, per pair
textColor = [0.92 0.92 0.92];
nPairs = numel(remPairs);
if nPairs == 0
    return
end
pairNames = {remPairs.pairName};

%% Figure 1: zero-lag correlation + DTW distance, REM vs Awake
fig = figure('Name', 'LR_REMvsAwake_zeroLagCorr_and_DTW');
t = tiledlayout(fig, 1, 2);

ax1 = nexttile(t); applyDarkAxes(ax1); hold(ax1, 'on');
for p = 1:nPairs
    vR = remPairs(p).zeroLagCorr;   vR = vR(isfinite(vR));
    vA = awakePairs(p).zeroLagCorr; vA = vA(isfinite(vA));
    jitterScatter(ax1, p-0.15, vR, remColor);
    jitterScatter(ax1, p+0.15, vA, awakeColor);
    if ~isempty(vR)
        errorbar(ax1, p-0.15, mean(vR), std(vR)/sqrt(numel(vR)), 'o', 'Color', remColor, 'MarkerFaceColor', remColor, 'LineWidth', 1.5, 'CapSize', 5);
    end
    if ~isempty(vA)
        errorbar(ax1, p+0.15, mean(vA), std(vA)/sqrt(numel(vA)), 'o', 'Color', awakeColor, 'MarkerFaceColor', awakeColor, 'LineWidth', 1.5, 'CapSize', 5);
    end
end
xlim(ax1, [0.5 nPairs+0.5]); xticks(ax1, 1:nPairs); xticklabels(ax1, strrep(pairNames,'_',' '));
ylabel(ax1, 'Zero-lag correlation (r)', 'Color', textColor);
title(ax1, 'L vs R zero-lag correlation: REM vs Awake', 'Color', textColor);
yline(ax1, 0, 'Color', [0.6 0.6 0.6]);
grid(ax1, 'on'); hold(ax1, 'off');

ax2 = nexttile(t); applyDarkAxes(ax2); hold(ax2, 'on');
for p = 1:nPairs
    vR = remPairs(p).dtwNormDist;   vR = vR(isfinite(vR));
    vA = awakePairs(p).dtwNormDist; vA = vA(isfinite(vA));
    jitterScatter(ax2, p-0.15, vR, remColor);
    jitterScatter(ax2, p+0.15, vA, awakeColor);
    if ~isempty(vR)
        errorbar(ax2, p-0.15, mean(vR), std(vR)/sqrt(numel(vR)), 'o', 'Color', remColor, 'MarkerFaceColor', remColor, 'LineWidth', 1.5, 'CapSize', 5);
    end
    if ~isempty(vA)
        errorbar(ax2, p+0.15, mean(vA), std(vA)/sqrt(numel(vA)), 'o', 'Color', awakeColor, 'MarkerFaceColor', awakeColor, 'LineWidth', 1.5, 'CapSize', 5);
    end
end
xlim(ax2, [0.5 nPairs+0.5]); xticks(ax2, 1:nPairs); xticklabels(ax2, strrep(pairNames,'_',' '));
ylabel(ax2, 'DTW distance / mean(n samples)', 'Color', textColor);
title(ax2, 'L vs R shape dissimilarity: REM vs Awake', 'Color', textColor);
grid(ax2, 'on'); hold(ax2, 'off');

legendDummy(ax1, remColor, awakeColor);

%% Figure 2: cross-correlation lag profile, REM vs Awake overlay
fig = figure('Name', 'LR_REMvsAwake_xcorr');
t = tiledlayout(fig, 1, nPairs);
for p = 1:nPairs
    ax = nexttile(t); applyDarkAxes(ax); hold(ax, 'on');
    plotMeanSEMLine(ax, remPairs(p).lagsSec, remPairs(p).xcorrProfiles, remColor);
    plotMeanSEMLine(ax, awakePairs(p).lagsSec, awakePairs(p).xcorrProfiles, awakeColor);
    xline(ax, 0, 'Color', [0.6 0.6 0.6]);
    yline(ax, 0, 'Color', [0.4 0.4 0.4]);
    xlabel(ax, 'Lag (s)', 'Color', textColor);
    ylabel(ax, 'Cross-corr (coeff)', 'Color', textColor);
    title(ax, pairNames{p}, 'Interpreter', 'none', 'Color', textColor);
    grid(ax, 'on'); hold(ax, 'off');
end
sgtitle(t, 'L vs R cross-correlation: REM (orange/red) vs Awake (blue), mean \pm SEM', 'Color', textColor);

%% Figure 3: coherence spectrum, REM vs Awake overlay
fig = figure('Name', 'LR_REMvsAwake_coherence');
t = tiledlayout(fig, 1, nPairs);
for p = 1:nPairs
    ax = nexttile(t); applyDarkAxes(ax); hold(ax, 'on');
    plotMeanSEMLine(ax, remPairs(p).coherenceFreq, remPairs(p).coherenceProfiles, remColor);
    plotMeanSEMLine(ax, awakePairs(p).coherenceFreq, awakePairs(p).coherenceProfiles, awakeColor);
    ylim(ax, [0 1]);
    xlabel(ax, 'Frequency (Hz)', 'Color', textColor);
    ylabel(ax, 'Coherence', 'Color', textColor);
    title(ax, pairNames{p}, 'Interpreter', 'none', 'Color', textColor);
    grid(ax, 'on'); hold(ax, 'off');
end
sgtitle(t, 'L vs R coherence: REM (orange/red) vs Awake (blue), mean \pm SEM', 'Color', textColor);

end % plotPairwiseComparison


%% =====================================================================
function plotMeanSEMLine(ax, x, M, color)
if isempty(M)
    return
end
mu = mean(M, 2, 'omitnan');
n  = sum(isfinite(M), 2);
sem = std(M, 0, 2, 'omitnan') ./ max(sqrt(n), 1);
valid = isfinite(mu);
fill(ax, [x(valid); flipud(x(valid))], [mu(valid)-sem(valid); flipud(mu(valid)+sem(valid))], ...
    color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(ax, x, mu, 'Color', color, 'LineWidth', 2);
end % plotMeanSEMLine


%% =====================================================================
function legendDummy(ax, remColor, awakeColor)
hold(ax, 'on');
h1 = plot(ax, nan, nan, 'o', 'Color', remColor, 'MarkerFaceColor', remColor, 'DisplayName', 'REM');
h2 = plot(ax, nan, nan, 'o', 'Color', awakeColor, 'MarkerFaceColor', awakeColor, 'DisplayName', 'Awake');
lg = legend(ax, [h1 h2], 'Location', 'best');
lg.TextColor = [0.92 0.92 0.92];
lg.Color = [0.12 0.12 0.12];
end % legendDummy


%% =====================================================================
function plotGlobalCorrComparison(remGlobal, awakeGlobal, haveLinkage)
% Cross-keypoint mean correlation matrices for REM and Awake, plus their
% difference, and (if available) dendrograms for each.
textColor = [0.92 0.92 0.92];
names = strrep(remGlobal.keypointNames, '_', ' ');
n = numel(names);

hasREM = ~isempty(remGlobal.meanMatrix);
hasAwake = ~isempty(awakeGlobal.meanMatrix);
if ~hasREM && ~hasAwake
    return
end

fig = figure('Name', 'Global_keypoint_corr_REM_vs_Awake');
set(fig, 'Color', [0.12 0.12 0.12]);
t = tiledlayout(fig, 1, 3);

ax1 = nexttile(t); applyDarkAxes(ax1);
if hasREM
    imagesc(ax1, remGlobal.meanMatrix); axis(ax1, 'square'); clim(ax1, [0 1]); colorbar(ax1); colormap(ax1, 'turbo');
    xticks(ax1, 1:n); yticks(ax1, 1:n); xticklabels(ax1, names); yticklabels(ax1, names); xtickangle(ax1, 45);
end
title(ax1, 'REM: mean cross-keypoint correlation', 'Color', textColor);

ax2 = nexttile(t); applyDarkAxes(ax2);
if hasAwake
    imagesc(ax2, awakeGlobal.meanMatrix); axis(ax2, 'square'); clim(ax2, [0 1]); colorbar(ax2); colormap(ax2, 'turbo');
    xticks(ax2, 1:n); yticks(ax2, 1:n); xticklabels(ax2, names); yticklabels(ax2, names); xtickangle(ax2, 45);
end
title(ax2, 'Awake: mean cross-keypoint correlation', 'Color', textColor);

ax3 = nexttile(t); applyDarkAxes(ax3);
if hasREM && hasAwake
    diffMat = remGlobal.meanMatrix - awakeGlobal.meanMatrix;
    imagesc(ax3, diffMat); axis(ax3, 'square'); clim(ax3, [-.25 .25]); colorbar(ax3); colormap(ax3, 'turbo');
    xticks(ax3, 1:n); yticks(ax3, 1:n); xticklabels(ax3, names); yticklabels(ax3, names); xtickangle(ax3, 45);
end
title(ax3, 'REM - Awake', 'Color', textColor);

sgtitle(t, 'Cross-keypoint correlation structure: REM vs Awake', 'Color', textColor);

if haveLinkage
    fig2 = figure('Name', 'Global_keypoint_dendrogram_REM_vs_Awake');
    set(fig2, 'Color', [0.12 0.12 0.12]);
    t2 = tiledlayout(fig2, 1, 2);

    axD1 = nexttile(t2); applyDarkAxes(axD1);
    if hasREM
        plotDendrogramOn(axD1, remGlobal.meanMatrix, names);
    end
    title(axD1, 'REM clustering (1 - correlation)', 'Color', textColor);

    axD2 = nexttile(t2); applyDarkAxes(axD2);
    if hasAwake
        plotDendrogramOn(axD2, awakeGlobal.meanMatrix, names);
    end
    title(axD2, 'Awake clustering (1 - correlation)', 'Color', textColor);

    sgtitle(t2, 'Hierarchical clustering of keypoints: REM vs Awake', 'Color', textColor);
end
end % plotGlobalCorrComparison


%% =====================================================================
function plotDendrogramOn(ax, M, names)
n = size(M,1);
D = 1 - M;
D = (D + D') / 2;
D(1:n+1:end) = 0;
Z = linkage(squareform(D, 'tovector'), 'average');
dendrogram(ax, Z, 0, 'Labels', names);
set(ax, 'XTickLabelRotation', 45);
end % plotDendrogramOn

% function results = analyzeREMKeypointDisplacement(excelPath, csvFolder, fps, options)
% % ANALYZEREMKEYPOINTDISPLACEMENT
% % For each REM event listed in an excel file, finds the matching DLC
% % keypoint CSV, computes a 10-second-prior baseline (x,y) position for
% % each keypoint, then computes the frame-by-frame displacement magnitude
% % from that baseline across the REM event.
% %
% % In addition to the per-keypoint displacement traces, this also computes
% % left-vs-right symmetry/dissimilarity metrics for any keypoint pairs
% % named "left_X" / "right_X":
% %   - zero-lag correlation + full cross-correlation lag profile
% %   - magnitude-squared coherence spectrum
% %   - DTW (dynamic time warping) distance, duration-normalized
% %   - an asymmetry-index trace: |zscore(left) - zscore(right)| over the
% %     event, at full resolution (no downsampling)
% %   - paired peak-displacement values (for L vs R scatter)
% % and a cross-keypoint correlation matrix + dendrogram across ALL
% % keypoints (not just L/R pairs), to compare against the higher-
% % dimensional/disjointed structure your ROI-based PCA/dendrogram analysis
% % already showed.
% %
% % USAGE:
% %   results = analyzeREMKeypointDisplacement(excelPath, csvFolder, fps)
% %   results = analyzeREMKeypointDisplacement(excelPath, csvFolder, fps, ...
% %                 'Sheet', 2, 'BaselineDuration', 10, ...
% %                 'MaxLagSec', 2, 'OutputPath', 'results.mat', 'MakePlots', true)
% %
% % INPUTS:
% %   excelPath  - path to REM excel file (Sheet 2, col 1 = tdms/recording
% %                path, cols 2:end = pairs of [start end] REM segment
% %                times in seconds, NaN-padded)
% %   csvFolder  - folder containing DLC keypoint CSVs. CSV filenames are
% %                expected to START WITH the same base filename as column 1
% %                of the excel file (extension/path stripped).
% %   fps        - frame rate (frames/sec) of the video/keypoint tracking.
% %                REQUIRED because it is not encoded in either file.
% %
% % NAME-VALUE OPTIONS:
% %   'Sheet'             - excel sheet to read (default 2)
% %   'BaselineDuration'  - seconds prior to event onset used as baseline (default 10)
% %   'MaxLagSec'         - max lag (seconds) for the L/R cross-correlation profile (default 2)
% %   'CoherenceNFFT'     - FFT length for L/R coherence (default 256)
% %   'OutputPath'        - if provided, saves `results` struct to this .mat path
% %   'MakePlots'         - true/false, whether to generate plots (default true)
% %   'LineAlpha'         - transparency of each individual event trace (default 0.15)
% %   'PairwiseAnalysis'  - true/false, compute L/R pairwise metrics (default true)
% %   'GlobalCorrAnalysis'- true/false, compute cross-keypoint corr matrix + dendrogram (default true)
% %
% % OUTPUT: results struct with three fields:
% %   results.keypoints  - struct array, one per keypoint:
% %       .name, .rawTraces (cell, full-res displacement magnitude per event),
% %       .xNorm (cell, that event's own 0-100% x-axis), .eventInfo
% %   results.pairwise   - struct array, one per detected L/R keypoint pair:
% %       .pairName, .leftName, .rightName
% %       .zeroLagCorr      [nEvents x 1] Pearson r at zero lag
% %       .lagsSec          [nLags x 1] fixed lag axis (same for all events/pairs)
% %       .xcorrProfiles    [nLags x nEvents] (NaN column if event too short for MaxLagSec)
% %       .coherenceFreq    [nFreq x 1] fixed frequency axis
% %       .coherenceProfiles[nFreq x nEvents]
% %       .dtwNormDist      [nEvents x 1] DTW distance / mean(nSamplesL,nSamplesR)
% %       .AItraces         {1 x nEvents} full-res |zscore(L)-zscore(R)| traces
% %       .AIxNorm          {1 x nEvents} that event's own 0-100% x-axis
% %       .peakLeft, .peakRight [nEvents x 1] peak displacement magnitude per side
% %       .eventInfo        struct array, same event metadata as .keypoints
% %   results.globalCorr - cross-keypoint correlation structure:
% %       .meanMatrix       [nKp x nKp] mean correlation across events
% %       .keypointNames    {1 x nKp}
% %       .perEventMatrices {1 x nEvents} per-event [nKp x nKp] correlation matrices
% %
% % NOTES / ASSUMPTIONS:
% %   - Segment times in the excel file are assumed to be in SECONDS.
% %   - Baseline = mean x,y position of the keypoint over BaselineDuration
% %     seconds immediately preceding event onset.
% %   - Displacement is Euclidean magnitude from that baseline. Left/right
% %     SIGN (which direction each side moved) is not computed here -- only
% %     whether the two sides' displacement magnitude time-courses are
% %     correlated/coherent/similarly-shaped, which is what's needed to ask
% %     "do they twitch together or separately."
% %   - No resampling/downsampling anywhere: within one event, L and R
% %     keypoint traces are naturally the same length (same frame range),
% %     so pairwise metrics use the full native resolution. The global
% %     cross-keypoint correlation matrix is likewise computed per-event at
% %     native resolution (correlation is a summary statistic, not a
% %     time-resolved trace, so no common time grid is needed), then
% %     averaged across events.
% %   - Events with insufficient baseline data, a missing/ambiguous CSV
% %     match, or too few frames for a given metric are skipped for that
% %     metric with a warning, not silently dropped.
% 
% arguments
%     excelPath (1,:) char
%     csvFolder (1,:) char
%     fps (1,1) double {mustBePositive}
%     options.Sheet (1,1) double = 2
%     options.BaselineDuration (1,1) double = 10
%     options.MaxLagSec (1,1) double = 2
%     options.CoherenceNFFT (1,1) double = 256
%     options.OutputPath (1,:) char = ''
%     options.MakePlots (1,1) logical = true
%     options.LineAlpha (1,1) double = 0.15
%     options.PairwiseAnalysis (1,1) logical = true
%     options.GlobalCorrAnalysis (1,1) logical = true
% end
% 
% haveDTW = (exist('dtw', 'file') == 2);
% if options.PairwiseAnalysis && ~haveDTW
%     warning('dtw() not found (Signal Processing Toolbox) -- DTW distance will be skipped.');
% end
% haveLinkage = (exist('linkage', 'file') == 2) && (exist('dendrogram', 'file') == 2);
% if options.GlobalCorrAnalysis && ~haveLinkage
%     warning('linkage()/dendrogram() not found (Statistics and Machine Learning Toolbox) -- will show the correlation heatmap only, no dendrogram.');
% end
% 
% %% ---- 1) Read REM excel file ----
% T = readtable(excelPath, 'Sheet', options.Sheet, 'ReadVariableNames', false);
% if strcmp(T{1,1}{1}, 'file name')
%     T = T(2:end,:);
% end
% nRec = height(T);
% 
% %% ---- 2) Fixed axes for pairwise metrics (same across all events/pairs) ----
% fullMaxLagSamp = max(1, round(options.MaxLagSec * fps));
% lagsSecFixed   = (-fullMaxLagSamp:fullMaxLagSamp)' / fps;
% nfftCoh        = options.CoherenceNFFT;
% fCohFixed      = (0:floor(nfftCoh/2))' * (fps / nfftCoh);
% 
% %% ---- 3) Containers ----
% kpNames  = {};
% kpData   = struct('rawTraces', {}, 'xNorm', {}, 'eventInfo', {});
% pairsList = struct('pairName', {}, 'leftName', {}, 'leftIdx', {}, 'rightName', {}, 'rightIdx', {});
% pairData  = struct('zeroLagCorr', {}, 'xcorrProfiles', {}, 'coherenceProfiles', {}, ...
%     'dtwNormDist', {}, 'AItraces', {}, 'AIxNorm', {}, 'peakLeft', {}, 'peakRight', {}, 'eventInfo', {});
% globalCorrEventMatrices = {};
% 
% %% ---- 4) Loop over recordings / REM segments ----
% for r = 1:nRec
%     tdmsPath = T{r,1}{1};
%     if isempty(tdmsPath) || (isstring(tdmsPath) && strlength(tdmsPath) == 0)
%         continue
%     end
% 
%     segRaw = table2array(T(r, 2:end));
%     segRaw = segRaw(~isnan(segRaw));
%     if isempty(segRaw)
%         continue
%     end
%     segments = reshape(segRaw, 2, [])';   % nSeg x 2 : [start end] in seconds
% 
%     [~, baseName, ~] = fileparts(tdmsPath);
%     csvMatches = dir(fullfile(csvFolder, [baseName '*.csv']));
%     if isempty(csvMatches)
%         warning('No CSV found for recording "%s" (row %d). Skipping.', baseName, r);
%         continue
%     elseif numel(csvMatches) > 1
%         warning('Multiple CSVs match "%s" (row %d): using "%s". Rename files if this is wrong.', ...
%             baseName, r, csvMatches(1).name);
%     end
%     csvFile = fullfile(csvMatches(1).folder, csvMatches(1).name);
% 
%     try
%         kp = loadDLCcsv(csvFile);
%     catch ME
%         warning('Failed to load CSV "%s" (row %d): %s. Skipping.', csvFile, r, ME.message);
%         continue
%     end
% 
%     if isempty(kpNames)
%         kpNames = {kp.name};
%         for k = 1:numel(kpNames)
%             kpData(k).rawTraces = {};
%             kpData(k).xNorm     = {};
%             kpData(k).eventInfo = struct('recording',{}, 'rowIndex',{}, 'segIndex',{}, ...
%                 'startTime',{}, 'endTime',{}, 'durationSec',{}, 'nFramesBaseline',{});
%         end
% 
%         pairsList = findSymmetricPairs(kpNames);
%         for p = 1:numel(pairsList)
%             pairData(p).zeroLagCorr       = [];
%             pairData(p).xcorrProfiles     = [];
%             pairData(p).coherenceProfiles = [];
%             pairData(p).dtwNormDist       = [];
%             pairData(p).AItraces          = {};
%             pairData(p).AIxNorm           = {};
%             pairData(p).peakLeft          = [];
%             pairData(p).peakRight         = [];
%             pairData(p).eventInfo         = struct('recording',{}, 'rowIndex',{}, 'segIndex',{}, ...
%                 'startTime',{}, 'endTime',{}, 'durationSec',{});
%         end
%         if isempty(pairsList)
%             warning('No left_X / right_X keypoint pairs detected in "%s" -- pairwise L/R analysis will be empty.', csvFile);
%         end
%     elseif ~isequal(kpNames, {kp.name})
%         warning('Keypoint set in "%s" differs from earlier CSVs. Skipping recording.', csvFile);
%         continue
%     end
% 
%     nFramesTotal = numel(kp(1).frameIdx);
% 
%     for s = 1:size(segments,1)
%         startT = segments(s,1);
%         endT   = segments(s,2);
%         if isnan(startT) || isnan(endT) || endT <= startT
%             continue
%         end
% 
%         startFrame = round(startT * fps);
%         endFrame   = round(endT   * fps);
%         baseStartFrame = startFrame - round(options.BaselineDuration * fps);
% 
%         baseStartRow  = baseStartFrame + 1;
%         eventStartRow = startFrame + 1;
%         eventEndRow   = endFrame + 1;
% 
%         if baseStartRow < 1
%             warning('Row %d, segment %d: insufficient baseline (event too close to recording start). Skipping.', r, s);
%             continue
%         end
%         if eventEndRow > nFramesTotal
%             warning('Row %d, segment %d: event extends past end of CSV data. Truncating to available frames.', r, s);
%             eventEndRow = nFramesTotal;
%         end
%         if eventStartRow > eventEndRow
%             warning('Row %d, segment %d: no valid frames in event window. Skipping.', r, s);
%             continue
%         end
% 
%         baselineRows = baseStartRow:(eventStartRow-1);
%         eventRows    = eventStartRow:eventEndRow;
% 
%         ei.recording       = csvMatches(1).name;
%         ei.rowIndex        = r;
%         ei.segIndex        = s;
%         ei.startTime       = startT;
%         ei.endTime         = endT;
%         ei.durationSec     = endT - startT;
%         ei.nFramesBaseline = numel(baselineRows);
% 
%         % ---- Per-keypoint baseline + displacement (also cached for pairwise/global use) ----
%         eventDispRawTemp = cell(1, numel(kp));
%         for k = 1:numel(kp)
%             x0 = mean(kp(k).x(baselineRows), 'omitnan');
%             y0 = mean(kp(k).y(baselineRows), 'omitnan');
% 
%             dx = kp(k).x(eventRows) - x0;
%             dy = kp(k).y(eventRows) - y0;
%             dispTrace = sqrt(dx.^2 + dy.^2);
% 
%             nRaw = numel(dispTrace);
%             if nRaw < 2
%                 continue
%             end
%             eventDispRawTemp{k} = dispTrace(:);
% 
%             xNormThisEvent = linspace(0, 100, nRaw);
%             kpData(k).rawTraces{end+1} = dispTrace(:);
%             kpData(k).xNorm{end+1}     = xNormThisEvent(:);
%             kpData(k).eventInfo(end+1) = ei;
%         end
% 
%         % ---- Cross-keypoint correlation matrix for this event ----
%         if options.GlobalCorrAnalysis
%             lens = cellfun(@numel, eventDispRawTemp);
%             if all(lens >= 5)
%                 M = cat(2, eventDispRawTemp{:});
%                 try
%                     Cmat = corrcoef(M);
%                     globalCorrEventMatrices{end+1} = Cmat; %#ok<AGROW>
%                 catch ME
%                     warning('Row %d, segment %d: correlation matrix failed (%s). Skipping.', r, s, ME.message);
%                 end
%             end
%         end
% 
%         % ---- Pairwise (left vs right) metrics ----
%         if options.PairwiseAnalysis
%             for p = 1:numel(pairsList)
%                 L = eventDispRawTemp{pairsList(p).leftIdx};
%                 R = eventDispRawTemp{pairsList(p).rightIdx};
%                 if isempty(L) || isempty(R) || numel(L) < 5 || numel(R) < 5
%                     continue
%                 end
%                 n = min(numel(L), numel(R));
%                 L = L(1:n); R = R(1:n);
% 
%                 % zero-lag correlation (always computable if n>=2)
%                 cc = corrcoef(L, R);
%                 zeroLagVal = cc(1,2);
% 
%                 % full cross-correlation lag profile (fixed length, NaN if event too short)
%                 if n - 1 >= fullMaxLagSamp
%                     [xc, ~] = xcorr(L - mean(L), R - mean(R), fullMaxLagSamp, 'coeff');
%                     xcVec = xc(:);
%                 else
%                     xcVec = nan(2*fullMaxLagSamp+1, 1);
%                 end
% 
%                 % coherence (fixed frequency axis regardless of event length)
%                 try
%                     [Cxy, ~] = mscohere(L, R, [], [], nfftCoh, fps);
%                     CxyVec = Cxy(:);
%                     if numel(CxyVec) ~= numel(fCohFixed)
%                         CxyVec = nan(size(fCohFixed));
%                     end
%                 catch
%                     CxyVec = nan(size(fCohFixed));
%                 end
% 
%                 % DTW distance, normalized by average sequence length
%                 if haveDTW
%                     try
%                         dtwDist = dtw(L, R) / n;
%                     catch
%                         dtwDist = NaN;
%                     end
%                 else
%                     dtwDist = NaN;
%                 end
% 
%                 % asymmetry-index trace, full resolution
%                 Lz = zscoreSafe(L);
%                 Rz = zscoreSafe(R);
%                 AItrace = abs(Lz - Rz);
%                 xNormAI = linspace(0, 100, n)';
% 
%                 pairData(p).zeroLagCorr(end+1,1)         = zeroLagVal;
%                 pairData(p).xcorrProfiles(:,end+1)       = xcVec;
%                 pairData(p).coherenceProfiles(:,end+1)   = CxyVec;
%                 pairData(p).dtwNormDist(end+1,1)         = dtwDist;
%                 pairData(p).AItraces{end+1}              = AItrace;
%                 pairData(p).AIxNorm{end+1}               = xNormAI;
%                 pairData(p).peakLeft(end+1,1)            = max(L);
%                 pairData(p).peakRight(end+1,1)           = max(R);
% 
%                 eiPair.recording   = ei.recording;
%                 eiPair.rowIndex    = ei.rowIndex;
%                 eiPair.segIndex    = ei.segIndex;
%                 eiPair.startTime   = ei.startTime;
%                 eiPair.endTime     = ei.endTime;
%                 eiPair.durationSec = ei.durationSec;
%                 pairData(p).eventInfo(end+1) = eiPair;
%             end
%         end
%     end
% end
% 
% %% ---- 5) Package results ----
% results.keypoints = struct('name', {}, 'rawTraces', {}, 'xNorm', {}, 'eventInfo', {});
% for k = 1:numel(kpNames)
%     results.keypoints(k).name       = kpNames{k};
%     results.keypoints(k).rawTraces  = kpData(k).rawTraces;
%     results.keypoints(k).xNorm      = kpData(k).xNorm;
%     results.keypoints(k).eventInfo  = kpData(k).eventInfo;
% end
% 
% results.pairwise = struct('pairName', {}, 'leftName', {}, 'rightName', {}, ...
%     'zeroLagCorr', {}, 'lagsSec', {}, 'xcorrProfiles', {}, ...
%     'coherenceFreq', {}, 'coherenceProfiles', {}, 'dtwNormDist', {}, ...
%     'AItraces', {}, 'AIxNorm', {}, 'peakLeft', {}, 'peakRight', {}, 'eventInfo', {});
% for p = 1:numel(pairsList)
%     results.pairwise(p).pairName          = pairsList(p).pairName;
%     results.pairwise(p).leftName          = pairsList(p).leftName;
%     results.pairwise(p).rightName         = pairsList(p).rightName;
%     results.pairwise(p).zeroLagCorr       = pairData(p).zeroLagCorr;
%     results.pairwise(p).lagsSec           = lagsSecFixed;
%     results.pairwise(p).xcorrProfiles     = pairData(p).xcorrProfiles;
%     results.pairwise(p).coherenceFreq     = fCohFixed;
%     results.pairwise(p).coherenceProfiles = pairData(p).coherenceProfiles;
%     results.pairwise(p).dtwNormDist       = pairData(p).dtwNormDist;
%     results.pairwise(p).AItraces          = pairData(p).AItraces;
%     results.pairwise(p).AIxNorm           = pairData(p).AIxNorm;
%     results.pairwise(p).peakLeft          = pairData(p).peakLeft;
%     results.pairwise(p).peakRight         = pairData(p).peakRight;
%     results.pairwise(p).eventInfo         = pairData(p).eventInfo;
% end
% 
% if ~isempty(globalCorrEventMatrices)
%     meanCorrMat = mean(cat(3, globalCorrEventMatrices{:}), 3, 'omitnan');
% else
%     meanCorrMat = [];
% end
% results.globalCorr.meanMatrix       = meanCorrMat;
% results.globalCorr.keypointNames    = kpNames;
% results.globalCorr.perEventMatrices = globalCorrEventMatrices;
% 
% if ~isempty(options.OutputPath)
%     save(options.OutputPath, 'results');
% end
% 
% %% ---- 6) Plot ----
% if options.MakePlots
%     plotKeypointTraces(results.keypoints, options.LineAlpha);
%     if options.PairwiseAnalysis && ~isempty(results.pairwise)
%         plotAItraces(results.pairwise, options.LineAlpha);
%         plotPairwiseSummary(results.pairwise);
%     end
%     if options.GlobalCorrAnalysis && ~isempty(results.globalCorr.meanMatrix)
%         plotGlobalCorrDendrogram(results.globalCorr, haveLinkage);
%     end
% end
% 
% end % analyzeREMKeypointDisplacement
% 
% 
% %% =====================================================================
% function kp = loadDLCcsv(csvFile)
% % Loads a DeepLabCut-style csv with 3 header rows:
% %   row1: scorer / row2: bodyparts (x,y,likelihood x3 per keypoint) / row3: coords
% % Returns kp: struct array, one per keypoint, fields .name, .frameIdx, .x, .y, .likelihood
% 
% fid = fopen(csvFile, 'r');
% if fid == -1
%     error('Could not open file: %s', csvFile);
% end
% line1 = fgetl(fid); %#ok<NASGU>
% line2 = fgetl(fid);
% line3 = fgetl(fid);
% fclose(fid);
% 
% if numel(strfind(line2, ',')) >= numel(strfind(line2, sprintf('\t')))
%     delim = ',';
% else
%     delim = sprintf('\t');
% end
% 
% bodypartsCells = strsplit(line2, delim);
% bodypartsCells(1) = [];
% 
% if mod(numel(bodypartsCells), 3) ~= 0
%     error('Unexpected header format in %s: bodyparts columns not a multiple of 3.', csvFile);
% end
% 
% keypointNames = bodypartsCells(1:3:end);
% 
% data = readmatrix(csvFile, 'NumHeaderLines', 3, 'Delimiter', delim);
% frameIdx = data(:,1);
% coordData = data(:,2:end);
% 
% nKp = numel(keypointNames);
% kp = struct('name', cell(1,nKp), 'frameIdx', cell(1,nKp), 'x', cell(1,nKp), ...
%     'y', cell(1,nKp), 'likelihood', cell(1,nKp));
% for k = 1:nKp
%     colBase = (k-1)*3;
%     kp(k).name       = keypointNames{k};
%     kp(k).frameIdx   = frameIdx;
%     kp(k).x          = coordData(:, colBase+1);
%     kp(k).y          = coordData(:, colBase+2);
%     kp(k).likelihood = coordData(:, colBase+3);
% end
% 
% end % loadDLCcsv
% 
% 
% %% =====================================================================
% function pairs = findSymmetricPairs(kpNames)
% % Detects "left_X" / "right_X" keypoint pairs by name.
% pairs = struct('pairName', {}, 'leftName', {}, 'leftIdx', {}, 'rightName', {}, 'rightIdx', {});
% for i = 1:numel(kpNames)
%     nm = kpNames{i};
%     if startsWith(nm, 'left_')
%         suffix = nm(6:end);
%         rightNm = ['right_' suffix];
%         j = find(strcmp(kpNames, rightNm), 1);
%         if ~isempty(j)
%             pairs(end+1).pairName  = suffix; %#ok<AGROW>
%             pairs(end).leftName    = nm;
%             pairs(end).leftIdx     = i;
%             pairs(end).rightName   = rightNm;
%             pairs(end).rightIdx    = j;
%         end
%     end
% end
% end % findSymmetricPairs
% 
% 
% %% =====================================================================
% function z = zscoreSafe(x)
% s = std(x);
% if s > 0
%     z = (x - mean(x)) / s;
% else
%     z = x - mean(x);
% end
% end % zscoreSafe
% 
% 
% %% =====================================================================
% function applyDarkAxes(ax)
% bgColor   = [0.12 0.12 0.12];
% textColor = [0.92 0.92 0.92];
% gridColor = [0.4 0.4 0.4];
% set(ancestor(ax,'figure'), 'Color', bgColor);
% set(ax, 'Color', bgColor, 'XColor', textColor, 'YColor', textColor, ...
%     'GridColor', gridColor, 'Box', 'on');
% end % applyDarkAxes
% 
% 
% %% =====================================================================
% function plotKeypointTraces(kpResults, lineAlpha)
% % One figure per keypoint. All events overlaid, semi-transparent, each
% % event plotted at its own native sample count on its own 0-100% x-axis.
% lineColor = [0.30 0.75 0.95];
% textColor = [0.92 0.92 0.92];
% 
% for k = 1:numel(kpResults)
%     traces = kpResults(k).rawTraces;
%     xNormAll = kpResults(k).xNorm;
%     if isempty(traces)
%         continue
%     end
% 
%     fig = figure('Name', kpResults(k).name);
%     ax = axes(fig);
%     applyDarkAxes(ax);
%     hold(ax, 'on');
% 
%     for e = 1:numel(traces)
%         h = plot(ax, xNormAll{e}, traces{e}, 'Color', lineColor);
%         h.Color(4) = lineAlpha;
%     end
% 
%     xlabel(ax, '% of REM event completed', 'Color', textColor);
%     ylabel(ax, 'Displacement from baseline (pixels)', 'Color', textColor);
%     title(ax, sprintf('%s (n = %d events)', kpResults(k).name, numel(traces)), ...
%         'Interpreter', 'none', 'Color', textColor);
%     xlim(ax, [0 100]);
%     grid(ax, 'on');
%     hold(ax, 'off');
% end
% end % plotKeypointTraces
% 
% 
% %% =====================================================================
% function plotAItraces(pairResults, lineAlpha)
% % One figure per L/R pair: asymmetry-index |zscore(L)-zscore(R)| traces,
% % overlaid semi-transparent, full resolution, each on its own 0-100% axis.
% lineColor = [0.95 0.55 0.25];
% textColor = [0.92 0.92 0.92];
% 
% for p = 1:numel(pairResults)
%     traces = pairResults(p).AItraces;
%     xNormAll = pairResults(p).AIxNorm;
%     if isempty(traces)
%         continue
%     end
% 
%     fig = figure('Name', ['AI_' pairResults(p).pairName]);
%     ax = axes(fig);
%     applyDarkAxes(ax);
%     hold(ax, 'on');
% 
%     for e = 1:numel(traces)
%         h = plot(ax, xNormAll{e}, traces{e}, 'Color', lineColor);
%         h.Color(4) = lineAlpha;
%     end
% 
%     xlabel(ax, '% of REM event completed', 'Color', textColor);
%     ylabel(ax, '|z(left) - z(right)|  (asymmetry index)', 'Color', textColor);
%     title(ax, sprintf('%s: L/R asymmetry (n = %d events)', pairResults(p).pairName, numel(traces)), ...
%         'Interpreter', 'none', 'Color', textColor);
%     xlim(ax, [0 100]);
%     grid(ax, 'on');
%     hold(ax, 'off');
% end
% end % plotAItraces
% 
% 
% %% =====================================================================
% function plotPairwiseSummary(pairResults)
% % Summary figures comparing L/R symmetry across pairs:
% %   1) zero-lag correlation + DTW distance (jittered scatter per pair)
% %   2) mean cross-correlation lag profile ± SEM per pair
% %   3) mean coherence spectrum ± SEM per pair
% %   4) paired peak displacement (L vs R) scatter with unity line, per pair
% textColor = [0.92 0.92 0.92];
% dotColor  = [0.30 0.75 0.95];
% meanColor = [0.95 0.55 0.25];
% nPairs = numel(pairResults);
% if nPairs == 0
%     return
% end
% pairNames = {pairResults.pairName};
% 
% %% Figure 1: zero-lag correlation + DTW distance
% fig = figure('Name', 'LR_zeroLagCorr_and_DTW');
% t = tiledlayout(fig, 1, 2);
% 
% ax1 = nexttile(t); applyDarkAxes(ax1); hold(ax1, 'on');
% for p = 1:nPairs
%     v = pairResults(p).zeroLagCorr;
%     v = v(isfinite(v));
%     jitterScatter(ax1, p, v, dotColor);
%     if ~isempty(v)
%         errorbar(ax1, p, mean(v), std(v)/sqrt(numel(v)), 'o', 'Color', meanColor, ...
%             'MarkerFaceColor', meanColor, 'LineWidth', 1.5, 'CapSize', 6);
%     end
% end
% xlim(ax1, [0.5 nPairs+0.5]); xticks(ax1, 1:nPairs); xticklabels(ax1, strrep(pairNames,'_',' '));
% ylabel(ax1, 'Zero-lag correlation (r)', 'Color', textColor);
% title(ax1, 'L vs R zero-lag correlation', 'Color', textColor);
% yline(ax1, 0, 'Color', [0.6 0.6 0.6]);
% grid(ax1, 'on'); hold(ax1, 'off');
% 
% ax2 = nexttile(t); applyDarkAxes(ax2); hold(ax2, 'on');
% for p = 1:nPairs
%     v = pairResults(p).dtwNormDist;
%     v = v(isfinite(v));
%     jitterScatter(ax2, p, v, dotColor);
%     if ~isempty(v)
%         errorbar(ax2, p, mean(v), std(v)/sqrt(numel(v)), 'o', 'Color', meanColor, ...
%             'MarkerFaceColor', meanColor, 'LineWidth', 1.5, 'CapSize', 6);
%     end
% end
% xlim(ax2, [0.5 nPairs+0.5]); xticks(ax2, 1:nPairs); xticklabels(ax2, strrep(pairNames,'_',' '));
% ylabel(ax2, 'DTW distance / mean(n samples)', 'Color', textColor);
% title(ax2, 'L vs R shape dissimilarity (DTW)', 'Color', textColor);
% grid(ax2, 'on'); hold(ax2, 'off');
% 
% %% Figure 2: mean cross-correlation lag profile per pair
% fig = figure('Name', 'LR_xcorr_profiles');
% t = tiledlayout(fig, 1, nPairs);
% for p = 1:nPairs
%     ax = nexttile(t); applyDarkAxes(ax); hold(ax, 'on');
%     M = pairResults(p).xcorrProfiles;
%     lagsSec = pairResults(p).lagsSec;
%     if ~isempty(M)
%         mu = mean(M, 2, 'omitnan');
%         n  = sum(isfinite(M), 2);
%         sem = std(M, 0, 2, 'omitnan') ./ max(sqrt(n), 1);
%         valid = isfinite(mu);
%         fill(ax, [lagsSec(valid); flipud(lagsSec(valid))], ...
%             [mu(valid)-sem(valid); flipud(mu(valid)+sem(valid))], ...
%             dotColor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%         plot(ax, lagsSec, mu, 'Color', dotColor, 'LineWidth', 2);
%     end
%     xline(ax, 0, 'Color', [0.6 0.6 0.6]);
%     yline(ax, 0, 'Color', [0.4 0.4 0.4]);
%     xlabel(ax, 'Lag (s)', 'Color', textColor);
%     ylabel(ax, 'Cross-corr (coeff)', 'Color', textColor);
%     title(ax, pairResults(p).pairName, 'Interpreter', 'none', 'Color', textColor);
%     grid(ax, 'on'); hold(ax, 'off');
% end
% sgtitle(t, 'L vs R cross-correlation (mean \pm SEM)', 'Color', textColor);
% 
% %% Figure 3: mean coherence spectrum per pair
% fig = figure('Name', 'LR_coherence');
% t = tiledlayout(fig, 1, nPairs);
% for p = 1:nPairs
%     ax = nexttile(t); applyDarkAxes(ax); hold(ax, 'on');
%     M = pairResults(p).coherenceProfiles;
%     f = pairResults(p).coherenceFreq;
%     if ~isempty(M)
%         mu = mean(M, 2, 'omitnan');
%         n  = sum(isfinite(M), 2);
%         sem = std(M, 0, 2, 'omitnan') ./ max(sqrt(n), 1);
%         valid = isfinite(mu);
%         fill(ax, [f(valid); flipud(f(valid))], ...
%             [mu(valid)-sem(valid); flipud(mu(valid)+sem(valid))], ...
%             dotColor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%         plot(ax, f, mu, 'Color', dotColor, 'LineWidth', 2);
%     end
%     ylim(ax, [0 1]);
%     xlabel(ax, 'Frequency (Hz)', 'Color', textColor);
%     ylabel(ax, 'Coherence', 'Color', textColor);
%     title(ax, pairResults(p).pairName, 'Interpreter', 'none', 'Color', textColor);
%     grid(ax, 'on'); hold(ax, 'off');
% end
% sgtitle(t, 'L vs R coherence (mean \pm SEM)', 'Color', textColor);
% 
% %% Figure 4: paired peak displacement scatter (L vs R), unity line
% fig = figure('Name', 'LR_peak_displacement_paired');
% t = tiledlayout(fig, 1, nPairs);
% for p = 1:nPairs
%     ax = nexttile(t); applyDarkAxes(ax); hold(ax, 'on');
%     pl = pairResults(p).peakLeft;
%     pr = pairResults(p).peakRight;
%     valid = isfinite(pl) & isfinite(pr);
%     scatter(ax, pl(valid), pr(valid), 30, dotColor, 'filled', 'MarkerFaceAlpha', 0.5);
%     if any(valid)
%         m = max([pl(valid); pr(valid)]);
%         plot(ax, [0 m], [0 m], '--', 'Color', [0.6 0.6 0.6]);
%         axis(ax, [0 m 0 m]);
%         axis(ax, 'square');
%     end
%     xlabel(ax, 'Peak displacement, left (px)', 'Color', textColor);
%     ylabel(ax, 'Peak displacement, right (px)', 'Color', textColor);
%     title(ax, pairResults(p).pairName, 'Interpreter', 'none', 'Color', textColor);
%     grid(ax, 'on'); hold(ax, 'off');
% end
% sgtitle(t, 'Peak displacement per event: left vs right', 'Color', textColor);
% 
% end % plotPairwiseSummary
% 
% 
% %% =====================================================================
% function jitterScatter(ax, x0, data, color)
% if isempty(data)
%     return
% end
% jitter = (rand(size(data)) - 0.5) * 0.3;
% scatter(ax, x0 + jitter, data, 25, color, 'filled', 'MarkerFaceAlpha', 0.5);
% end % jitterScatter
% 
% 
% %% =====================================================================
% function plotGlobalCorrDendrogram(globalCorr, haveLinkage)
% % Cross-keypoint mean correlation matrix (heatmap) + optional dendrogram,
% % for comparing structure against the ROI-based PCA/dendrogram analysis.
% textColor = [0.92 0.92 0.92];
% M = globalCorr.meanMatrix;
% names = strrep(globalCorr.keypointNames, '_', ' ');
% n = size(M, 1);
% 
% if haveLinkage
%     fig = figure('Name', 'Global_keypoint_corr_and_dendrogram');
%     t = tiledlayout(fig, 1, 2);
% else
%     fig = figure('Name', 'Global_keypoint_corr');
%     t = tiledlayout(fig, 1, 1);
% end
% set(fig, 'Color', [0.12 0.12 0.12]);
% 
% ax1 = nexttile(t); applyDarkAxes(ax1);
% imagesc(ax1, M); axis(ax1, 'square'); clim(ax1, [-1 1]); colorbar(ax1);
% colormap(ax1, 'turbo');
% xticks(ax1, 1:n); yticks(ax1, 1:n);
% xticklabels(ax1, names); yticklabels(ax1, names);
% xtickangle(ax1, 45);
% title(ax1, 'Mean cross-keypoint correlation (REM events)', 'Color', textColor);
% 
% if haveLinkage
%     ax2 = nexttile(t); applyDarkAxes(ax2);
%     D = 1 - M;
%     D = (D + D') / 2;
%     D(1:n+1:end) = 0;
%     Z = linkage(squareform(D, 'tovector'), 'average');
%     [~, ~, outperm] = dendrogram(ax2, Z, 0, 'Labels', names);
%     set(ax2, 'XTickLabelRotation', 45);
%     title(ax2, 'Hierarchical clustering (1 - correlation)', 'Color', textColor); %#ok<NASGU>
% end
% 
% sgtitle(t, 'Cross-keypoint structure during REM', 'Color', textColor);
% end % plotGlobalCorrDendrogram

% function results = analyzeREMKeypointDisplacement(excelPath, csvFolder, fps, options)
% % ANALYZEREMKEYPOINTDISPLACEMENT
% % For each REM event listed in an excel file, finds the matching DLC
% % keypoint CSV, computes a 10-second-prior baseline (x,y) position for
% % each keypoint, then computes the frame-by-frame displacement magnitude
% % from that baseline across the REM event. Traces are time-normalized to
% % 0-100% of event duration and plotted overlaid (semi-transparent) per
% % keypoint, one figure per keypoint.
% %
% % USAGE:
% %   results = analyzeREMKeypointDisplacement(excelPath, csvFolder, fps)
% %   results = analyzeREMKeypointDisplacement(excelPath, csvFolder, fps, ...
% %                 'Sheet', 2, 'BaselineDuration', 10, ...
% %                 'OutputPath', 'results.mat', 'MakePlots', true)
% %
% % INPUTS:
% %   excelPath  - path to REM excel file (same format as prior wake/REM code:
% %                Sheet 2, col 1 = tdms/recording path, cols 2:end = pairs of
% %                [start end] REM segment times in seconds, NaN-padded)
% %   csvFolder  - folder containing DLC keypoint CSVs. CSV filenames are
% %                expected to START WITH the same base filename as column 1
% %                of the excel file (extension/path stripped).
% %   fps        - frame rate (frames/sec) of the video/keypoint tracking.
% %                REQUIRED because it is not encoded in either file.
% %
% % NAME-VALUE OPTIONS:
% %   'Sheet'             - excel sheet to read (default 2)
% %   'BaselineDuration'  - seconds prior to event onset used as baseline (default 10)
% %   'OutputPath'        - if provided, saves `results` struct to this .mat path
% %   'MakePlots'         - true/false, whether to generate the overlay plots (default true)
% %   'LineAlpha'         - transparency of each individual event trace (default 0.15)
% %
% % OUTPUT:
% %   results - struct array, one element per keypoint, with fields:
% %       .name          - keypoint name (from DLC header)
% %       .rawTraces     - {1 x nEvents} cell array of displacement magnitude
% %                        traces, one per event, at FULL original sample
% %                        resolution (no downsampling -- each event keeps
% %                        every frame it had).
% %       .xNorm         - {1 x nEvents} cell array, same size as rawTraces,
% %                        giving each trace's own x-axis rescaled to 0-100%
% %                        (so traces with different frame counts still plot
% %                        on the same 0-100% axis without resampling values).
% %       .eventInfo     - struct array, one per event, with fields
% %                        .recording (csv file name), .rowIndex (row in
% %                        excel), .segIndex (segment # within that row),
% %                        .startTime, .endTime, .durationSec, .nFramesBaseline
% %
% % NOTES / ASSUMPTIONS:
% %   - Segment times in the excel file are assumed to be in SECONDS.
% %   - Baseline = mean x,y position of the keypoint over the
% %     BaselineDuration seconds immediately preceding event onset.
% %   - Displacement is the Euclidean (magnitude) distance from that
% %     baseline (x0,y0): sqrt((x-x0)^2 + (y-y0)^2). No left/right sign
% %     information is retained here -- for the eventual L-vs-R symmetry
% %     analysis you'll likely want signed x (and/or y) displacement too,
% %     which can be pulled from the same rawTraces machinery later.
% %   - Traces are NOT resampled/downsampled -- every frame in every event
% %     is kept, so brief twitches aren't smoothed or lost. Each event's
% %     x-axis is independently rescaled to span 0-100%, so events of
% %     different durations/frame counts still overlay on a common axis.
% %   - Events with insufficient baseline data (too close to start of
% %     recording) or a missing/ambiguous CSV match are skipped with a
% %     warning, not silently dropped.
% 
% arguments
%     excelPath (1,:) char
%     csvFolder (1,:) char
%     fps (1,1) double {mustBePositive}
%     options.Sheet (1,1) double = 2
%     options.BaselineDuration (1,1) double = 10
%     options.OutputPath (1,:) char = ''
%     options.MakePlots (1,1) logical = true
%     options.LineAlpha (1,1) double = 0.15
% end
% 
% %% ---- 1) Read REM excel file (same logic as prior wake/REM script) ----
% T = readtable(excelPath, 'Sheet', options.Sheet, 'ReadVariableNames', false);
% if strcmp(T{1,1}{1}, 'file name')
%     T = T(2:end,:);
% end
% nRec = height(T);
% 
% %% ---- 2) Containers keyed by keypoint name ----
% kpNames = {};          % will be populated from first successfully-read CSV
% kpData  = struct('rawTraces', {}, 'xNorm', {}, 'eventInfo', {});
% 
% %% ---- 3) Loop over recordings / REM segments ----
% for r = 1:nRec
%     tdmsPath = T{r,1}{1};
%     if isempty(tdmsPath) || (isstring(tdmsPath) && strlength(tdmsPath)==0)
%         continue
%     end
% 
%     segRaw = table2array(T(r, 2:end));
%     segRaw = segRaw(~isnan(segRaw));
%     if isempty(segRaw)
%         continue
%     end
%     segments = reshape(segRaw, 2, [])';   % nSeg x 2 : [start end] in seconds
% 
%     % ---- Find matching CSV ----
%     [~, baseName, ~] = fileparts(tdmsPath);
%     csvMatches = dir(fullfile(csvFolder, [baseName '*.csv']));
%     if isempty(csvMatches)
%         warning('No CSV found for recording "%s" (row %d). Skipping.', baseName, r);
%         continue
%     elseif numel(csvMatches) > 1
%         warning('Multiple CSVs match "%s" (row %d): using "%s". Rename files if this is wrong.', ...
%             baseName, r, csvMatches(1).name);
%     end
%     csvFile = fullfile(csvMatches(1).folder, csvMatches(1).name);
% 
%     % ---- Load DLC csv (cached per row so we don't re-read for every segment) ----
%     try
%         kp = loadDLCcsv(csvFile);
%     catch ME
%         warning('Failed to load CSV "%s" (row %d): %s. Skipping.', csvFile, r, ME.message);
%         continue
%     end
% 
%     if isempty(kpNames)
%         kpNames = {kp.name};
%         for k = 1:numel(kpNames)
%             kpData(k).rawTraces  = {};
%             kpData(k).xNorm      = {};
%             kpData(k).eventInfo  = struct('recording',{}, 'rowIndex',{}, 'segIndex',{}, ...
%                 'startTime',{}, 'endTime',{}, 'durationSec',{}, 'nFramesBaseline',{});
%         end
%     elseif ~isequal(kpNames, {kp.name})
%         warning('Keypoint set in "%s" differs from earlier CSVs. Skipping recording.', csvFile);
%         continue
%     end
% 
%     nFramesTotal = numel(kp(1).frameIdx);
% 
%     % ---- Loop over REM segments within this recording ----
%     for s = 1:size(segments,1)
%         startT = segments(s,1);
%         endT   = segments(s,2);
%         if isnan(startT) || isnan(endT) || endT <= startT
%             continue
%         end
% 
%         startFrame = round(startT * fps);          % 0-based frame number
%         endFrame   = round(endT   * fps);
%         baseStartFrame = startFrame - round(options.BaselineDuration * fps);
% 
%         % Convert 0-based DLC frame numbers -> 1-based row indices,
%         % assuming frameIdx is contiguous starting at 0.
%         baseStartRow = baseStartFrame + 1;
%         eventStartRow = startFrame + 1;
%         eventEndRow   = endFrame + 1;
% 
%         if baseStartRow < 1
%             warning('Row %d, segment %d: insufficient baseline (event too close to recording start). Skipping.', r, s);
%             continue
%         end
%         if eventEndRow > nFramesTotal
%             warning('Row %d, segment %d: event extends past end of CSV data. Truncating to available frames.', r, s);
%             eventEndRow = nFramesTotal;
%         end
%         if eventStartRow > eventEndRow
%             warning('Row %d, segment %d: no valid frames in event window. Skipping.', r, s);
%             continue
%         end
% 
%         baselineRows = baseStartRow:(eventStartRow-1);
%         eventRows    = eventStartRow:eventEndRow;
% 
%         % ---- Per-keypoint baseline + displacement ----
%         for k = 1:numel(kp)
%             x0 = mean(kp(k).x(baselineRows), 'omitnan');
%             y0 = mean(kp(k).y(baselineRows), 'omitnan');
% 
%             dx = kp(k).x(eventRows) - x0;
%             dy = kp(k).y(eventRows) - y0;
%             dispTrace = sqrt(dx.^2 + dy.^2);   % raw displacement magnitude trace, full resolution
% 
%             nRaw = numel(dispTrace);
%             if nRaw < 2
%                 continue
%             end
% 
%             % No resampling: keep every sample. Just rescale this event's
%             % own x-axis to span 0-100%, so events with different frame
%             % counts still overlay on a common axis without altering values.
%             xNormThisEvent = linspace(0, 100, nRaw);
% 
%             kpData(k).rawTraces{end+1} = dispTrace(:);
%             kpData(k).xNorm{end+1}     = xNormThisEvent(:);
% 
%             ei.recording       = csvMatches(1).name;
%             ei.rowIndex        = r;
%             ei.segIndex        = s;
%             ei.startTime       = startT;
%             ei.endTime         = endT;
%             ei.durationSec     = endT - startT;
%             ei.nFramesBaseline = numel(baselineRows);
%             kpData(k).eventInfo(end+1) = ei;
%         end
%     end
% end
% 
% %% ---- 4) Package results ----
% results = struct('name', {}, 'rawTraces', {}, 'xNorm', {}, 'eventInfo', {});
% for k = 1:numel(kpNames)
%     results(k).name       = kpNames{k};
%     results(k).rawTraces  = kpData(k).rawTraces;
%     results(k).xNorm      = kpData(k).xNorm;
%     results(k).eventInfo  = kpData(k).eventInfo;
% end
% 
% if ~isempty(options.OutputPath)
%     save(options.OutputPath, 'results');
% end
% 
% %% ---- 5) Plot ----
% if options.MakePlots
%     plotKeypointTraces(results, options.LineAlpha);
% end
% 
% end % analyzeREMKeypointDisplacement
% 
% 
% %% =====================================================================
% function kp = loadDLCcsv(csvFile)
% % Loads a DeepLabCut-style csv with 3 header rows:
% %   row1: scorer
% %   row2: bodyparts (repeated 3x per keypoint: x,y,likelihood columns)
% %   row3: coords (x/y/likelihood labels)
% % Returns kp: struct array, one element per keypoint, with fields
% %   .name, .frameIdx, .x, .y, .likelihood
% 
% fid = fopen(csvFile, 'r');
% if fid == -1
%     error('Could not open file: %s', csvFile);
% end
% line1 = fgetl(fid); %#ok<NASGU> % scorer row, unused
% line2 = fgetl(fid); % bodyparts row
% line3 = fgetl(fid); % coords row
% fclose(fid);
% 
% % detect delimiter (DLC default is comma; be tolerant of tabs too)
% if numel(strfind(line2, ',')) >= numel(strfind(line2, sprintf('\t')))
%     delim = ',';
% else
%     delim = sprintf('\t');
% end
% 
% bodypartsCells = strsplit(line2, delim);
% bodypartsCells(1) = [];   % drop leading "bodyparts" label / empty cell
% 
% if mod(numel(bodypartsCells), 3) ~= 0
%     error('Unexpected header format in %s: bodyparts columns not a multiple of 3.', csvFile);
% end
% 
% keypointNames = bodypartsCells(1:3:end);
% 
% data = readmatrix(csvFile, 'NumHeaderLines', 3, 'Delimiter', delim);
% frameIdx = data(:,1);
% coordData = data(:,2:end);
% 
% nKp = numel(keypointNames);
% kp = struct('name', cell(1,nKp), 'frameIdx', cell(1,nKp), 'x', cell(1,nKp), ...
%     'y', cell(1,nKp), 'likelihood', cell(1,nKp));
% for k = 1:nKp
%     colBase = (k-1)*3;
%     kp(k).name       = keypointNames{k};
%     kp(k).frameIdx   = frameIdx;
%     kp(k).x          = coordData(:, colBase+1);
%     kp(k).y          = coordData(:, colBase+2);
%     kp(k).likelihood = coordData(:, colBase+3);
% end
% 
% end % loadDLCcsv
% 
% 
% %% =====================================================================
% function plotKeypointTraces(results, lineAlpha)
% % One figure per keypoint, styled for a dark background. All events for
% % that keypoint overlaid, each as a semi-transparent line. No resampling:
% % every event is plotted at its own native sample count, with its own
% % x-axis independently rescaled to span 0-100%.
% 
% bgColor   = [0.12 0.12 0.12];
% axColor   = [0.12 0.12 0.12];
% textColor = [0.92 0.92 0.92];
% gridColor = [0.4 0.4 0.4];
% lineColor = [0.30 0.75 0.95];   % bright cyan-blue, visible on dark bg
% 
% for k = 1:numel(results)
%     traces = results(k).rawTraces;
%     xNormAll = results(k).xNorm;
%     if isempty(traces)
%         continue
%     end
% 
%     fig = figure('Name', results(k).name, 'Color', bgColor);
%     ax = axes(fig);
%     set(ax, 'Color', axColor, 'XColor', textColor, 'YColor', textColor, ...
%         'GridColor', gridColor, 'Box', 'on');
%     hold(ax, 'on');
% 
%     for e = 1:numel(traces)
%         h = plot(ax, xNormAll{e}, traces{e}, 'Color', lineColor);
%         h.Color(4) = lineAlpha;   % set alpha (transparency) on the line
%     end
% 
%     xlabel(ax, '% of REM event completed', 'Color', textColor);
%     ylabel(ax, 'Displacement from baseline (pixels)', 'Color', textColor);
%     title(ax, sprintf('%s (n = %d events)', results(k).name, numel(traces)), ...
%         'Interpreter', 'none', 'Color', textColor);
%     xlim(ax, [0 100]);
%     grid(ax, 'on');
%     hold(ax, 'off');
% end
% 
% end % plotKeypointTraces