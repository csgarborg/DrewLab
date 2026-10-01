function eventTableAll = aggregateREMEventsTiming(folderPath, timingType, binSizeSeconds, durationBinSizeSeconds, rasterResolutionSeconds)
%AGGREGATEREMEVENTS  Combine per-mouse REM event CSVs into pooled figures.
%
%   eventTableAll = aggregateREMEvents(folderPath, timingType, ...
%       binSizeSeconds, durationBinSizeSeconds, rasterResolutionSeconds)
%
%   Reads every pair of CSVs in folderPath that were produced by
%   analyzeREMEvents (one pair per mouse/run):
%     '<name>_eventTIming_exact.csv' or '<name>_eventTIming_simple.csv'
%         (REM event timing; which one is loaded is controlled by
%         timingType)
%     '<name>_videoChain.csv'
%         (that mouse's video chain/offsets -- required for the raster
%         plot; a mouse missing this file is still included in the other
%         six plots, just skipped in the raster with a warning)
%
%   Stacks them into pooled data and produces the same seven figures as
%   analyzeREMEvents, but pooling every event/day from every mouse:
%     1. Histogram of REM event start times (s since session start).
%     2. Scatter: duration vs. time since previous REM event *start*.
%     3. Scatter: duration vs. time since previous REM event *end*.
%     4. Histogram of REM event durations.
%     5. Scatter: duration vs. time until next REM event *start*.
%     6. Scatter: duration vs. time until next REM event, measured from
%        this event's *end*.
%     7. Raster plot: one row per (mouse, day) *that has at least one REM
%        event* (mouse-days with none are omitted), black/blue/yellow as
%        in analyzeREMEvents, stacked across every mouse so you can
%        compare REM timing across days and animals at a glance.
%   Figures 1-6 are plotted in remColor = [0.85 0.33 0.10]; figure 7 uses
%   its own fixed black/blue/yellow scheme. Figures 2, 3, 5, and 6 (the
%   four duration-vs-gap scatter plots) each get an overlaid least-squares
%   linear fit with its equation, R^2, and p-value (for the slope)
%   annotated on the plot, computed on the pooled data across all mice.
%
%   INPUTS
%     folderPath               - folder containing the per-mouse CSVs
%     timingType                - 'exact' or 'simple' (default 'exact');
%                                  selects which event-timing CSV variant
%                                  to load (a folder may contain both per
%                                  mouse, and they shouldn't be pooled
%                                  together)
%     binSizeSeconds            - start-time histogram bin width in
%                                  seconds (default 300)
%     durationBinSizeSeconds    - duration histogram bin width in seconds
%                                  (default 30)
%     rasterResolutionSeconds   - time-bin width (s) for the raster plot
%                                  (default 5)
%
%   OUTPUT
%     eventTableAll - all per-mouse event tables vertically concatenated,
%                     with an added 'mouse' column identifying the source
%                     file (filename with the suffix stripped)
%
%   Example:
%     allEvents = aggregateREMEvents('C:\data\rem_csvs', 'exact');

if nargin < 5 || isempty(rasterResolutionSeconds), rasterResolutionSeconds = 5;  end
if nargin < 4 || isempty(durationBinSizeSeconds),  durationBinSizeSeconds = 30;  end
if nargin < 3 || isempty(binSizeSeconds),          binSizeSeconds = 300;        end
if nargin < 2 || isempty(timingType),              timingType = 'exact';       end

if ~any(strcmpi(timingType, {'exact','simple'}))
    error('timingType must be ''exact'' or ''simple''.');
end

remColor = [0.85 0.33 0.10];
offsetCol = 'offset_actual_s';
if strcmpi(timingType, 'simple')
    offsetCol = 'offset_nominal_s';
end

suffix = sprintf('_eventTIming_%s', lower(timingType));
files = dir(fullfile(folderPath, ['*' suffix '.csv']));

if isempty(files)
    error('No files matching *%s.csv found in %s', suffix, folderPath);
end

tables = cell(numel(files), 1);
rasterRows = {};
rasterLabels = {};
skippedRasterMice = {};

for i = 1:numel(files)
    fPath = fullfile(files(i).folder, files(i).name);
    t = readtable(fPath);
    mouseID = erase(files(i).name, [suffix '.csv']);

    % Keep a reliable datetime copy of 'day' for raster matching before
    % standardizing the column to string (needed so concatenation across
    % mice never fails on a type mismatch -- not needed for the plots).
    dayDatetime = t.day;
    if ~isdatetime(dayDatetime)
        try
            dayDatetime = datetime(dayDatetime);
        catch
            dayDatetime = NaT(height(t), 1);
        end
    end

    % ---- raster contribution for this mouse ----
    vcPath = fullfile(folderPath, [mouseID '_videoChain.csv']);
    if exist(vcPath, 'file')
        vc = readtable(vcPath);
        if ~isdatetime(vc.date)
            try
                vc.date = datetime(vc.date);
            catch
                vc.date = NaT(height(vc), 1);
            end
        end
        etForRaster = t;
        etForRaster.day = dayDatetime;
        [mouseRows, mouseLabels] = buildMouseDayRaster(vc, etForRaster, offsetCol, rasterResolutionSeconds, mouseID);
        rasterRows = [rasterRows; mouseRows]; %#ok<AGROW>
        rasterLabels = [rasterLabels; mouseLabels]; %#ok<AGROW>
    else
        skippedRasterMice{end+1} = mouseID; %#ok<AGROW>
    end

    % ---- standardize for pooled scatter/histogram plots ----
    t.day = string(dayDatetime);
    t.mouse = repmat(string(mouseID), height(t), 1);
    tables{i} = t;
end

if ~isempty(skippedRasterMice)
    warning('No matching _videoChain.csv found for %d mouse/mice (skipped in raster only): %s', ...
        numel(skippedRasterMice), strjoin(skippedRasterMice, ', '));
end

eventTableAll = vertcat(tables{:});

fprintf('Loaded %d file(s) / mice from %s (%s timing).\n', numel(files), folderPath, timingType);
fprintf('Pooled %d total REM events.\n', height(eventTableAll));

% ---- Plots (kept as MATLAB figures, not saved as images) ----
plotHistogramFcn(eventTableAll, binSizeSeconds, remColor);
plotDurationVsGapFcn(eventTableAll, remColor);
plotDurationVsGapEndFcn(eventTableAll, remColor);
plotDurationHistogramFcn(eventTableAll, durationBinSizeSeconds, remColor);
plotDurationVsNextStartFcn(eventTableAll, remColor);
plotDurationVsNextEndFcn(eventTableAll, remColor);

if ~isempty(rasterRows)
    plotPooledRasterFcn(rasterRows, rasterLabels, rasterResolutionSeconds);
else
    warning('No _videoChain.csv files found for any mouse -- skipping raster plot.');
end

plotDurationVsStartFcn(eventTableAll, remColor);

end


% ===================== local functions =====================
% (same plotting logic as analyzeREMEvents.m, operating on pooled data)

function fig = plotHistogramFcn(et, binSize, remColor)
maxT = max(et.abs_start_s);
edges = 0:binSize:(maxT + binSize);

fig = figure('Name', 'REM event histogram (all mice)');
subplot(1,2,1)
histogram(et.abs_start_s, edges, 'EdgeColor', 'black', 'FaceColor', remColor);
xlabel(sprintf('Time since recording session start (s) \x2014 %g s bins', binSize));
ylabel('REM event count');
title(sprintf('REM events by time since start of recording session (n = %d mice)', numel(unique(et.mouse))));

% ax = gca;
% ax2 = axes('Position', ax.Position, 'XAxisLocation', 'top', ...
%     'YAxisLocation', 'right', 'Color', 'none');
% ax2.YTick = [];
% ax2.XLim = ax.XLim / 3600;
% ax2.XLabel.String = 'Hours since start';

subplot(1,2,2)
histogram(et.abs_start_s, edges, 'Normalization', 'pdf', 'EdgeColor', 'black', 'FaceColor', remColor, 'DisplayName', 'PDF');
xlabel(sprintf('Time since recording session start (s) \x2014 %g s bins', binSize));
ylabel('REM event PDF');
title(sprintf('REM events by time since start of recording session (n = %d mice)', numel(unique(et.mouse))));

% ax = gca;
% ax2 = axes('Position', ax.Position, 'XAxisLocation', 'top', ...
%     'YAxisLocation', 'right', 'Color', 'none');
% ax2.YTick = [];
% ax2.XLim = ax.XLim / 3600;
% ax2.XLabel.String = 'Hours since start';

% Normal fit
hold on
pdREM = fitdist(et.abs_start_s,'Normal');
x = linspace(edges(1),edges(end),1000);
t = seconds(pdREM.mu);
t.Format = 'hh:mm:ss';
plot(x,pdf(pdREM,x),...
    'Color','white',...
    'LineWidth',3,...
    'DisplayName',sprintf(['PDF fit (\\mu=%.2fs, ' char(t) ')'],pdREM.mu));
legend show
end


function fig = plotDurationVsGapFcn(et, remColor)
mask = ~isnan(et.gap_from_prev_s);

fig = figure('Name', 'REM duration vs. gap from previous start (all mice)');
scatter(et.gap_from_prev_s(mask), et.duration_s(mask), 36, remColor, 'filled', ...
    'MarkerFaceAlpha', 0.6);
xlabel('Time since previous REM event start (s)');
ylabel('REM event duration (s)');
title({sprintf('REM event duration vs. time since previous REM event start (n = %d mice)', numel(unique(et.mouse))), ...
       '(first event of each day excluded)'});
addFitLine(et.gap_from_prev_s(mask), et.duration_s(mask));
end


function fig = plotDurationVsGapEndFcn(et, remColor)
mask = ~isnan(et.gap_from_prev_end_s);

fig = figure('Name', 'REM duration vs. gap from previous end (all mice)');
scatter(et.gap_from_prev_end_s(mask), et.duration_s(mask), 36, remColor, 'filled', ...
    'MarkerFaceAlpha', 0.6);
xlabel('Time since previous REM event end (s)');
ylabel('REM event duration (s)');
title({sprintf('REM event duration vs. time since previous REM event end (n = %d mice)', numel(unique(et.mouse))), ...
       '(first event of each day excluded)'});
addFitLine(et.gap_from_prev_end_s(mask), et.duration_s(mask));
end


function fig = plotDurationHistogramFcn(et, binSize, remColor)
maxD = max(et.duration_s);
edges = 0:binSize:(maxD + binSize);

fig = figure('Name', 'REM event duration histogram (all mice)');
subplot(1,2,1)
histogram(et.duration_s, edges, 'EdgeColor', 'black', 'FaceColor', remColor);
xlabel(sprintf('REM event duration (s) \x2014 %g s bins', binSize));
ylabel('REM event count');
title(sprintf('REM event durations (n = %d mice)', numel(unique(et.mouse))));

subplot(1,2,2)
histogram(et.duration_s, edges, 'Normalization', 'pdf', 'EdgeColor', 'black', 'FaceColor', remColor, 'DisplayName', 'PDF');
xlabel(sprintf('REM event duration (s) \x2014 %g s bins', binSize));
ylabel('REM event PDF');
title(sprintf('REM event durations (n = %d mice)', numel(unique(et.mouse))));

% Normal fit
hold on
pdREM = fitdist(et.duration_s,'Normal');
x = linspace(edges(1),edges(end),1000);
t = seconds(pdREM.mu);
t.Format = 'mm:ss';
plot(x,pdf(pdREM,x),...
    'Color','white',...
    'LineWidth',3,...
    'DisplayName',sprintf(['PDF fit (\\mu=%.2fs, ' char(t) ')'],pdREM.mu));
legend show
end


function fig = plotDurationVsNextStartFcn(et, remColor)
mask = ~isnan(et.gap_to_next_start_s);

fig = figure('Name', 'REM duration vs. gap until next start (all mice)');
scatter(et.duration_s(mask), et.gap_to_next_start_s(mask), 36, remColor, 'filled', ...
    'MarkerFaceAlpha', 0.6);
ylabel('Time until next REM event start (s)');
xlabel('REM event duration (s)');
title({sprintf('REM event duration vs. time until next REM event start (n = %d mice)', numel(unique(et.mouse))), ...
       '(last event of each day excluded)'});
addFitLine(et.duration_s(mask), et.gap_to_next_start_s(mask));
end


function fig = plotDurationVsNextEndFcn(et, remColor)
mask = ~isnan(et.gap_to_next_end_s);

fig = figure('Name', 'REM duration vs. gap (end) until next start (all mice)');
scatter(et.duration_s(mask), et.gap_to_next_end_s(mask), 36, remColor, 'filled', ...
    'MarkerFaceAlpha', 0.6);
ylabel('Time from this event''s end until next REM event start (s)');
xlabel('REM event duration (s)');
title({sprintf('REM event duration vs. time until next REM event (n = %d mice)', numel(unique(et.mouse))), ...
       '(measured from this event''s end; last event of each day excluded)'});
addFitLine(et.duration_s(mask), et.gap_to_next_end_s(mask));
end


function [rows, labels] = buildMouseDayRaster(vc, et, offsetCol, resolutionSeconds, mouseID)
% Build one raster row per day present in this mouse's video chain.
maxSeconds = 5 * 3600;

days = unique(vc.date, 'stable');
days = days(~isnat(days));
hasREM = ismember(days, unique(et.day));
days = days(hasREM);
nDays = numel(days);

edges = 0:resolutionSeconds:maxSeconds;
nBins = numel(edges) - 1;

rows = zeros(nDays, nBins);
labels = cell(nDays, 1);

for d = 1:nDays
    dayVideos = vc(vc.date == days(d), :);
    dayEvents = et(et.day == days(d), :);
    row = zeros(1, nBins);

    for v = 1:height(dayVideos)
        o = dayVideos.(offsetCol)(v);
        binStart = max(1, floor(o / resolutionSeconds) + 1);
        binEnd   = min(nBins, ceil((o + 3600) / resolutionSeconds));
        if binEnd >= binStart
            row(binStart:binEnd) = 1;
        end
    end

    for e = 1:height(dayEvents)
        s = dayEvents.abs_start_s(e);
        p = dayEvents.abs_stop_s(e);
        binStart = max(1, floor(s / resolutionSeconds) + 1);
        binEnd   = min(nBins, ceil(p / resolutionSeconds));
        if binEnd >= binStart
            row(binStart:binEnd) = 2;
        end
    end

    rows(d,:) = row;
    labels{d} = sprintf('%s: %s', mouseID, string(days(d), 'yyyy-MM-dd'));
end

rows = num2cell(rows, 2);
end


function fig = plotPooledRasterFcn(rasterRows, rasterLabels, resolutionSeconds)
maxSeconds = 5 * 3600;
edges = 0:resolutionSeconds:maxSeconds;
nBins = numel(edges) - 1;
binCenters = edges(1:end-1) + resolutionSeconds/2;

nRows = numel(rasterRows);
catMat = zeros(nRows, nBins);
for r = 1:nRows
    catMat(r,:) = rasterRows{r};
end

rgbImg = categoriesToRGB(catMat);

fig = figure('Name', 'REM raster by day and mouse');
image(binCenters / 3600, 1:nRows, rgbImg);
set(gca, 'YDir', 'normal');
yticks(1:nRows);
yticklabels(rasterLabels);
xlabel('Time since recording session start (hours)');
ylabel('Mouse: day');
title(sprintf('REM events across recording days (n = %d mouse-days)', nRows));
addRasterLegend();
end


function fig = plotDurationVsStartFcn(et, remColor)
mask = ~isnan(et.abs_start_s);

fig = figure('Name', 'REM duration vs. REM event start time');
scatter(et.abs_start_s(mask), et.duration_s(mask), 36, remColor, 'filled', ...
    'MarkerFaceAlpha', 0.6);
xlabel('Time from recording start to REM start (s)');
ylabel('REM event duration (s)');
title('REM event duration vs. REM event start time');
addFitLine(et.abs_start_s(mask), et.duration_s(mask));
xlim([0 3600*5])
end


function rgb = categoriesToRGB(catMat)
noDataColor = [0 0 0];
nonREMColor = [0.15 0.35 0.80];
remColorRaster = [1.00 0.84 0.00];

[nR, nC] = size(catMat);
rgb = zeros(nR, nC, 3);
for ch = 1:3
    layer = zeros(nR, nC);
    layer(catMat == 0) = noDataColor(ch);
    layer(catMat == 1) = nonREMColor(ch);
    layer(catMat == 2) = remColorRaster(ch);
    rgb(:,:,ch) = layer;
end
end


function addRasterLegend()
hold on;
h1 = patch(NaN, NaN, [0 0 0]);
h2 = patch(NaN, NaN, [0.15 0.35 0.80]);
h3 = patch(NaN, NaN, [1.00 0.84 0.00]);
legend([h1 h2 h3], {'No data','Not REM','REM'}, 'Location', 'eastoutside');
hold off;
end


function addFitLine(x, y)
% Overlay a least-squares linear fit with R^2 / p-value / n annotation.
% Uses only base MATLAB (polyfit + betainc), no toolbox required.
valid = ~isnan(x) & ~isnan(y);
x = x(valid);
y = y(valid);
n = numel(x);
if n < 3 || range(x) == 0
    return
end

p = polyfit(x, y, 1);
xFit = linspace(min(x), max(x), 100);
yFit = polyval(p, xFit);

hold on;
plot(xFit, yFit, '--w', 'LineWidth', 1.5);
hold off;

yPred = polyval(p, x);
ssRes = sum((y - yPred).^2);
ssTot = sum((y - mean(y)).^2);
if ssTot > 0
    r2 = 1 - ssRes / ssTot;
else
    r2 = NaN;
end

r = pearsonR(x, y);
df = n - 2;
if df > 0 && abs(r) < 1
    tstat = r * sqrt(df / (1 - r^2));
    pval = betainc(df / (df + tstat^2), df/2, 0.5);
else
    pval = NaN;
end

txt = sprintf('y = %.3gx + %.3g\nR^2 = %.3f\np %s\nn = %d', ...
    p(1), p(2), r2, formatPValue(pval), n);
text(0.05, 0.95, txt, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'Color', 'white', 'FontSize', 9);
end


function r = pearsonR(x, y)
r = sum((x - mean(x)) .* (y - mean(y))) / ...
    sqrt(sum((x - mean(x)).^2) * sum((y - mean(y)).^2));
end


function s = formatPValue(p)
if isnan(p)
    s = '= NA';
elseif p < 0.001
    s = '< 0.001';
else
    s = sprintf('= %.3f', p);
end
end