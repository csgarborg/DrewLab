function [videoTable, eventTable] = analyzeREMEventsTiming(excelPath, binSizeSeconds, useActualTimestamps, outputDir, durationBinSizeSeconds, rasterResolutionSeconds)
%ANALYZEREMEVENTS  Analyze REM event timing from a two-sheet Excel file.
%
%   [videoTable, eventTable] = analyzeREMEvents(excelPath, binSizeSeconds,
%                                                useActualTimestamps, outputDir)
%
%   Reads a two-sheet Excel workbook:
%     Sheet 1 (first column, header optional): all video file paths, e.g.
%         E:\SGNMData\260622\260622_10_43_3522.mp4
%     Sheet 2 (header row expected): REM events per recording, e.g.
%         file name | start1_s | stop1_s | ... | start5_s | stop5_s
%         E:\SGNMData\260624\260624_14_00_2424.tdms | 1180 | 1365 | ...
%
%   Filenames are assumed to encode: yymmdd_hh_mm_ss[xx]
%     - yymmdd   : recording date (20yy assumed)
%     - hh_mm_ss : the recording's actual start time
%     - trailing 2 digits (if present) are a sub-second/sequence tag,
%       ignored for timing purposes.
%
%   For each calendar day, all videos on Sheet 1 are chained together in
%   start-time order (up to 5 videos/day; not all have REM events, so not
%   all appear on Sheet 2 -- but Sheet 1 has the complete chain).
%
%   Each REM event's absolute time (seconds since that day's recording
%   session started) is:
%
%       abs_time = offset_of_its_video_within_the_day + start/stop_s_in_video
%
%   The video's offset within the day can be computed two ways:
%     - actual  (default): (this video's start datetime - first video's
%                start datetime). This automatically captures the ~25-30s
%                save gap between files, since it's measured directly from
%                the filename timestamps rather than assumed.
%     - nominal : (chain_index - 1) * 3600  [assumes exactly 1 hr/video]
%
%   Seven figures are produced (left open as normal MATLAB figures, not
%   saved as image files):
%     1. Histogram of REM event start times (seconds since the day's
%        recording session began), with a user-chosen bin size (seconds).
%        [remColor]
%     2. Scatter of REM event duration vs. time since the *start* of the
%        previous REM event, computed within each day's chain, skipping
%        the first event of each day (it has no "previous" event). [remColor]
%     3. Same as #2, but x-axis is time since the *end* (stop) of the
%        previous REM event instead of its start. [remColor]
%     4. Histogram of REM event durations, with a user-chosen bin size
%        (seconds). [remColor]
%     5. Scatter of REM event duration vs. time *until* the *start* of the
%        next REM event, skipping the last event of each day (it has no
%        "next" event). Tests whether duration predicts the upcoming gap
%        (the reverse of #2). [remColor]
%     6. Same as #5, but the gap is measured from this event's *end*
%        (stop) until the next event's start. [remColor]
%     7. Raster plot: one row per recording day *that has at least one
%        REM event* (days with none are omitted), x-axis = time since
%        that day's recording session started (0-5 hours), colored black
%        where there is no recording, blue where a video is recording but
%        no REM event is occurring, and yellow where a REM event is
%        occurring. Built with image() (true-color), not imagesc.
%
%   remColor = [0.85 0.33 0.10] is used for figures 1-4. Figure 7 uses its
%   own fixed black/blue/yellow scheme (unrelated to remColor).
%
%   INPUTS
%     excelPath              - path to the .xlsx file
%     binSizeSeconds          - start-time histogram bin width in seconds
%                                (default 300)
%     useActualTimestamps     - true/false (default true, see above)
%     outputDir                - folder for the output CSV (default '.')
%     durationBinSizeSeconds  - duration histogram bin width in seconds
%                                (default 30)
%     rasterResolutionSeconds - time-bin width (s) for the raster plot
%                                (default 5)
%
%   OUTPUTS
%     videoTable  - parsed, per-day-chained video list
%     eventTable  - one row per REM event with computed timing
%
%   Two CSVs are written to outputDir:
%     '<excel base name>_eventTIming_exact.csv' or '_simple.csv'
%         (per-event timing; suffix depends on useActualTimestamps)
%     '<excel base name>_videoChain.csv'
%         (the video chain/offsets, needed by aggregateREMEvents to build
%         its pooled raster plot)
%
%   Figures 2, 3, 5, and 6 (the four duration-vs-gap scatter plots) each
%   get an overlaid least-squares linear fit with its equation, R^2, and
%   p-value (for the slope) annotated on the plot.
%
%   Example:
%     [v, e] = analyzeREMEvents('my_data.xlsx', 300, true, '.');

if nargin < 6 || isempty(rasterResolutionSeconds), rasterResolutionSeconds = 5;  end
if nargin < 5 || isempty(durationBinSizeSeconds), durationBinSizeSeconds = 30; end
if nargin < 4 || isempty(outputDir),           outputDir = '.';       end
if nargin < 3 || isempty(useActualTimestamps), useActualTimestamps = true; end
if nargin < 2 || isempty(binSizeSeconds),      binSizeSeconds = 300;  end

remColor = [0.85 0.33 0.10];

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% ---- Sheet 1: video paths ----
raw1 = readcell(excelPath, 'Sheet', 1);
videoTable = buildVideoChain(raw1(:,1));

% ---- Sheet 2: REM events (assumes a header row; it is dropped) ----
raw2 = readcell(excelPath, 'Sheet', 2);
raw2(1,:) = [];
eventTable = buildEvents(raw2, videoTable, useActualTimestamps);

% ---- Plots (kept as MATLAB figures, not saved as images) ----
plotHistogramFcn(eventTable, binSizeSeconds, remColor);
plotDurationVsGapFcn(eventTable, remColor);
plotDurationVsGapEndFcn(eventTable, remColor);
plotDurationHistogramFcn(eventTable, durationBinSizeSeconds, remColor);
plotDurationVsNextStartFcn(eventTable, remColor);
plotDurationVsNextEndFcn(eventTable, remColor);
plotDayRasterFcn(videoTable, eventTable, useActualTimestamps, rasterResolutionSeconds);
plotDurationVsStartFcn(eventTable, remColor);

% ---- CSVs ----
[~, baseName, ~] = fileparts(excelPath);
if useActualTimestamps
    suffix = '_eventTIming_exact';
else
    suffix = '_eventTIming_simple';
end
csvPath = fullfile(outputDir, [baseName suffix '.csv']);
writetable(eventTable, csvPath);

videoChainPath = fullfile(outputDir, [baseName '_videoChain.csv']);
writetable(videoTable, videoChainPath);

fprintf('Parsed %d videos across %d day(s).\n', height(videoTable), numel(unique(videoTable.date)));
fprintf('Matched %d REM events.\n', height(eventTable));
fprintf('Saved CSVs to:\n  %s\n  %s\n', csvPath, videoChainPath);

end


% ===================== local functions =====================

function info = parseFilename(pathStr)
% Extract date/time info from a filename or path. info.valid = false if
% it doesn't match the expected pattern (e.g. a header cell).
info = struct('name', '', 'date', NaT, 'dt', NaT, 'xx', '', 'valid', false);

if ~(ischar(pathStr) || isstring(pathStr))
    return
end
pathStr = char(pathStr);
if isempty(strtrim(pathStr))
    return
end

parts = regexp(strtrim(pathStr), '[\\/]+', 'split');
name = parts{end};

tok = regexp(name, '(\d{6})_(\d{2})_(\d{2})_(\d{2})(\d{2})?', 'tokens', 'once');
if isempty(tok)
    return
end

dateStr = tok{1}; hh = tok{2}; mm = tok{3}; ss = tok{4};
if numel(tok) >= 5 && ~isempty(tok{5})
    xx = tok{5};
else
    xx = '';
end

yy = str2double(dateStr(1:2));
mo = str2double(dateStr(3:4));
dd = str2double(dateStr(5:6));

try
    dt = datetime(2000 + yy, mo, dd, str2double(hh), str2double(mm), str2double(ss));
catch
    return
end

info.name  = name;
info.date  = dateshift(dt, 'start', 'day');
info.dt    = dt;
info.xx    = xx;
info.valid = true;
end


function vt = buildVideoChain(pathCol)
% Parse Sheet 1 paths into a per-day, time-ordered chain with offsets.
n = numel(pathCol);
pathOut = strings(n,1);
nameOut = strings(n,1);
dateOut = NaT(n,1);
dtOut   = NaT(n,1);
xxOut   = strings(n,1);
valid   = false(n,1);

for i = 1:n
    info = parseFilename(pathCol{i});
    if info.valid
        valid(i)   = true;
        pathOut(i) = string(pathCol{i});
        nameOut(i) = info.name;
        dateOut(i) = info.date;
        dtOut(i)   = info.dt;
        xxOut(i)   = string(info.xx);
    end
end

if ~any(valid)
    error('No parseable video filenames found on Sheet 1.');
end

vt = table(pathOut(valid), nameOut(valid), dateOut(valid), dtOut(valid), xxOut(valid), ...
    'VariableNames', {'path','name','date','dt','xx'});
vt = sortrows(vt, {'date','dt'});

[udays, ~, ic] = unique(vt.date);
chainIndex    = zeros(height(vt),1);
offsetActual  = zeros(height(vt),1);
offsetNominal = zeros(height(vt),1);

for d = 1:numel(udays)
    idx = find(ic == d);          % already in dt order within the day
    dayStart = vt.dt(idx(1));
    for k = 1:numel(idx)
        chainIndex(idx(k))    = k;
        offsetActual(idx(k))  = seconds(vt.dt(idx(k)) - dayStart);
        offsetNominal(idx(k)) = (k-1) * 3600;
    end
end

vt.chain_index      = chainIndex;
vt.offset_actual_s  = offsetActual;
vt.offset_nominal_s = offsetNominal;
end


function vrow = matchVideo(vt, pathStr)
% Find the Sheet-1 video row corresponding to a Sheet-2 filename.
vrow = [];
info = parseFilename(pathStr);
if ~info.valid
    return
end

mask = (vt.date == info.date) & ...
       (hour(vt.dt) == hour(info.dt)) & ...
       (minute(vt.dt) == minute(info.dt)) & ...
       (second(vt.dt) == second(info.dt));

if ~any(mask)
    return
end

maskExact = mask & (vt.xx == string(info.xx));
if any(maskExact)
    idx = find(maskExact, 1);
else
    idx = find(mask, 1);
end
vrow = vt(idx,:);
end


function tf = isNumericValid(v)
tf = isnumeric(v) && isscalar(v) && ~isnan(v);
end


function et = buildEvents(raw2, vt, useActual)
nRows = size(raw2, 1);

dayCol     = NaT(0,1);
videoPath  = strings(0,1);
chainIdx   = zeros(0,1);
eventNum   = zeros(0,1);
startInVid = zeros(0,1);
stopInVid  = zeros(0,1);
absStart   = zeros(0,1);
absStop    = zeros(0,1);
dur        = zeros(0,1);

unmatchedCount = 0;

for r = 1:nRows
    fileCell = raw2{r,1};
    vrow = matchVideo(vt, fileCell);
    if isempty(vrow)
        if ischar(fileCell) || isstring(fileCell)
            unmatchedCount = unmatchedCount + 1;
        end
        continue
    end

    if useActual
        offset = vrow.offset_actual_s;
    else
        offset = vrow.offset_nominal_s;
    end

    for i = 1:5
        sCol = 2 + (i-1)*2;   % 1=file, 2/3=start1/stop1, 4/5=start2/stop2, ...
        eCol = sCol + 1;
        if eCol > size(raw2, 2)
            break
        end
        sVal = raw2{r, sCol};
        eVal = raw2{r, eCol};
        if ~isNumericValid(sVal) || ~isNumericValid(eVal)
            continue
        end
        s = double(sVal);
        e = double(eVal);

        dayCol(end+1,1)     = vrow.date;      %#ok<AGROW>
        videoPath(end+1,1)  = vrow.path;      %#ok<AGROW>
        chainIdx(end+1,1)   = vrow.chain_index; %#ok<AGROW>
        eventNum(end+1,1)   = i;              %#ok<AGROW>
        startInVid(end+1,1) = s;              %#ok<AGROW>
        stopInVid(end+1,1)  = e;              %#ok<AGROW>
        absStart(end+1,1)   = offset + s;     %#ok<AGROW>
        absStop(end+1,1)    = offset + e;     %#ok<AGROW>
        dur(end+1,1)        = e - s;          %#ok<AGROW>
    end
end

if unmatchedCount > 0
    warning('%d row(s) on Sheet 2 could not be matched to a video on Sheet 1 and were skipped.', unmatchedCount);
end
if isempty(dayCol)
    error('No REM events could be matched/parsed from Sheet 2.');
end

et = table(dayCol, videoPath, chainIdx, eventNum, startInVid, stopInVid, absStart, absStop, dur, ...
    'VariableNames', {'day','video_path','chain_index','event_num', ...
                       'start_s_in_video','stop_s_in_video', ...
                       'abs_start_s','abs_stop_s','duration_s'});
et = sortrows(et, {'day','abs_start_s'});

et.gap_from_prev_s = nan(height(et),1);
et.gap_from_prev_end_s = nan(height(et),1);
et.gap_to_next_start_s = nan(height(et),1);
et.gap_to_next_end_s = nan(height(et),1);
[udays, ~, ic] = unique(et.day);
for d = 1:numel(udays)
    idx = find(ic == d);
    starts = et.abs_start_s(idx);
    stops  = et.abs_stop_s(idx);
    gaps = [NaN; diff(starts)];
    gapsFromEnd = [NaN; starts(2:end) - stops(1:end-1)];
    gapsToNextStart = [diff(starts); NaN];
    gapsToNextEnd = [starts(2:end) - stops(1:end-1); NaN];
    et.gap_from_prev_s(idx) = gaps;
    et.gap_from_prev_end_s(idx) = gapsFromEnd;
    et.gap_to_next_start_s(idx) = gapsToNextStart;
    et.gap_to_next_end_s(idx) = gapsToNextEnd;
end
end


function fig = plotHistogramFcn(et, binSize, remColor)
maxT = max(et.abs_start_s);
edges = 0:binSize:(maxT + binSize);

fig = figure('Name', 'REM event histogram');
subplot(1,2,1)
histogram(et.abs_start_s, edges, 'EdgeColor', 'black', 'FaceColor', remColor);
xlabel(sprintf('Time since recording session start (s) \x2014 %g s bins', binSize));
ylabel('REM event count');
title('REM events by time since start of recording session');

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
title('REM event PDF by time since start of recording session');

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

fig = figure('Name', 'REM duration vs. gap from previous start');
scatter(et.gap_from_prev_s(mask), et.duration_s(mask), 36, remColor, 'filled', ...
    'MarkerFaceAlpha', 0.6);
xlabel('Time since previous REM event start (s)');
ylabel('REM event duration (s)');
title({'REM event duration vs. time since previous REM event start', ...
       '(first event of each day excluded)'});
addFitLine(et.gap_from_prev_s(mask), et.duration_s(mask));
end


function fig = plotDurationVsGapEndFcn(et, remColor)
mask = ~isnan(et.gap_from_prev_end_s);

fig = figure('Name', 'REM duration vs. gap from previous end');
scatter(et.gap_from_prev_end_s(mask), et.duration_s(mask), 36, remColor, 'filled', ...
    'MarkerFaceAlpha', 0.6);
xlabel('Time since previous REM event end (s)');
ylabel('REM event duration (s)');
title({'REM event duration vs. time since previous REM event end', ...
       '(first event of each day excluded)'});
addFitLine(et.gap_from_prev_end_s(mask), et.duration_s(mask));
end


function fig = plotDurationHistogramFcn(et, binSize, remColor)
maxD = max(et.duration_s);
edges = 0:binSize:(maxD + binSize);

fig = figure('Name', 'REM event duration histogram');
subplot(1,2,1)
histogram(et.duration_s, edges, 'EdgeColor', 'black', 'FaceColor', remColor);
xlabel(sprintf('REM event duration (s) \x2014 %g s bins', binSize));
ylabel('REM event count');
title('REM event durations');

subplot(1,2,2)
histogram(et.duration_s, edges, 'Normalization', 'pdf', 'EdgeColor', 'black', 'FaceColor', remColor, 'DisplayName', 'PDF');
xlabel(sprintf('REM event duration (s) \x2014 %g s bins', binSize));
ylabel('REM event PDF');
title('REM event durations');

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

fig = figure('Name', 'REM duration vs. gap until next start');
scatter(et.duration_s(mask), et.gap_to_next_start_s(mask), 36, remColor, 'filled', ...
    'MarkerFaceAlpha', 0.6);
ylabel('Time until next REM event start (s)');
xlabel('REM event duration (s)');
title({'REM event duration vs. time until next REM event start', ...
       '(last event of each day excluded)'});
addFitLine(et.duration_s(mask), et.gap_to_next_start_s(mask));
end


function fig = plotDurationVsNextEndFcn(et, remColor)
mask = ~isnan(et.gap_to_next_end_s);

fig = figure('Name', 'REM duration vs. gap (end) until next start');
scatter(et.duration_s(mask), et.gap_to_next_end_s(mask), 36, remColor, 'filled', ...
    'MarkerFaceAlpha', 0.6);
ylabel('Time from this event''s end until next REM event start (s)');
xlabel('REM event duration (s)');
title({'REM event duration vs. time until next REM event', ...
       '(measured from this event''s end; last event of each day excluded)'});
addFitLine(et.duration_s(mask), et.gap_to_next_end_s(mask));
end


function fig = plotDayRasterFcn(vt, et, useActual, resolutionSeconds)
% One row per recording day: black = no recording, blue = recording with
% no REM, yellow = REM event. Assumes each video covers a nominal 3600s
% starting at its offset (see header comment on offset computation).
maxSeconds = 5 * 3600;

if useActual
    offsetCol = 'offset_actual_s';
else
    offsetCol = 'offset_nominal_s';
end

days = unique(vt.date, 'stable'); % vt is already sorted chronologically
hasREM = ismember(days, unique(et.day));
days = days(hasREM);
nDays = numel(days);

edges = 0:resolutionSeconds:maxSeconds;
nBins = numel(edges) - 1;
binCenters = edges(1:end-1) + resolutionSeconds/2;

catMat = zeros(nDays, nBins);
for d = 1:nDays
    dayVideos = vt(vt.date == days(d), :);
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

    catMat(d,:) = row;
end

rgbImg = categoriesToRGB(catMat);

fig = figure('Name', 'REM raster by day');
image(binCenters / 3600, 1:nDays, rgbImg);
set(gca, 'YDir', 'normal');
yticks(1:nDays);
yticklabels(string(days, 'yyyy-MM-dd'));
xlabel('Time since recording session start (hours)');
ylabel('Day');
title('REM events across recording days');
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