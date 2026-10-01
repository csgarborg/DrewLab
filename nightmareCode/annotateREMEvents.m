function annotateREMEvents(xlsxFile, varargin)
%ANNOTATEREMEVENTS  Step through candidate REM events listed in an Excel file
%   and log the start/stop time of each one.
%
%   annotateREMEvents('C:\data\REM_candidates.xlsx')
%   annotateREMEvents(xlsxFile, 'ProcDataDir','D:\ProcData', 'RedoSkipped',true)
%
% EXCEL LAYOUT
%   Sheet 1 : column A = path to .mp4 file. Columns B, C, ... = approximate
%             REM timestamps in MINUTES. "x" = no REM event. Anything after
%             the first 1-2 digits in a cell (notes) is ignored.
%   Sheet 2 : file name | start1_s | stop1_s | start2_s | stop2_s | ...
%             (file name = mp4 name with .tdms extension, times in whole s).
%             Created automatically if missing. Events are kept sorted by
%             start time within each row.
%
% For every timestamp that has not been handled yet, generateSleepScorePlotRL
% is called for [t-PreMin, t+PostMin] (clipped to the recording) and you can:
%   * click START then STOP on the plot (a third click starts over),
%     then press Enter / "Save"  -> times go to sheet 2, sheet 1 cell = GREEN
%   * close the figure / press "Skip" -> sheet 1 cell = RED, nothing saved
%   * "Extend start / end" -> replot with that side widened by ExtendMin
%   * "Quit session" -> stop now, current cell left unmarked
% Cells already GREEN (or RED, unless 'RedoSkipped' is true) are skipped, so
% you can add new timestamps to sheet 1 later and just rerun.
%
% NAME-VALUE OPTIONS
%   'PreMin'        2      minutes shown before the timestamp
%   'PostMin'       5      minutes shown after the timestamp
%   'ExtendMin'     1      minutes added by each extend button
%   'ProcDataDir'   ''     folder holding ProcData .mat files ('' = same
%                          folder as the mp4)
%   'ProcDataFcn'   []     @(mp4Path) -> ProcData .mat path (overrides the
%                          default lookup)
%   'RedoSkipped'   false  also revisit RED cells
%   'DarkMode'      true   dark figure background
%   'ConfirmClicks' true   require Enter/Save after the two clicks
%
% REQUIREMENTS: Windows + desktop Excel (cell colours are read/written
% through COM), generateSleepScorePlotRL.m on the MATLAB path, and the
% workbook must NOT be open in Excel while this runs.

p = inputParser;
addParameter(p,'PreMin',2);
addParameter(p,'PostMin',5);
addParameter(p,'ExtendMin',1);
addParameter(p,'ProcDataDir','');
addParameter(p,'ProcDataFcn',[]);
addParameter(p,'RedoSkipped',false);
addParameter(p,'DarkMode',true);
addParameter(p,'ConfirmClicks',true);
parse(p,varargin{:});
o = p.Results;

if isempty(o.ProcDataFcn)
    procFcn = @(mp4) findProcData(mp4, o.ProcDataDir);
else
    procFcn = o.ProcDataFcn;
end
viewOpts = struct('PreSec',o.PreMin*60, 'PostSec',o.PostMin*60, ...
    'ExtendSec',o.ExtendMin*60, 'Dark',o.DarkMode, 'Confirm',o.ConfirmClicks);

GREEN = xlrgb(146,208,80);
RED   = xlrgb(255,80,80);

%% Open workbook through COM
if ~ispc
    error('annotateREMEvents:platform', ...
        'Reading/writing cell colours needs Windows with desktop Excel (COM).');
end
[ok, att] = fileattrib(xlsxFile);
if ~ok, error('annotateREMEvents:nofile','File not found: %s', xlsxFile); end
xlsxFile = att.Name;

xl = actxserver('Excel.Application');
cleaner = onCleanup(@() shutdownExcel(xl)); %#ok<NASGU>
xl.Visible = false;
xl.DisplayAlerts = false;
wb = xl.Workbooks.Open(xlsxFile);
if wb.ReadOnly
    error('annotateREMEvents:readonly', ...
        'Workbook opened read-only. Close it in Excel and try again.');
end

%% Read sheet 1
sh1 = wb.Sheets.Item(1);
ur = sh1.UsedRange;
lastRow = ur.Row + ur.Rows.Count - 1;
lastCol = ur.Column + ur.Columns.Count - 1;
raw = sh1.Range(rangeAddr(1,1,lastRow,lastCol)).Value;
if ~iscell(raw), raw = {raw}; end

%% Sheet 2: create if needed, read existing results
if wb.Sheets.Count < 2
    sh2 = wb.Sheets.Add([], wb.Sheets.Item(1));
else
    sh2 = wb.Sheets.Item(2);
end
if isMissingCell(sh2.Range('A1').Value)
    sh2.Range('A1').Value = 'file name';
end
ur2 = sh2.UsedRange;
lr2 = ur2.Row + ur2.Rows.Count - 1;
lc2 = ur2.Column + ur2.Columns.Count - 1;
raw2 = sh2.Range(rangeAddr(1,1,lr2,lc2)).Value;
if ~iscell(raw2), raw2 = {raw2}; end

rowOf   = containers.Map('KeyType','char','ValueType','double');
pairsOf = containers.Map('KeyType','char','ValueType','any');
for r = 2:size(raw2,1)
    nm = raw2{r,1};
    if ~ischar(nm) || isempty(strtrim(nm)), continue; end
    key = lower(strtrim(nm));
    v = nan(1, max(size(raw2,2)-1, 0));
    for k = 2:size(raw2,2)
        c = raw2{r,k};
        if isnumeric(c) && isscalar(c), v(k-1) = c; end
    end
    if mod(numel(v),2), v(end+1) = NaN; end %#ok<AGROW>
    pr = [v(1:2:end)', v(2:2:end)'];
    pr = pr(~any(isnan(pr),2),:);
    rowOf(key)   = r;
    pairsOf(key) = pr;
end
nextRow2 = max(lr2,1) + 1;

%% Collect events that still need attention
ev = struct('row',{},'col',{},'mp4',{},'min',{});
nDone = 0;
for r = 1:lastRow
    a = raw{r,1};
    if ~ischar(a) || isempty(regexpi(strtrim(a),'\.mp4$','once')), continue; end
    mp4 = strtrim(a);
    for c = 2:lastCol
        cv = raw{r,c};
        if isMissingCell(cv), continue; end
        txt = cellText(cv);
        if isempty(txt) || ~isempty(regexpi(txt,'^\s*x','once')), continue; end
        tok = regexp(txt,'\d{1,2}','match','once');
        if isempty(tok)
            warning('annotateREMEvents:noDigits', ...
                'No timestamp found in %s ("%s") - ignored.', cellAddr(r,c), txt);
            continue
        end
        col = sh1.Range(cellAddr(r,c)).Interior.Color;
        if col == GREEN || (col == RED && ~o.RedoSkipped)
            nDone = nDone + 1;
            continue
        end
        ev(end+1) = struct('row',r,'col',c,'mp4',mp4,'min',str2double(tok)); %#ok<AGROW>
    end
end
fprintf('%d event(s) to review, %d already done.\n', numel(ev), nDone);

%% Review loop
cache  = containers.Map('KeyType','char','ValueType','any');
nSaved = 0; nSkipped = 0;
for i = 1:numel(ev)
    e = ev(i);
    addr = cellAddr(e.row, e.col);

    if ~isKey(cache, e.mp4)
        info = struct('proc','','dur',NaN);
        try
            info.proc = procFcn(e.mp4);
            if ~isempty(info.proc), info.dur = getTrialDur(info.proc); end
        catch ME
            warning('annotateREMEvents:procdata', ...
                'Could not load ProcData for %s: %s', e.mp4, ME.message);
            info.proc = '';
        end
        cache(e.mp4) = info;
    end
    info = cache(e.mp4);
    if isempty(info.proc)
        warning('annotateREMEvents:skipfile','No ProcData for %s - %s left unmarked.', e.mp4, addr);
        continue
    end

    tCenter = e.min * 60;
    if tCenter - viewOpts.PreSec >= info.dur
        warning('annotateREMEvents:beyond', ...
            '%s: %d min is past the end of the recording (%.1f min) - marked red.', ...
            addr, e.min, info.dur/60);
        sh1.Range(addr).Interior.Color = RED;
        wb.Save;
        nSkipped = nSkipped + 1;
        continue
    end

    [mp4Path, base] = fileparts(e.mp4);
    fprintf('[%d/%d] %s  ~%d min  (%s)\n', i, numel(ev), base, e.min, addr);
    titleStr = sprintf('[%d/%d] %s  |  ~%d min  |  %s', i, numel(ev), base, e.min, addr);

    [action, s, t] = reviewEvent(info.proc, tCenter, info.dur, titleStr, viewOpts);

    switch action
        case 'save'
            name = [fullfile(mp4Path,base) '.tdms'];
            key  = lower(name);
            if isKey(rowOf,key)
                rr = rowOf(key);
                pr = pairsOf(key);
            else
                rr = nextRow2; nextRow2 = nextRow2 + 1;
                pr = zeros(0,2);
                rowOf(key) = rr;
                sh2.Range(cellAddr(rr,1)).Value = name;
            end
            pr = sortrows([pr; s t]);
            pairsOf(key) = pr;
            n = size(pr,1);
            hdr = cell(1,2*n);
            for k = 1:n
                hdr{2*k-1} = sprintf('start%d_s',k);
                hdr{2*k}   = sprintf('stop%d_s',k);
            end
            sh2.Range(rangeAddr(1,2,1,1+2*n)).Value  = hdr;
            sh2.Range(rangeAddr(rr,2,rr,1+2*n)).Value = reshape(pr.',1,[]);
            sh1.Range(addr).Interior.Color = GREEN;
            wb.Save;
            nSaved = nSaved + 1;
            fprintf('    saved: start %d s, stop %d s\n', s, t);

        case 'skip'
            sh1.Range(addr).Interior.Color = RED;
            wb.Save;
            nSkipped = nSkipped + 1;
            fprintf('    skipped\n');

        case 'quit'
            fprintf('Session ended early.\n');
            break
    end
end

wb.Save;
fprintf('Done. Saved %d, skipped %d.\n', nSaved, nSkipped);
end


%% ========================================================================
function [action, startS, stopS] = reviewEvent(procFile, tCenter, trialDur, titleStr, o)
% Show the plot for one event and run the click / extend / skip interaction.

tStart = max(0, tCenter - o.PreSec);
tEnd   = min(trialDur, tCenter + o.PostSec);
clicks = [];
action = 'skip'; startS = NaN; stopS = NaN;
fig = []; axAll = []; ax1 = []; ax6 = [];
hLines = gobjects(0); hStatus = []; figPos = []; figState = 'normal'; state = '';
cleaner = onCleanup(@() delete(findall(groot,'Type','figure','Tag','REMEventReview'))); %#ok<NASGU>

while true
    state = '';
    [fig, ax1, ax2, ax3, ax4, ax5, ax6] = generateSleepScorePlotRL(procFile, tStart, tEnd);
    axAll = [ax1 ax2 ax3 ax4 ax5 ax6];
    delete(findobj(fig,'Type','scatter'));
    hLines = gobjects(0);

    set(fig, 'NumberTitle','off', 'Name',titleStr, 'Tag','REMEventReview', ...
        'CloseRequestFcn',    @(~,~) finish('skip'), ...
        'WindowButtonDownFcn',@onClick, ...
        'KeyPressFcn',        @onKey);
    if o.Dark, applyDarkTheme(fig, axAll); end
    buildControls();

    % dashed marker at the approximate timestamp
    for k = 1:numel(axAll)
        xline(axAll(k), tCenter, '--', 'Color',[0.7 0.7 0.7], 'LineWidth',0.75);
    end
    redrawClicks();
    updateStatus();

    if isempty(figPos)
        try fig.WindowState = 'maximized'; catch, end
    else
        fig.Position = figPos;
        try fig.WindowState = figState; catch, end
    end

    uiwait(fig);

    if ~isempty(fig) && isvalid(fig)
        figPos = fig.Position;
        try figState = fig.WindowState; catch, end
        delete(fig);
    end

    switch state
        case 'replot'
            continue
        case 'save'
            action = 'save';
            startS = round(min(clicks));
            stopS  = round(max(clicks));
            return
        case 'quit'
            action = 'quit';
            return
        otherwise
            action = 'skip';
            return
    end
end

    % ---------------------------------------------------------------- UI
    function buildControls()
        ext = o.ExtendSec/60;
        addButton(sprintf('Extend start -%g min',ext), [0.005 0.945 0.13 0.045], @(~,~) extendWindow('start'), tStart > 0);
        addButton(sprintf('Extend end +%g min',ext),   [0.140 0.945 0.13 0.045], @(~,~) extendWindow('end'),   tEnd < trialDur);
        addButton('Reset clicks (R)',  [0.275 0.945 0.11 0.045], @(~,~) resetClicks(), true);
        addButton('Save (Enter)',      [0.390 0.945 0.10 0.045], @(~,~) trySave(),     true);
        addButton('Skip (mark red)',   [0.495 0.945 0.11 0.045], @(~,~) finish('skip'),true);
        addButton('Quit session',      [0.610 0.945 0.09 0.045], @(~,~) finish('quit'),true);
        txtArgs = {};
        if o.Dark
            txtArgs = {'BackgroundColor',[0.10 0.10 0.10], 'ForegroundColor',[1 0.9 0.4]};
        end
        hStatus = uicontrol(fig, 'Style','text', 'Units','normalized', ...
            'Position',[0.705 0.945 0.29 0.045], 'HorizontalAlignment','left', ...
            'FontSize',10, txtArgs{:});
    end

    function addButton(label, pos, cb, enabled)
        en = 'on';
        if ~enabled, en = 'off'; end
        btnArgs = {};
        if o.Dark
            btnArgs = {'BackgroundColor',[0.25 0.25 0.25], 'ForegroundColor',[0.95 0.95 0.95]};
        end
        uicontrol(fig, 'Style','pushbutton', 'String',label, 'Units','normalized', ...
            'Position',pos, 'Callback',cb, 'Enable',en, 'FontSize',10, btnArgs{:});
    end

    function updateStatus(msg)
        if nargin == 0
            switch numel(clicks)
                case 0
                    msg = 'Click REM START, then STOP';
                case 1
                    msg = sprintf('Mark at %.0f s - click the second point', clicks);
                otherwise
                    msg = sprintf('Start %d s | Stop %d s  (Enter = save)', ...
                        round(min(clicks)), round(max(clicks)));
            end
        end
        if ~isempty(hStatus) && isgraphics(hStatus), hStatus.String = msg; end
    end

    function redrawClicks()
        delete(hLines(isvalid(hLines)));
        hLines = gobjects(0);
        n = numel(clicks);
        if n == 0, return; end
        if n == 1
            xs = clicks;       cols = [1 0.85 0.2];
        else
            xs = sort(clicks); cols = [0.3 1 0.3; 1 0.35 0.35];
        end
        for m = 1:numel(xs)
            for k = 1:numel(axAll)
                hLines(end+1) = xline(axAll(k), xs(m), '-', ...
                    'Color',cols(m,:), 'LineWidth',1.5); %#ok<AGROW>
            end
        end
    end

    % -------------------------------------------------------- callbacks
    function onClick(~,~)
        if ~strcmp(fig.SelectionType,'normal'), return; end
        pT = getpixelposition(ax1,true);      % top axes
        pB = getpixelposition(ax6,true);      % bottom (spectrogram) axes
        cp = fig.CurrentPoint;
        if cp(1) < pT(1) || cp(1) > pT(1)+pT(3) || ...
           cp(2) < pB(2) || cp(2) > pT(2)+pT(4)
            return
        end
        x = ax1.CurrentPoint(1,1);            % all axes share the same x-range
        x = min(max(x, tStart), tEnd);
        if numel(clicks) >= 2
            clicks = x;                       % third click starts over
        else
            clicks(end+1) = x;
        end
        redrawClicks();
        updateStatus();
        if numel(clicks) == 2 && ~o.Confirm
            trySave();
        end
    end

    function onKey(~,evt)
        switch evt.Key
            case {'return','enter'}
                trySave();
            case 'r'
                resetClicks();
        end
    end

    function resetClicks()
        clicks = [];
        redrawClicks();
        updateStatus();
    end

    function trySave()
        if numel(clicks) < 2
            updateStatus('Need two clicks (start and stop) before saving');
            return
        end
        if round(max(clicks)) <= round(min(clicks))
            resetClicks();
            updateStatus('Start and stop round to the same second - click again');
            return
        end
        finish('save');
    end

    function extendWindow(which)
        if strcmp(which,'start')
            tStart = max(0, tStart - o.ExtendSec);
        else
            tEnd = min(trialDur, tEnd + o.ExtendSec);
        end
        finish('replot');
    end

    function finish(s)
        state = s;
        uiresume(fig);
    end

end


%% ========================================================================
function applyDarkTheme(fig, axList)
bg = [0.10 0.10 0.10];
fg = [0.92 0.92 0.92];
fig.Color = bg;
for k = 1:numel(axList)
    ax = axList(k);
    ax.Color     = [0.13 0.13 0.13];
    ax.XColor    = fg;
    ax.YColor    = fg;
    ax.GridColor = fg;
    ax.GridAlpha = 0.25;
    ax.XLabel.Color = fg;
    ax.YLabel.Color = fg;
    ax.Title.Color  = fg;
    set(findobj(ax,'Type','text'), 'Color', fg);
    ln = findobj(ax,'Type','line');
    for j = 1:numel(ln)
        ln(j).Color = lightenColor(ln(j).Color);
    end
    sc = findobj(ax,'Type','scatter');
    for j = 1:numel(sc)
        set(sc(j), 'MarkerFaceColor',[1 0.75 0.3], 'MarkerEdgeColor',[1 0.75 0.3]);
    end
end
end

function c = lightenColor(c)
% Brighten traces that would be hard to see on a dark background.
lum = 0.2126*c(1) + 0.7152*c(2) + 0.0722*c(3);
if lum < 0.45
    c = c + (1 - c)*0.55;
end
end


%% ========================================================================
function procFile = findProcData(mp4Path, procDir)
% Default mp4 -> ProcData lookup. EDIT HERE (or pass 'ProcDataFcn') if your
% ProcData files are named differently.
[folder, base] = fileparts(mp4Path);
if ~isempty(procDir), folder = procDir; end
procFile = '';
cand = fullfile(folder, [base '_ProcData.mat']);
if isfile(cand)
    procFile = cand;
    return
end
d = dir(fullfile(folder, [base '*ProcData*.mat']));
if ~isempty(d)
    procFile = fullfile(d(1).folder, d(1).name);
    return
end
[f, pth] = uigetfile(fullfile(folder,'*ProcData*.mat'), ...
    sprintf('Locate ProcData file for %s', base));
if ~isequal(f,0), procFile = fullfile(pth,f); end
end

function dur = getTrialDur(procFile)
% Recording length in seconds, computed the same way as generateSleepScorePlotRL.
S = load(procFile,'ProcData');
P = S.ProcData;
Fs = 60;
if isfield(P,'notes') && isfield(P.notes,'dsFs') && ~isempty(P.notes.dsFs) && isfinite(P.notes.dsFs)
    Fs = P.notes.dsFs;
end
if isfield(P,'ECoG_DS_norm')
    n = numel(P.ECoG_DS_norm);
elseif isfield(P,'ECoG_DS')
    n = numel(P.ECoG_DS);
else
    error('Cannot determine trial duration.');
end
dur = n / Fs;
end


%% ========================================================================
function tf = isMissingCell(c)
tf = isempty(c) || (isnumeric(c) && isscalar(c) && isnan(c));
end

function txt = cellText(c)
if ischar(c)
    txt = strtrim(c);
elseif isnumeric(c) && isscalar(c) && isfinite(c)
    txt = sprintf('%g', c);
else
    txt = '';
end
end

function v = xlrgb(r,g,b)
% Excel stores colours as BGR integers.
v = r + 256*g + 65536*b;
end

function s = colLetter(n)
s = '';
while n > 0
    m = mod(n-1,26);
    s = [char(65+m) s]; %#ok<AGROW>
    n = floor((n-1)/26);
end
end

function a = cellAddr(r,c)
a = sprintf('%s%d', colLetter(c), r);
end

function a = rangeAddr(r1,c1,r2,c2)
a = [cellAddr(r1,c1) ':' cellAddr(r2,c2)];
end

function shutdownExcel(xl)
try xl.DisplayAlerts = false; xl.Quit; catch, end
try delete(xl); catch, end
end