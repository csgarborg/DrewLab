function [figHandle,ax1,ax2,ax3,ax4,ax6] = generateSleepScorePlot(procDataFileID,tStart,tEnd)

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
tPlot = tVec(tMask);

%% Force
if isfield(ProcData,'forceSensor_norm')
    force = ProcData.forceSensor_norm(:);
    forceLabel = 'Force (norm)';
elseif isfield(ProcData,'forceSensor')
    force = ProcData.forceSensor(:);
    forceLabel = 'Force (V)';
else
    force = zeros(refLen,1);
    forceLabel = 'Force';
end
force = force(tMask);

%% Binary force
if isfield(ProcData,'binForceSensor')
    binForce = logical(ProcData.binForceSensor(:));
else
    binForce = false(refLen,1);
end
binForce = binForce(tMask);

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
emgPower = emgPower(tMask);

%% Raw EMG (FIXED)
rawEMG = [];
if isfield(ProcData,'EMG') && isfield(ProcData.EMG,'emgSignal')
    rawEMG = ProcData.EMG.emgSignal(:);
end

if isempty(rawEMG)
    rawEMG = zeros(refLen,1);
elseif numel(rawEMG) < refLen
    rawEMG(end+1:refLen) = rawEMG(end);
elseif numel(rawEMG) > refLen
    rawEMG = rawEMG(1:refLen);
end
rawEMG = rawEMG(tMask);

%% Load spectrogram
[folder, base, ~] = fileparts(procDataFileID);
prefix = regexprep(base,'_ProcData$','');

specFile = '';
cand = dir(fullfile(folder, [prefix '*Spec*.mat']));
if ~isempty(cand)
    specFile = fullfile(folder, cand(1).name);
end

Sspec = []; Fspec = []; Tspec = [];
isNormSpec = false;

if ~isempty(specFile) && isfile(specFile)
    L = load(specFile);

    if isfield(L,'SpecData') && isfield(L.SpecData,'ECoG')
        SD = L.SpecData.ECoG;
        if isfield(SD,'normS')
            Sspec = SD.normS;
            Fspec = SD.F;
            Tspec = SD.T;
            isNormSpec = true;
        elseif all(isfield(SD,{'S','F','T'}))
            Sspec = SD.S;
            Fspec = SD.F;
            Tspec = SD.T;
        end
    end
end

%% --- SPECTROGRAM SAFETY FIXES ---

if ~isempty(Sspec) && ~isempty(Fspec) && ~isempty(Tspec)

    % Remove invalid freqs for log scale
    validIdx = Fspec > 0;
    Fspec = Fspec(validIdx);
    Sspec = Sspec(validIdx,:);

    % Fix NaN/Inf/negative values
    Sspec(~isfinite(Sspec)) = eps;
    Sspec(Sspec <= 0) = eps;

    % Dimension check
    if size(Sspec,1) ~= numel(Fspec) || size(Sspec,2) ~= numel(Tspec)
        warning('Spectrogram dimension mismatch. Skipping.');
        Sspec = [];
    end
end

%% Build figure
figHandle = figure;

%% Force
ax1 = subplot(5,1,1);
plot(tPlot, force,'LineWidth',1);
ylabel(forceLabel);
xlim([0 trialDur]);
set(gca,'XTickLabel',[]);
ax1.XGrid = 'on'; 
ax1.YGrid = 'off';
xlim([tStart tEnd])

%% EMG power + binForce
ax2 = subplot(5,1,2); hold on;
plot(tPlot, emgPower,'r');
ax2.XGrid = 'on'; 
ax2.YGrid = 'off';

if any(binForce)
    yMax = max(emgPower);
    ySpan = max(emgPower)-min(emgPower);
    if ySpan <= 0, ySpan = 1; end
    scatter(tPlot(binForce), yMax + 0.05*ySpan,'filled');
end

ylabel(emgLabel);
xlim([0 trialDur]);
set(gca,'XTickLabel',[]);

thetaBand = Fspec >= 6 & Fspec <= 10;
thetaPower = mean(Sspec(thetaBand,:),1);

hold(ax2,'on')
plot(ax2, Tspec, thetaPower/max(thetaPower)*max(emgPower), 'g')
xlim([tStart tEnd])

%% Raw EMG
ax3 = subplot(5,1,3);
plot(tPlot, rawEMG);
ylabel('EMG raw');
xlim([0 trialDur]);
set(gca,'XTickLabel',[]);
ax3.XGrid = 'on'; 
ax3.YGrid = 'off';
xlim([tStart tEnd])

%% Spectrogram
axSpec = subplot(5,1,4:5);
ax6 = axSpec;

tMaskSpec = Tspec >= tStart & Tspec <= tEnd;
Tspec = Tspec(tMaskSpec);
Sspec = Sspec(:, tMaskSpec);

if ~isempty(Sspec)

    if isNormSpec
        imagesc(Tspec, Fspec, 100*Sspec);
        % clim([-100 100]);
        Splot = 100*Sspec;
        p = prctile(Splot(:), [2 98]);  % or Sdb if using dB
        if all(isfinite(p))
            clim(p);
        end
        % colorbar;
    else
        Sdb = 10*log10(Sspec);
        imagesc(Tspec, Fspec, Sdb);
        p = prctile(Sdb(:),[2 98]);
        if all(isfinite(p)), clim(p); end
        % colorbar;
    end

    axis xy;
    set(gca,'YScale','log'); % safe now
    ylabel('Freq (Hz)');
    xlabel('Time (s)');
    xlim([0 trialDur]);

else
    text(0.5,0.5,'No spectrogram','Units','normalized');
end
xlim([tStart tEnd])

% axList = [ax1, ax2, ax3, axSpec];

% setTimeWindow(axList, tStart, tEnd)  % look at 100–120 sec

%% IMPORTANT: disable linkaxes (can re-enable later)
% linkaxes([ax1,ax2,ax3,axSpec],'x');

ax4 = ax2;

end

% function setTimeWindow(axList, tStart, tEnd)
%     for i = 1:length(axList)
%         if isvalid(axList(i))
%             xlim(axList(i), [tStart tEnd]);
%         end
%     end
% end

% function [figHandle,ax1,ax2,ax3,ax4,ax6] = generateSleepScorePlot(procDataFileID,saveFigs)
% 
% %   generateSleepScorePlotData
% %   Build a per-trial summary figure for manual sleep scoring for KECK analog data.
% %   Uses:
% %       • forceSensor_norm (if present) or forceSensor
% %       • binForceSensor (now overlaid on EMG subplot)
% %       • EMG.emgPower_norm (if present) or EMG.emgPower
% %       • ECoG_DS_norm (if present) or ECoG_DS
% %       • pre-computed ECoG spectrogram (SpecData.ECoG.normS if present, else S)
% %
% %   [figHandle,ax1,ax2,ax3,ax4,ax6] = generateSleepScorePlotData(procDataFileID, saveFigs)
% %
% %   RETURNS
% %       figHandle : figure handle
% %       ax1       : force axis
% %       ax2       : EMG axis (normalized if available) + binForce overlay
% %       ax3       : ECoG DS axis (normalized if available)
% %       ax4       : EMG axis (duplicate of ax2; used by scoring code for xlines)
% %       ax6       : spectrogram axis (used by scoring code for xlines/view)
% %
% %   NOTE
% %       This function is designed to be API-compatible with your
% %       CreateTrainingDataSet_* code, which expects ax4 and ax6 for markers.
% 
% %% Load ProcData
% S = load(procDataFileID,'ProcData');
% ProcData = S.ProcData;
% 
% %% Determine sampling rate and duration from available data (no trialDuration_sec)
% % Try to use a downsampled sampling rate if present
% Fs = [];
% if isfield(ProcData,'notes')
%     if isfield(ProcData.notes,'dsFs')
%         Fs = ProcData.notes.dsFs;
%     end
% end
% if isempty(Fs) || ~isnumeric(Fs) || ~isfinite(Fs)
%     % Fallback: assume 60 Hz if nothing is specified
%     Fs = 60;
% end
% 
% % Pick a reference signal to infer length
% refLen = [];
% % Prefer normalized ECoG if available, then raw
% if isfield(ProcData,'ECoG_DS_norm')
%     refLen = numel(ProcData.ECoG_DS_norm);
% elseif isfield(ProcData,'ECoG_DS')
%     refLen = numel(ProcData.ECoG_DS);
% else
%     error('Generate_SleepScoringFigure_Analog_SF:CannotDetermineDuration', ...
%         'Cannot determine trial duration from available ProcData fields.');
% end
% 
% trialDur = refLen / Fs;
% tVec     = (0:refLen-1)/Fs;
% 
% %% Pull signals with defensive field checks
% % --- Force sensor (prefer normalized) ---
% if isfield(ProcData,'forceSensor_norm')
%     force = ProcData.forceSensor_norm(:);
%     forceLabel = 'Force (norm)';
% else
%     if isfield(ProcData,'forceSensor')
%         force = ProcData.forceSensor(:);
%     else
%         force = nan(refLen,1);
%     end
%     forceLabel = 'Force (V)';
% end
% force = force(:);
% 
% % --- Binary force (movement) ---
% if isfield(ProcData,'binForceSensor')
%     binForce = ProcData.binForceSensor(:);
% else
%     binForce = false(refLen,1);
% end
% binForce = logical(binForce(:));
% 
% % --- EMG power (prefer normalized) ---
% emgPower = [];
% if isfield(ProcData,'EMG')
%     if isfield(ProcData.EMG,'emgPower_norm')
%         emgPower = ProcData.EMG.emgPower_norm(:);
%         emgLabel = 'EMG power (norm)';
%     elseif isfield(ProcData.EMG,'emgPower')
%         emgPower = ProcData.EMG.emgPower(:);
%         emgLabel = 'EMG power';
%     else
%         emgPower = nan(refLen,1);
%         emgLabel = 'EMG power';
%     end
% else
%     emgPower = nan(refLen,1);
%     emgLabel = 'EMG power';
% end
% 
% % --- ECoG downsampled (prefer normalized) ---
% ecogDS = [];
% if isfield(ProcData,'ECoG_DS_norm')
%     ecogDS = ProcData.ECoG_DS_norm(:);
%     ecogLabel = 'ECoG (DS norm)';
% elseif isfield(ProcData,'ECoG_DS')
%     ecogDS = ProcData.ECoG_DS(:);
%     ecogLabel = 'ECoG (DS)';
% else
%     ecogDS = [];
%     ecogLabel = 'ECoG (DS)';
% end
% 
% if ~isempty(ecogDS)
%     % If lengths mismatch slightly, trim/pad to refLen
%     if numel(ecogDS) > refLen
%         ecogDS = ecogDS(1:refLen);
%     elseif numel(ecogDS) < refLen
%         ecogDS(end+1:refLen) = ecogDS(end);
%     end
% end
% 
% %% Load spectrogram from a pre-computed file if possible
% [folder, base, ~] = fileparts(procDataFileID);
% % Drop the '_ProcData' suffix if present
% prefix = regexprep(base,'_ProcData$','');
% 
% specFile = '';
% % Try common patterns: <prefix>*Spec*.mat
% cand = dir(fullfile(folder, [prefix '*Spec*.mat']));
% if ~isempty(cand)
%     specFile = fullfile(folder, cand(1).name);
% end
% 
% Tspec = [];
% Fspec = [];
% Sspec = [];
% isNormSpec = false;
% 
% if ~isempty(specFile) && isfile(specFile)
%     try
%         L = load(specFile);
%         % Prefer SpecData.ECoG.normS if present
%         if isfield(L,'SpecData') && isfield(L.SpecData,'ECoG')
%             SD = L.SpecData;
%             if isfield(SD.ECoG,'normS')
%                 Sspec = SD.ECoG.normS;
%                 Fspec = SD.ECoG.F;
%                 Tspec = SD.ECoG.T;
%                 isNormSpec = true;
%             elseif all(isfield(SD.ECoG,{'S','F','T'}))
%                 Sspec = SD.ECoG.S;
%                 Fspec = SD.ECoG.F;
%                 Tspec = SD.ECoG.T;
%             end
%         end
% 
%         % Fallback heuristics if above failed
%         if isempty(Sspec)
%             fn = fieldnames(L);
%             for i = 1:numel(fn)
%                 val = L.(fn{i});
%                 if isstruct(val) && all(isfield(val,{'S','F','T'}))
%                     Sspec = val.S.*100;
%                     Fspec = val.F;
%                     Tspec = val.T;
%                     break
%                 end
%             end
%         end
%     catch ME
%         warning('Generate_SleepScoringFigure_Analog_SF:SpecLoadFailed', ...
%             'Failed to load spectrogram from %s (%s).', specFile, ME.message);
%     end
% end
% 
% %% Build figure
% figHandle = figure('Color','w');
% figHandle.Name = sprintf('Sleep scoring analog: %s', base);
% 
% % Try to expand the figure nicely
% try
%     set(figHandle,'Units','pixels');
%     mons = get(groot,'MonitorPositions');
%     fp   = get(figHandle,'OuterPosition');
%     ctr  = [fp(1)+fp(3)/2, fp(2)+fp(4)/2];
%     mIdx = 1;
%     for i = 1:size(mons,1)
%         m = mons(i,:);
%         if ctr(1) >= m(1) && ctr(1) <= m(1)+m(3) && ...
%            ctr(2) >= m(2) && ctr(2) <= m(2)+m(4)
%             mIdx = i; break
%         end
%     end
%     m = mons(mIdx,:);
%     newW = m(3);
%     newH = round(0.65 * m(4));
%     newX = m(1);
%     newY = m(2) + round(0.05 * m(4));
%     set(figHandle,'OuterPosition',[newX newY newW newH]);
% catch
% end
% 
% %% Subplot 1: Force only
% ax1 = subplot(5,1,1);
% hold(ax1,'on');
% 
% plot(ax1, tVec, force, 'LineWidth', 1);
% ylabel(ax1, forceLabel);
% xlim(ax1,[0 trialDur]);
% set(ax1,'TickLength',[0 0],'XTickLabel',[]);
% box(ax1,'off');
% title(ax1, strrep(base,'_',' '));
% 
% %% Subplot 2: EMG power (main scoring axis) + binForce overlay
% ax2 = subplot(5,1,2);
% hold(ax2,'on');
% plot(ax2, tVec, emgPower, 'LineWidth', 1,'Color','k');
% ylabel(ax2, emgLabel);
% xlim(ax2,[0 trialDur]);
% set(ax2,'TickLength',[0 0],'XTickLabel',[]);
% box(ax2,'off');
% 
% % ---- Move binForce overlay here (on EMG axis) ----
% if any(binForce)
%     goodEMG = emgPower(~isnan(emgPower));
%     if isempty(goodEMG)
%         ySpan = 1;
%         yMax  = 1;
%     else
%         ySpan = max(goodEMG) - min(goodEMG);
%         if ySpan <= 0 || ~isfinite(ySpan)
%             ySpan = 1;
%         end
%         yMax = max(goodEMG);
%     end
%     yBase  = yMax + 0.05*ySpan;
%     yMarks = yBase * ones(size(binForce));
%     tMarks = tVec(binForce);
%     yMarks = yMarks(binForce);
%     scatter(ax2, tMarks, yMarks, 10, 'filled');
% end
% 
% % %% Subplot 3: ECoG DS (if available)
% % ax3 = subplot(4,1,3);
% % hold(ax3,'on');
% % if ~isempty(ecogDS)
% %     plot(ax3, tVec, ecogDS, 'LineWidth', 1);
% %     ylabel(ax3, ecogLabel);
% % else
% %     ylabel(ax3, ecogLabel);
% %     text(ax3, trialDur*0.5, 0.5, 'ECoG DS not available', ...
% %         'HorizontalAlignment','center','VerticalAlignment','middle', ...
% %         'Color',[0.5 0.5 0.5],'FontAngle','italic');
% % end
% % xlim(ax3,[0 trialDur]);
% % set(ax3,'TickLength',[0 0],'XTickLabel',[]);
% % box(ax3,'off');
% %% Subplot 3: Raw EMG signal (replacing ECoG)
% ax3 = subplot(5,1,3);
% hold(ax3,'on');
% 
% % --- Pull raw EMG signal ---
% rawEMG = [];
% 
% % Preferred field: ProcData.EMG.emgSignal
% if isfield(ProcData,'EMG') && isfield(ProcData.EMG,'emgSignal_norm')
%     rawEMG = ProcData.EMG.emgSignal(:);
% end
% 
% % Final fallback
% if isempty(rawEMG)
%     rawEMG = nan(refLen,1);
% end
% 
% % Ensure correct length
% if numel(rawEMG) > refLen
%     rawEMG = rawEMG(1:refLen);
% elseif numel(rawEMG) < refLen
%     rawEMG(end+1:refLen) = rawEMG(end);
% end
% 
% % --- Plot ---
% plot(ax3, tVec, rawEMG, 'LineWidth',1);
% ylabel(ax3,'EMG (raw)');
% xlim(ax3,[0 trialDur]);
% set(ax3,'TickLength',[0 0],'XTickLabel',[]);
% box(ax3,'off');
% 
% %% Subplot 4: Spectrogram (log y, colorbar)
% axSpec = subplot(5,1,4:5);
% ax6    = axSpec;  % For compatibility with CreateTrainingDataSet_* code
% hold(axSpec,'on');
% 
% if ~isempty(Sspec) && ~isempty(Fspec) && ~isempty(Tspec)
% 
%     if isNormSpec
%         % Kevin-style normalization: normS = (S - baseline)./baseline
%         % Plot as percent change and use caxis([-100 100]) like Kevin
%         Splot = 100*Sspec;                          % % change from baseline
%         imagesc(axSpec, Tspec, Fspec, Splot);
%         axis(axSpec,'xy');
%         set(axSpec,'CLim',[-100 100]);
% 
%         cb = colorbar(axSpec);
%         ylabel(cb,'\Delta power (% of baseline)');
%     else
%         % Raw spectrogram: convert to dB, use robust CLim
%         Sdb = 10*log10(Sspec + eps);
%         imagesc(axSpec, Tspec, Fspec, Sdb);
%         axis(axSpec,'xy');
% 
%         p = prctile(Sdb(:),[5 95]);
%         if all(isfinite(p))
%             set(axSpec,'CLim',p);
%         end
% 
%         cb = colorbar(axSpec);
%         ylabel(cb,'Power (dB)');
%     end
% 
%     % ----- Accurate log-frequency ticks -----
%     set(axSpec,'YScale','log');
% 
%     % Candidate ticks; keep only those within Fspec range
%     candTicks = [1 2 4 8 15 30 60 80];
%     candTicks = candTicks(candTicks >= min(Fspec) & candTicks <= max(Fspec));
%     if isempty(candTicks)
%         candTicks = [min(Fspec) max(Fspec)];
%     end
%     set(axSpec,'YTick',candTicks, ...
%                'YTickLabel',arrayfun(@num2str,candTicks,'UniformOutput',false));
% 
%     ylabel(axSpec,'Freq (Hz)');
%     xlabel(axSpec,'Time (s)');
%     xlim(axSpec,[0 trialDur]);
% 
% else
%     ylabel(axSpec,'Freq (Hz)');
%     xlabel(axSpec,'Time (s)');
%     xlim(axSpec,[0 trialDur]);
%     text(axSpec, trialDur*0.5, 1, 'Spectrogram not available', ...
%         'HorizontalAlignment','center','VerticalAlignment','top', ...
%         'Color',[0.5 0.5 0.5],'FontAngle','italic');
% end
% set(axSpec,'TickLength',[0 0]);
% box(axSpec,'off');
% 
% 
% %% Link axes and make widths consistent
% linkaxes([ax1,ax2,ax3,axSpec],'x');
% xlim(ax1,[0 trialDur]);
% 
% pos1    = get(ax1,'Position');
% pos2    = get(ax2,'Position');
% pos3    = get(ax3,'Position');
% posSpec = get(axSpec,'Position');
% 
% % Match widths and left positions to ax1
% pos2(1)    = pos1(1); pos2(3)    = pos1(3);
% pos3(1)    = pos1(1); pos3(3)    = pos1(3);
% posSpec(1) = pos1(1); posSpec(3) = pos1(3);
% 
% set(ax2,'Position',pos2);
% set(ax3,'Position',pos3);
% set(axSpec,'Position',posSpec);
% %%x-ticks every 10 s, labels every 20 s on spectrogram axis ---
% tickStep = 10;  % tick spacing in seconds
% maxTick  = floor(trialDur/tickStep)*tickStep;
% xt = 0:tickStep:maxTick;
% 
% xtLbl = cell(size(xt));
% for i = 1:numel(xt)
%     if mod(xt(i),50) == 0       % label only 50-s multiples
%         xtLbl{i} = num2str(xt(i));
%     else
%         xtLbl{i} = '';
%     end
% end
% set(axSpec,'XTick',xt,'XTickLabel',xtLbl);
% 
% %% Outputs for scoring code
% ax4 = ax2;   % Use EMG axis as the "main" scoring axis
% % ax6 already set to spectrogram axis
% 
% %% Optional save
% if nargin > 1 && ischar(saveFigs) && strcmpi(saveFigs,'y')
%     [pDir,~,~] = fileparts(folder);
%     outDir = fullfile(pDir,'Figures','Sleep Scoring Figures');
%     if ~exist(outDir,'dir')
%         mkdir(outDir);
%     end
%     savefig(figHandle, fullfile(outDir, [base '_SleepScoring']));
%     saveas(figHandle, fullfile(outDir, [base '_SleepScoring']), 'tiff');
% end
% 
% end