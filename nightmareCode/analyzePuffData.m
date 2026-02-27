function analyzePuffData(dateList)

tdmsFilePaths = findTDMSFiles(dateList);

% for n = 1:length(tdmsFilePaths)
%     tdmsDataStruct = openTDMS(tdmsFilePaths{n});
%     puffComparisonPlot(tdmsDataStruct);
% end

puffMeanStruct = [];
for n = 1:length(tdmsFilePaths)
    if ~exist(strrep(tdmsFilePaths{n},'.tdms','_puffDataStruct.mat'))
        continue
    end
    load(strrep(tdmsFilePaths{n},'.tdms','_puffDataStruct.mat'))
    puffTimeSeconds = regexprep(puffDataStruct.puffDurationS, '(\.\d*?)0+$', '$1');
    puffTimeSeconds = ['Sec' regexprep(puffTimeSeconds, '\.$', '')];
    digitalFields = fieldnames(puffDataStruct.digitalData);
    analogFields = fieldnames(puffDataStruct.analogData);
    if ~isfield(puffMeanStruct,puffTimeSeconds)
        puffMeanStruct = buildEmptyPuffStruct(puffMeanStruct,puffTimeSeconds,digitalFields,analogFields);
        puffMeanStruct.(puffTimeSeconds).digitalTimescaleS = puffDataStruct.digitalTimescaleS;
        puffMeanStruct.(puffTimeSeconds).analogTimescaleS = puffDataStruct.analogTimescaleS;
    end
    for i = 1:length(digitalFields)
        puffMeanStruct.(puffTimeSeconds).digitalData.(digitalFields{i})(end+1,:) = puffDataStruct.digitalData.(digitalFields{i});
    end
    for i = 1:length(analogFields)
        puffMeanStruct.(puffTimeSeconds).analogData.(analogFields{i})(end+1,:) = puffDataStruct.analogData.(analogFields{i});
    end
    puffMeanStruct.(puffTimeSeconds).name{end+1} = [puffDataStruct.AnimalID ' - ' strrep(puffDataStruct.name, '_', ' ') ' - ' num2str(puffDataStruct.numPuffs) ' Puffs'];
end

close all
secondsFields = fieldnames(puffMeanStruct);
for s = 1:length(secondsFields)
    fields = fieldnames(puffMeanStruct.(secondsFields{s}).digitalData);
    fieldDFFAvgs = [];
    for i = 1:length(fields)
        figure(numel(findall(0, 'Type', 'figure')) + 1)
        X_dFF = normalizePixelDiff(puffMeanStruct.(secondsFields{s}).digitalData.(fields{i}));
        avgCalcMat = [];
        for n = 1:length(puffMeanStruct.(secondsFields{s}).name)
            subplot(3,1,1)
            plot(puffMeanStruct.(secondsFields{s}).digitalTimescaleS,puffMeanStruct.(secondsFields{s}).digitalData.(fields{i})(n,:))
            hold on
            subplot(3,1,2)
            plot(puffMeanStruct.(secondsFields{s}).digitalTimescaleS,X_dFF(n,:))
            avgCalcMat(n,:) = puffMeanStruct.(secondsFields{s}).digitalData.(fields{i})(n,:);
            hold on
        end
        subplot(3,1,1)
        ylimAxisVals = ylim;
        rectangle('Position', [0, 0, str2double(strrep(secondsFields{s}, 'Sec', '')), ylimAxisVals(2)], ...
            'EdgeColor', 'r', ...
            'LineWidth', 1, ...
            'FaceColor', [1 0 0], ...
            'FaceAlpha', 0.1);
        plot(puffMeanStruct.(secondsFields{s}).digitalTimescaleS,mean(avgCalcMat),'k')
        xlabel('Time (s)')
        ylabel(strrep(fields{i}, '_', ' '))
        xlim([puffMeanStruct.(secondsFields{s}).digitalTimescaleS(1) puffMeanStruct.(secondsFields{s}).digitalTimescaleS(end)])
        title([secondsFields(s) '-' strrep(fields{i}, '_', ' ')])
        legend(puffMeanStruct.(secondsFields{s}).name)
        
        subplot(3,1,2)
        ylimAxisVals = ylim;
        rectangle('Position', [0, -.5, str2double(strrep(secondsFields{s}, 'Sec', '')), ylimAxisVals(2)+.5], ...
            'EdgeColor', 'r', ...
            'LineWidth', 1, ...
            'FaceColor', [1 0 0], ...
            'FaceAlpha', 0.1);
        plot(puffMeanStruct.(secondsFields{s}).digitalTimescaleS,mean(X_dFF),'k')
        fieldDFFAvgs(i,:) = mean(X_dFF);
        xlabel('Time (s)')
        ylabel('\DeltaP/P')
        xlim([puffMeanStruct.(secondsFields{s}).digitalTimescaleS(1) puffMeanStruct.(secondsFields{s}).digitalTimescaleS(end)])
        title([secondsFields(s) '-' strrep(fields{i}, '_', ' ')])
        legend(puffMeanStruct.(secondsFields{s}).name)

        subplot(3,1,3)
        [avgSpec, avgSpec_dB, F, T] = averageSpectrogram(X_dFF', 1/(puffMeanStruct.(secondsFields{s}).digitalTimescaleS(2)-puffMeanStruct.(secondsFields{s}).digitalTimescaleS(1)));
        imagesc(T+puffMeanStruct.(secondsFields{s}).digitalTimescaleS(1), F, avgSpec_dB);
        axis xy;
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title('Average Spectrogram');
        colorbar;
        colormap turbo;
    end

    figure(numel(findall(0, 'Type', 'figure')) + 1)
    for n = 1:size(fieldDFFAvgs,1)
        plot(puffMeanStruct.(secondsFields{s}).digitalTimescaleS,fieldDFFAvgs(n,:))
        hold on
    end
    xlabel('Time (s)')
    ylabel('Mean \DeltaP/P')
    xlim([puffMeanStruct.(secondsFields{s}).digitalTimescaleS(1) puffMeanStruct.(secondsFields{s}).digitalTimescaleS(end)])
    rectangle('Position', [0, -.5, str2double(strrep(secondsFields{s}, 'Sec', '')), ylimAxisVals(2)+.5], ...
            'EdgeColor', 'r', ...
            'LineWidth', 1, ...
            'FaceColor', [1 0 0], ...
            'FaceAlpha', 0.1);
    title([secondsFields(s) ' - Mean Percent Pixel Diff Changes'])
    legend(cellfun(@(s) strtrim(regexprep(strrep(strrep(s,'mean',''),'_',' '),'\s+',' ')), ...
            fields, 'UniformOutput', false))

    fields = fieldnames(puffMeanStruct.(secondsFields{s}).analogData);
    for i = 1:length(fields)
        figure(numel(findall(0, 'Type', 'figure')) + 1)
        avgCalcMat = [];
        for n = 1:length(puffMeanStruct.(secondsFields{s}).name)
            plot(puffMeanStruct.(secondsFields{s}).analogTimescaleS,puffMeanStruct.(secondsFields{s}).analogData.(fields{i})(n,:))
            avgCalcMat(n,:) = puffMeanStruct.(secondsFields{s}).analogData.(fields{i})(n,:);
            hold on
        end
        ylimAxisVals = ylim;
        rectangle('Position', [0, 0, str2double(strrep(secondsFields{s}, 'Sec', '')), ylimAxisVals(2)], ...
            'EdgeColor', 'r', ...
            'LineWidth', 1, ...
            'FaceColor', [1 0 0], ...
            'FaceAlpha', 0.1);
        plot(puffMeanStruct.(secondsFields{s}).analogTimescaleS,mean(avgCalcMat),'k')
        xlabel('Time (s)')
        ylabel(strrep(fields{i}, '_', ' '))
        xlim([puffMeanStruct.(secondsFields{s}).analogTimescaleS(1) puffMeanStruct.(secondsFields{s}).analogTimescaleS(end)])
        title([secondsFields(s) '-' strrep(fields{i}, '_', ' ')])
        legend(puffMeanStruct.(secondsFields{s}).name)
    end

    % saveFileString = [tdmsDataStruct.filePath(1:end-6) strrep(tdmsDataStruct.name, '_converted', '')];
    % 
    % figs = findall(0, 'Type', 'figure');
    % savefig(figs,[saveFileString '_puffMeanFigs.fig'])
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function puffMeanStruct = buildEmptyPuffStruct(puffMeanStruct,puffTimeSeconds,digitalFields,analogFields)
for i = 1:length(digitalFields)
    puffMeanStruct.(puffTimeSeconds).digitalData.(digitalFields{i}) = [];
end
for i = 1:length(analogFields)
    puffMeanStruct.(puffTimeSeconds).analogData.(analogFields{i}) = [];
end
puffMeanStruct.(puffTimeSeconds).name = {};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X_dFF = normalizePixelDiff(X)
X = X';
baselineIdx = 1:450;
F0 = mean(X(baselineIdx, :), 1);
F0(F0 == 0) = eps;
X_dFF = (X - F0) ./ F0;
X_dFF = X_dFF';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [avgSpec, avgSpec_dB, F, T] = averageSpectrogram(clips, Fs, varargin)
% averageSpectrogram
% Computes the average spectrogram across equal-length clips
%
% USAGE:
%   [avgSpec, avgSpec_dB, F, T] = averageSpectrogram(clips, Fs)
%
% INPUTS:
%   clips : cell array {N} of vectors OR matrix [samples x N]
%   Fs    : sampling rate (Hz)
%
% OPTIONAL NAME-VALUE PAIRS:
%   'WindowLength' : window length in samples (default: 256)
%   'Overlap'      : overlap in samples (default: 200)
%   'NFFT'         : FFT length (default: 512)
%   'Window'       : window function (default: hann)
%
% OUTPUTS:
%   avgSpec     : average power spectrogram
%   avgSpec_dB  : average spectrogram in dB
%   F           : frequency axis (Hz)
%   T           : time axis (s)

% ----------------------------
% Parse inputs
% ----------------------------
p = inputParser;
addRequired(p, 'clips');
addRequired(p, 'Fs');
addParameter(p, 'WindowLength', 256);
addParameter(p, 'Overlap', 200);
addParameter(p, 'NFFT', 512);
addParameter(p, 'Window', []);
parse(p, clips, Fs, varargin{:});

winLength = p.Results.WindowLength;
noverlap  = p.Results.Overlap;
nfft      = p.Results.NFFT;

if isempty(p.Results.Window)
    window = hann(winLength);
else
    window = p.Results.Window;
end

% ----------------------------
% Convert clips to matrix if needed
% ----------------------------
if iscell(clips)
    numClips = numel(clips);
    x = clips{1};
    numSamples = numel(x);
    X = zeros(numSamples, numClips);

    for k = 1:numClips
        X(:,k) = clips{k}(:);
    end
else
    X = clips;
    numClips = size(X,2);
end

% ----------------------------
% First spectrogram (for sizing)
% ----------------------------
[S, F, T] = spectrogram(X(:,1), window, noverlap, nfft, Fs);
specSum = abs(S).^2;

% ----------------------------
% Loop over clips
% ----------------------------
for k = 2:numClips
    S = spectrogram(X(:,k), window, noverlap, nfft, Fs);
    specSum = specSum + abs(S).^2;
end

% ----------------------------
% Average and convert to dB
% ----------------------------
avgSpec = specSum / numClips;
avgSpec_dB = 10*log10(avgSpec + eps);

end