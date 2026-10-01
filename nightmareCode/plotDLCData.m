function dlcData = plotDLCData(csvFile, movieFile, likelihoodThresh)
% PLOT_DLC_TRACES Load a DeepLabCut CSV + source movie and plot traces.
%
%   dlcData = plot_dlc_traces(csvFile, movieFile) loads DLC tracking data
%   from csvFile (standard DLC CSV export: scorer/bodyparts/coords header
%   rows) and the corresponding movieFile, then produces:
%
%     1. One figure per keypoint, with x, y, and likelihood vs. frame as
%        stacked subplots.
%     2. One figure with the 2D trajectories of all keypoints, axis-
%        matched to the movie's pixel width/height: trajectories alone
%        on the left, trajectories overlaid semi-transparently on the
%        first movie frame on the right.
%
%   dlcData = plot_dlc_traces(csvFile, movieFile, likelihoodThresh) lets
%   you set the confidence threshold used to mask low-confidence points
%   as NaN (default 0.9).
%
%   Figures are left on whatever renderer/colors your MATLAB root
%   defaults specify (e.g. dark theme) -- nothing here forces a white
%   background.

if nargin < 3 || isempty(likelihoodThresh)
    likelihoodThresh = 0.6;
end

dlcData = load_dlc_data(csvFile);

bodyparts = fieldnames(dlcData);
bodyparts(strcmp(bodyparts, 'frame')) = [];
nParts = numel(bodyparts);
frame  = dlcData.frame;

% Load the movie once up front so we have frame rate + dimensions + the
% first frame available for both the time axis and the trajectory plots.
vidObj     = VideoReader(movieFile);
firstFrame = read(vidObj, 1);
vidWidth   = vidObj.Width;
vidHeight  = vidObj.Height;
frameRate  = vidObj.FrameRate;

timeSec = frame / frameRate;

%% ---- Figure(s) 1: per-keypoint x / y / likelihood subplots ----
for i = 1:nParts
    part = bodyparts{i};
    x   = dlcData.(part).x;
    y   = dlcData.(part).y;
    lik = dlcData.(part).likelihood;

    xMasked = x;
    yMasked = y;
    xMasked(lik < likelihoodThresh) = NaN;
    yMasked(lik < likelihoodThresh) = NaN;

    figure('Name', sprintf('DLC trace: %s', part));

    subplot(3, 1, 1);
    plot(timeSec, xMasked, 'LineWidth', 1);
    title(sprintf('%s - x position', strrep(part, '_', ' ')), 'Interpreter', 'none');
    xlabel('Time (s)'); ylabel('X (px)');
    grid on;

    subplot(3, 1, 2);
    plot(timeSec, yMasked, 'LineWidth', 1);
    title(sprintf('%s - y position', strrep(part, '_', ' ')), 'Interpreter', 'none');
    xlabel('Time (s)'); ylabel('Y (px)');
    grid on;

    subplot(3, 1, 3);
    plot(timeSec, lik, 'LineWidth', 1);
    yline(likelihoodThresh, '--r');
    ylim([0 1]);
    title('Likelihood', 'Interpreter', 'none');
    xlabel('Time (s)'); ylabel('Likelihood');
    grid on;
end

%% ---- Figure 2: trajectories alone + trajectories over first frame ----
colors = lines(nParts);

figure('Name', 'DLC trajectories');

% --- left: trajectories only, axis matched to movie dimensions ---
subplot(1, 2, 1);
hold on;
for i = 1:nParts
    part = bodyparts{i};
    x   = dlcData.(part).x;
    y   = dlcData.(part).y;
    lik = dlcData.(part).likelihood;
    x(lik < likelihoodThresh) = NaN;
    y(lik < likelihoodThresh) = NaN;

    scatter(x, y, 8, colors(i, :), 'filled', ...
        'DisplayName', strrep(part, '_', ' '));
end
set(gca, 'YDir', 'reverse');
axis equal;
xlim([0 vidWidth]);
ylim([0 vidHeight]);
xlabel('X (px)'); ylabel('Y (px)');
title('Trajectories');
legend('show', 'Interpreter', 'none', 'Location', 'bestoutside');
grid on;

% --- right: trajectories overlaid semi-transparently on first frame ---
subplot(1, 2, 2);
imshow(firstFrame);
hold on;
for i = 1:nParts
    part = bodyparts{i};
    x   = dlcData.(part).x;
    y   = dlcData.(part).y;
    lik = dlcData.(part).likelihood;
    x(lik < likelihoodThresh) = NaN;
    y(lik < likelihoodThresh) = NaN;

    scatter(x, y, 8, colors(i, :), 'filled', ...
        'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0.35, ...
        'DisplayName', strrep(part, '_', ' '));
end
title('Trajectories over first frame');

end


%% ===================== LOCAL FUNCTION =====================
function dlcData = load_dlc_data(csvFile)
% LOAD_DLC_DATA Load a DeepLabCut CSV into a struct of per-bodypart traces.
%
%   dlcData = load_dlc_data(csvFile) returns a struct where:
%       dlcData.frame                 - frame index (column vector)
%       dlcData.<bodypart>.x          - x position over frames
%       dlcData.<bodypart>.y          - y position over frames
%       dlcData.<bodypart>.likelihood - confidence over frames

fid = fopen(csvFile, 'r');
if fid == -1
    error('Could not open file: %s', csvFile);
end
fgetl(fid);                   % scorer line (unused)
bodypartLine = fgetl(fid);    % bodyparts line
fclose(fid);

bodyparts = strsplit(bodypartLine, ',');
bodyparts = bodyparts(2:end);   % drop the leading empty/frame-index cell

% Read numeric data, skipping the 3 header rows
numData  = readmatrix(csvFile, 'NumHeaderLines', 3);
frameIdx = numData(:, 1);
numData  = numData(:, 2:end);

dlcData = struct();
dlcData.frame = frameIdx;

uniqueParts = unique(bodyparts, 'stable');
for i = 1:numel(uniqueParts)
    part = uniqueParts{i};
    cols = find(strcmp(bodyparts, part));   % should be [x, y, likelihood] columns
    fieldName = matlab.lang.makeValidName(part);
    dlcData.(fieldName).x          = numData(:, cols(1));
    dlcData.(fieldName).y          = numData(:, cols(2));
    dlcData.(fieldName).likelihood = numData(:, cols(3));
end

end