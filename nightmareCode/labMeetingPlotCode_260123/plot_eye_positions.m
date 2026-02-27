function plot_eye_positions(csvFile)
% PLOT_EYE_POSITIONS
% Reads a CSV file and plots eye feature positions segmented into
% 12 seconds + 5 frames long chunks, overlaid in each subplot.
%
% Also draws a red box from 2–4 seconds in segment time.
%
% Expected headers:
%   eye_top.x
%   eye_top.y
%   eye_bottom.x
%   eye_bottom.y
%
% Sampling rate: 60 FPS

    fps = 60;

    % Segment definition
    segLenFrames = 12*fps + 5;      % 12 s + 5 frames = 725 frames
    segLenSec    = segLenFrames / fps;

    % Read table
    T = readtable(csvFile, 'PreserveVariableNames', true);

    headers = { ...
        'eye_top.x', ...
        'eye_top.y', ...
        'eye_bottom.x', ...
        'eye_bottom.y'};

    % Check columns exist
    for i = 1:numel(headers)
        if ~ismember(headers{i}, T.Properties.VariableNames)
            error('Missing column: %s', headers{i});
        end
    end

    nFrames = height(T);
    nSeg = floor(nFrames / segLenFrames);

    % Time axis within each segment
    tSeg = (0:segLenFrames-1) / fps;

    % Colors for segments
    colors = lines(nSeg);

    figure;

    for i = 1:numel(headers)
        ax = subplot(4,1,i);
        hold on

        y = T.(headers{i});

        % Plot each segment with a different color
        for s = 1:nSeg
            idxStart = (s-1)*segLenFrames + 1;
            idxEnd   = idxStart + segLenFrames - 1;

            ySeg = y(idxStart:idxEnd);
            plot(tSeg, ySeg, ...
                 'Color', colors(s,:), ...
                 'LineWidth', 1);
        end

        % Fix x-axis range
        xlim([0 segLenSec])

        % Y limits for red box
        yl = ylim;

        % Draw red box from 2 to 4 seconds
        rectangle('Position', [2, yl(1), 2, yl(2)-yl(1)], ...
                  'FaceColor', [1 0.4 0.4 0.15], ...
                  'EdgeColor', 'none');

        % Keep data on top
        uistack(findobj(ax,'Type','line'),'top')

        hold off
        xlabel('Time within segment (s)')
        ylabel(headers{i}, 'Interpreter', 'none')
        grid on
    end

    sgtitle('Eye Feature Positions (Segmented, Color-Coded)')
end

