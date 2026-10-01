function plotMovieStdDevResults(resultsPath)
    load(resultsPath);
    figure('Name', 'Mouse Motion StdDev Analysis');
    % imagesc(results.stdMap);
    imagesc((results.stdMap).^0.9);
    axis image off; colormap(gca, hot); colorbar;
    title('Std-Dev Projection (whole video)');

    figure('Name', 'Mouse Motion Energy Analysis');
    % imagesc(results.motionEnergyMap);
    imagesc((results.motionEnergyMap).^0.55);
    axis image off; colormap(gca, hot); colorbar;
    title('Motion Energy Map (std of frame diffs)');

    figure('Name', 'Frame Difference Analysis');
    t = (1:length(results.motionTrace)) / results.fps;
    % plot(t, results.motionTrace, 'w', 'LineWidth', 1);
    plot(t, medfilt1(results.motionTrace,3), 'w', 'LineWidth', 1);
    xlabel('Time (s)'); ylabel('Mean |\Deltapixel|');
    title('Frame-Differencing Motion Trace');
    xlim([0 t(end)]);
    grid on;

    if isfield(results, 'windowOutputDir')
        fprintf('Windowed std-dev maps saved as individual .mat files in:\n  %s\n', ...
            results.windowOutputDir);
        fprintf('Load one with: load(fullfile(results.windowOutputDir, ''window_0001_...mat''))\n');

        % Show a quick preview montage of a subsample of windows so you
        % don't have to load hundreds of files by hand
        d = dir(fullfile(results.windowOutputDir, 'window_*.mat'));
        if ~isempty(d)
            nShow = min(16, numel(d));
            showIdx = round(linspace(1, numel(d), nShow));
            figure('Name', 'Windowed Std-Dev Preview', 'Color', 'w');
            for i = 1:nShow
                s = load(fullfile(d(showIdx(i)).folder, d(showIdx(i)).name));
                subplot(ceil(sqrt(nShow)), ceil(sqrt(nShow)), i);
                imagesc(s.windowStd); axis image off; colormap(hot);
                title(sprintf('win %d', showIdx(i)), 'FontSize', 8);
            end
            sgtitle('Sample of windowed std-dev projections over time');
        end
    end
end