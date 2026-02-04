function plot_nose_positions_6s(csvFile)

    fps = 60;

    segLenFrames = 16*fps + 5;
    segLenSec    = segLenFrames / fps;

    T = readtable(csvFile, 'PreserveVariableNames', true);

    headers = {
        'nose.x'
        'nose.y'
    };

    for i = 1:numel(headers)
        if ~ismember(headers{i}, T.Properties.VariableNames)
            error('Missing column: %s', headers{i});
        end
    end

    nFrames = height(T);
    nSeg = floor(nFrames / segLenFrames);
    tSeg = (0:segLenFrames-1) / fps;

    colors = lines(nSeg);

    figure;

    for i = 1:numel(headers)
        ax = subplot(2,1,i);
        hold on

        y = T.(headers{i});

        for s = 1:nSeg
            idx = (s-1)*segLenFrames + (1:segLenFrames);
            plot(tSeg, y(idx), 'Color', colors(s,:), 'LineWidth', 1);
        end

        xlim([0 segLenSec])
        yl = ylim;

        rectangle('Position', [2, yl(1), 6, yl(2)-yl(1)], ...
                  'FaceColor', [1 0.4 0.4 0.15], ...
                  'EdgeColor', 'none');

        uistack(findobj(ax,'Type','line'),'top')

        hold off
        ylabel(headers{i}, 'Interpreter', 'none')
        xlabel('Time within segment (s)')
        grid on
    end

    sgtitle('Nose Position (16 s + 5 frames, segmented)')
end
