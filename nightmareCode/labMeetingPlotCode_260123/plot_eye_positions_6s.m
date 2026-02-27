function plot_eye_positions_6s(csvFile)

    fps = 60;

    % Segment definition
    segLenFrames = 16*fps + 5;   % 16 s + 5 frames = 965 frames
    segLenSec    = segLenFrames / fps;

    T = readtable(csvFile, 'PreserveVariableNames', true);

    headers = {
        'eye_top.x'
        'eye_top.y'
        'eye_bottom.x'
        'eye_bottom.y'
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
        ax = subplot(4,1,i);
        hold on

        y = T.(headers{i});

        for s = 1:nSeg
            idx = (s-1)*segLenFrames + (1:segLenFrames);
            plot(tSeg, y(idx), 'Color', colors(s,:), 'LineWidth', 1);
        end

        xlim([0 segLenSec])
        yl = ylim;

        % Red box: 2–8 seconds
        rectangle('Position', [2, yl(1), 6, yl(2)-yl(1)], ...
                  'FaceColor', [1 0.4 0.4 0.15], ...
                  'EdgeColor', 'none');

        uistack(findobj(ax,'Type','line'),'top')

        hold off
        ylabel(headers{i}, 'Interpreter', 'none')
        grid on

        if i == numel(headers)
            xlabel('Time within segment (s)')
        else
            set(gca,'XTickLabel',[])
        end
    end

    sgtitle('Eye Positions (16 s + 5 frames, segmented)')
end
