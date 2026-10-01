function pairwiseBarPlotCrossCorrelation
% Example data - replace with your actual vectors
group2 = [.3558, .3358, .3919];   % vector 1 (3 values) (REM mean norm xcorr)
group1 = [.4854, .5458, .5267];   % vector 2 (3 values, paired with group1) (awake mean norm xcorr)

figure; hold on;

awakeColor = [0 0.4470 0.7410];
remColor   = [0.8500 0.3250 0.0980];

% Define custom RGB colors for each bar
color2 = remColor;
color1 = awakeColor;

% Bar plot of means
means = [mean(group1), mean(group2)];
b = bar([1 2], means, 0.5, 'FaceColor', 'flat', 'EdgeColor', 'k');
b.CData(1,:) = color1;
b.CData(2,:) = color2;

% Overlay individual paired points connected by lines
n = length(group1);
for i = 1:n
    plot([1 2], [group1(i) group2(i)], '-o', ...
        'Color', [0.3 0.3 0.3], ...
        'MarkerFaceColor', [0.3 0.3 0.3], ...
        'MarkerSize', 5, ...
        'LineWidth', 1);
end

% Formatting
xlim([0.5 2.5]);
set(gca, 'XTick', [1 2], 'XTickLabel', {'Awake', 'REM'});
ylabel('Mean Normalized Cross-Correlation');
title('Paired Comparison of Cross Correlations Between Mice (n=3)');
box off;
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Example data - replace with your actual vectors
group2 = [44.467, 38.0907, 40.1195];   % vector 1 (3 values) (REM PC1)
group1 = [56.3279, 57.7395, 52.0629];   % vector 2 (3 values, paired with group1) (awake PC1)

figure; hold on;

awakeColor = [0 0.4470 0.7410];
remColor   = [0.8500 0.3250 0.0980];

% Define custom RGB colors for each bar
color2 = remColor;
color1 = awakeColor;

% Bar plot of means
means = [mean(group1), mean(group2)];
b = bar([1 2], means, 0.5, 'FaceColor', 'flat', 'EdgeColor', 'k');
b.CData(1,:) = color1;
b.CData(2,:) = color2;

% Overlay individual paired points connected by lines
n = length(group1);
for i = 1:n
    plot([1 2], [group1(i) group2(i)], '-o', ...
        'Color', [0.3 0.3 0.3], ...
        'MarkerFaceColor', [0.3 0.3 0.3], ...
        'MarkerSize', 5, ...
        'LineWidth', 1);
end

% Formatting
xlim([0.5 2.5]);
ylim([0 60]);
set(gca, 'XTick', [1 2], 'XTickLabel', {'Awake', 'REM'});
ylabel('% Variance');
title('Paired Comparison of PC1 Strength Between Mice (n=3)');
box off;
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Example data - replace with your actual vectors
group2 = [12.3017, 10.5318, 11.2082];   % vector 1 (3 values) (REM PC2)
group1 = [8.08689, 7.60844, 11.2972];   % vector 2 (3 values, paired with group1) (Awake PC2)

figure; hold on;

awakeColor = [0 0.4470 0.7410];
remColor   = [0.8500 0.3250 0.0980];

% Define custom RGB colors for each bar
color2 = remColor;
color1 = awakeColor;

% Bar plot of means
means = [mean(group1), mean(group2)];
b = bar([1 2], means, 0.5, 'FaceColor', 'flat', 'EdgeColor', 'k');
b.CData(1,:) = color1;
b.CData(2,:) = color2;

% Overlay individual paired points connected by lines
n = length(group1);
for i = 1:n
    plot([1 2], [group1(i) group2(i)], '-o', ...
        'Color', [0.3 0.3 0.3], ...
        'MarkerFaceColor', [0.3 0.3 0.3], ...
        'MarkerSize', 5, ...
        'LineWidth', 1);
end

% Formatting
xlim([0.5 2.5]);
ylim([0 60]);
set(gca, 'XTick', [1 2], 'XTickLabel', {'Awake', 'REM'});
ylabel('% Variance');
title('Paired Comparison of PC2 Strength Between Mice (n=3)');
box off;
hold off;