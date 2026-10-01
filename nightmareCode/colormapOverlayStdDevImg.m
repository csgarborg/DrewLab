figure;
imshow(photo);
hold on;

h = imagesc(stdImg);

% Hot colormap
colormap(gca, hot);

% Threshold
stdThreshold = 5;   % <-- change this

% Transparency mask
alphaData = zeros(size(stdImg));
alphaData(stdImg > stdThreshold) = 0.6;

% Apply transparency
h.AlphaData = alphaData;

% Color range
clim([stdThreshold max(stdImg(:),[],'omitnan')]);

colorbar;
title('Pixel Standard Deviation');

hold off;