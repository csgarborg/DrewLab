close all
load('calibrationValues.mat')
x = calibrationValues.file_210408_001;
widthData = [x.pixelMidpointXWidth x.pixelMidpointYWidth x.width x.leftBoundWidth x.rightBoundWidth x.polarityWidth];
heightData = [x.pixelMidpointXHeight x.pixelMidpointYHeight x.height x.topBoundHeight x.bottomBoundHeight x.polarityHeight];
barSize = 6;
holeSize = 19;
micronLengthValWidth = ones(size(widthData,1),1);
micronLengthValWidth(widthData(:,6) == 1,1) = barSize;
micronLengthValWidth(widthData(:,6) == -1,1) = holeSize;
micronLengthValHeight = ones(size(heightData,1),1);
micronLengthValHeight(heightData(:,6) == 1,1) = barSize;
micronLengthValHeight(heightData(:,6) == -1,1) = holeSize;

h(1) = figure('Color','White');
scatter(widthData(:,1),widthData(:,2),50,micronLengthValWidth./widthData(:,3),'filled');
colorbar;
k = colorbar;
ylabel(k,'\mum/pixel','FontSize',15)
hold off
axis equal
axis([1 512 1 512])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title(['\fontsize{20pt}\bf{Position of X Calibration Measurements}'])
xlabel('Pixels')
ylabel('Pixels')

h(2) = figure('Color','White');
scatter(heightData(:,1),heightData(:,2),50,micronLengthValHeight./heightData(:,3),'filled');
colorbar;
k = colorbar;
ylabel(k,'\mum/pixel','FontSize',15)
hold off
axis equal
axis([1 512 1 512])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title(['\fontsize{20pt}\bf{Position of Y Calibration Measurements}'])
xlabel('Pixels')
ylabel('Pixels')

% for n = 1:11
%     if n < 10
%         fieldName = ['file_210408_00' num2str(n)];
%     else
%         fieldName = ['file_210408_0' num2str(n)];
%     end
%     widthData = [calibrationValues.(fieldName).pixelMidpointXWidth calibrationValues.(fieldName).pixelMidpointYWidth calibrationValues.(fieldName).width calibrationValues.(fieldName).leftBoundWidth calibrationValues.(fieldName).rightBoundWidth calibrationValues.(fieldName).polarityWidth];
%     heightData = [calibrationValues.(fieldName).pixelMidpointXHeight calibrationValues.(fieldName).pixelMidpointYHeight calibrationValues.(fieldName).height calibrationValues.(fieldName).topBoundHeight calibrationValues.(fieldName).bottomBoundHeight calibrationValues.(fieldName).polarityHeight];
%     barSize = 6;
%     holeSize = 19;
%     micronLengthValWidth = ones(size(widthData,1),1);
%     micronLengthValWidth(widthData(:,6) == 1,1) = barSize;
%     micronLengthValWidth(widthData(:,6) == -1,1) = holeSize;
%     micronLengthValHeight = ones(size(heightData,1),1);
%     micronLengthValHeight(heightData(:,6) == 1,1) = barSize;
%     micronLengthValHeight(heightData(:,6) == -1,1) = holeSize;