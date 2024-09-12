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
W = x.imageWidthPixels;
H = x.imageHeightPixels;

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

% plot 3d surface plot of x calibration values across image
h(3) = figure('Color','White');
% micronLengthValWidth = ones(size(widthData,1),1);
% micronLengthValWidth(widthData(:,6) == 1,1) = barSize;
% micronLengthValWidth(widthData(:,6) == -1,1) = holeSize;
% surfaceCalibFitX = fit(widthData(:,1:2),micronLengthValWidth./widthData(:,3),'poly55');
% plot(x.surfaceCalibFitX,widthData(:,1:2),micronLengthValWidth./widthData(:,3));
plot(x.surfaceCalibFitX);
axis([1 W 1 H 0.5 ceil(max(micronLengthValWidth./widthData(:,3)))])
title(['\fontsize{20pt}\bf{Surface Calibration Fit (X)}' 10 '\fontsize{10pt}\rm{210408_001}' 10 '\fontsize{10pt}\rm{No Comments}'])
xlabel('Pixels (X)')
ylabel('Pixels (Y)')
zlabel('\mum/Pixel')

% plot 3d surface plot of y calibration values across image
h(4) = figure('Color','White');
% micronLengthValHeight = ones(size(heightData,1),1);
% micronLengthValHeight(heightData(:,6) == 1,1) = barSize;
% micronLengthValHeight(heightData(:,6) == -1,1) = holeSize;
% surfaceCalibFitY = fit(heightData(:,1:2),micronLengthValHeight./heightData(:,3),'poly55');
% plot(x.surfaceCalibFitY,heightData(:,1:2),micronLengthValHeight./heightData(:,3));
plot(x.surfaceCalibFitY);
axis([1 W 1 H 0.5 ceil(max(micronLengthValHeight./heightData(:,3)))])
title(['\fontsize{20pt}\bf{Surface Calibration Fit (Y)}' 10 '\fontsize{10pt}\rm{210408_001}' 10 '\fontsize{10pt}\rm{No Comments}'])
xlabel('Pixels (X)')
ylabel('Pixels (Y)')
zlabel('\mum/Pixel')

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