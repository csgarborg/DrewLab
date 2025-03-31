function makeCalibrationComparePlot(fileStr1,fileStr2)
figure
subplot(3,1,1)
load('calibrationValues.mat');
plot(ones(length(calibrationValues.(['file_' fileStr1]).pixelDiff)),calibrationValues.(['file_' fileStr1]).pixelDiff,'ro')
hold on
plot(2*ones(length(calibrationValues.(['file_' fileStr2]).pixelDiff)),calibrationValues.(['file_' fileStr2]).pixelDiff,'bo')
axis([0 3 0 max([calibrationValues.(['file_' fileStr1]).pixelDiff calibrationValues.(['file_' fileStr2]).pixelDiff])+10])

subplot(3,1,2)
plot(1:length(calibrationValues.(['file_' fileStr1]).pixelDiff),calibrationValues.(['file_' fileStr1]).pixelDiff,'ro')
axis([0 length(calibrationValues.(['file_' fileStr1]).pixelDiff)+2 0 max([calibrationValues.(['file_' fileStr1]).pixelDiff calibrationValues.(['file_' fileStr2]).pixelDiff])+10])

subplot(3,1,3)
plot(1:length(calibrationValues.(['file_' fileStr2]).pixelDiff),calibrationValues.(['file_' fileStr2]).pixelDiff,'ro')
axis([0 length(calibrationValues.(['file_' fileStr2]).pixelDiff)+2 0 max([calibrationValues.(['file_' fileStr1]).pixelDiff calibrationValues.(['file_' fileStr2]).pixelDiff])+10])

figure(2)
