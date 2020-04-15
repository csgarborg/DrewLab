close all
load('calibrationValues.mat')
ax=calibrationValues.file_191210_015.surfaceCalibFitX;
ay=calibrationValues.file_191210_015.surfaceCalibFitY;
bx=calibrationValues.file_191210_016.surfaceCalibFitX;
by=calibrationValues.file_191210_016.surfaceCalibFitY;
[x,y] = meshgrid(100:10:400,100:10:400);
axz = ax(x,y);
ayz = ay(x,y);
bxz = bx(x,y);
byz = by(x,y);
xDiff = abs(axz-bxz);
yDiff = abs(ayz-byz);
xDiffMean = mean(mean(xDiff));
yDiffMean = mean(mean(yDiff));
xMean = mean(mean(axz));
yMean = mean(mean(ayz));
xPercentChange = round((xDiffMean/xMean)*100,2);
yPercentChange = round((yDiffMean/yMean)*100,2);

figure(1)
pcolor(x,y,xDiff)
shading interp
axis([0 512 0 512])
colorbar
title('Absolute Difference in X Calibration Between ETL Planes (-2.1 and 1 Diopters)')
xlabel('Pixels (X)')
ylabel('Pixels (Y)')
h=colorbar;
ylabel(h,'|\Delta| (\mum/pixel)','FontSize',12)
text(50,450,['Average percent change in \mum/pixel between planes: ' num2str(xPercentChange) '%']);

figure(2)
pcolor(x,y,yDiff)
shading interp
axis([0 512 0 512])
colorbar
title('Absolute Difference in Y Calibration Between ETL Planes (-2.1 and 1 Diopters)')
xlabel('Pixels (X)')
ylabel('Pixels (Y)')
h=colorbar;
ylabel(h,'|\Delta| (\mum/pixel)','FontSize',12)
text(50,450,['Average percent change in \mum/pixel between planes: ' num2str(yPercentChange) '%']);