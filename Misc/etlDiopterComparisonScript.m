close all
load('calibrationValues.mat');
diopterValsPower = -1.27:.3:1.73;
powerVals = [48 48 49 49 49 49 49 49 49 50 50]; % 16x, 30% power, mW
diopterVals = [.23 .53 .83 1.13 1.43 1.73 -.07 -0.37 -0.67 -0.97 -1.27];
widthVals = [];
heightVals = [];
borderDistance = 220;
lowBoundPixel = borderDistance;
highBoundPixel = 512-borderDistance;
for n = 1:11
    if n < 10
        fieldName = ['file_210408_00' num2str(n)];
    else
        fieldName = ['file_210408_0' num2str(n)];
    end
    widthLocX = calibrationValues.(fieldName).pixelMidpointXWidth;
    widthLocY = calibrationValues.(fieldName).pixelMidpointYWidth;
    widthInd = lowBoundPixel <= widthLocX & widthLocX <= highBoundPixel & lowBoundPixel <= widthLocY & widthLocY <= highBoundPixel;
    meanWidth = calibrationValues.(fieldName).negativePolarityLengthInMicrons/median(calibrationValues.(fieldName).width(widthInd));
    widthVals(end+1) = meanWidth;
    
    heightLocX = calibrationValues.(fieldName).pixelMidpointXHeight;
    heightLocY = calibrationValues.(fieldName).pixelMidpointYHeight;
    heightInd = lowBoundPixel <= heightLocX & heightLocX <= highBoundPixel & lowBoundPixel <= heightLocY & heightLocY <= highBoundPixel;
    meanHeight = calibrationValues.(fieldName).negativePolarityLengthInMicrons/median(calibrationValues.(fieldName).height(heightInd));
    heightVals(end+1) = meanHeight;
end
[diopterVals,i] = sort(diopterVals);
widthVals = widthVals(i);
heightVals = heightVals(i);

figure(1)
subplot(3,1,1)
plot(diopterVals,widthVals,'mo')
xlabel('Diopter Input (1/m)')
ylabel('\mum/Pixel in X')
lsline
subplot(3,1,2)
plot(diopterVals,heightVals,'go')
xlabel('Diopter Input (1/m)')
ylabel('\mum/Pixel in Y')
lsline
subplot(3,1,3)
plot(diopterValsPower,powerVals,'ko')
xlabel('Diopter Input (1/m)')
ylabel('Measured Power from Objective')