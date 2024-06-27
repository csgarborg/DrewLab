load('H:\24-06-11_MouseExp\240611_001_processed_Layer1_1.mat')
uiopen('H:\24-06-11_MouseExp\Values.csv',1)
close all
valuesInterp=interp1(Values(:,1)*.03,Values(:,2),movementData.emgData(:,1),'linear');
valuesInterp(1) = 57.8;
plot(movementData.emgData(:,1),movementData.emgData(:,2));hold on;plot((Values(:,1)*.03)-1.7,Values(:,2)*.02);xlim([100 115])
figure(2);[c,lags] = xcorr(movementData.emgData(:,2),valuesInterp*.02,100,'coeff');stem(lags,c)
figure(3);plot(movementData.emgData(:,2));hold on;plot(valuesInterp*0.02)