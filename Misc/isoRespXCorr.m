load('H:\24-06-11_MouseExp\240611_001_processed_Layer1_1.mat')
% uiopen('H:\24-06-11_MouseExp\Values.csv',1)
close all
valuesInterp=interp1(Values(:,1)*.03,Values(:,2),movementData.emgData(:,1),'linear');
valuesInterp(1) = 57.8;
subplot(2,1,1);plot(movementData.emgData(:,1),movementData.emgData(:,2));xlim([15 225]);subplot(2,1,2);plot((Values(:,1)*.03)-1.7,Values(:,2)*.01);xlim([15 225]);ylim([.5 1]);legend('EMG','behavioral camera');xlabel('seconds');title('EMG vs behavior camera - 1.7s');
figure(2);[c,lags] = xcorr(movementData.emgData(:,2)-mean(movementData.emgData(:,2)),valuesInterp*.02-mean(valuesInterp*.02),100,'coeff');stem(lags/30,c);ylim([0 1]);title('EMG vs behavior camera xcorr');xlabel('seconds')
figure(3);plot(movementData.emgData(:,2));hold on;plot(valuesInterp*0.02);legend('EMG','behavioral camera');xlabel('frames')