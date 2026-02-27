clear
close all
microsphereValues = csvread('H:\221207_010_microsphereValues.csv',1,0);
microsphereDarkValues = csvread('H:\221207_010_microsphereDarkValues.csv',1,0);
brainValues = csvread('H:\221207_010_brainValues.csv',1,0);
microsphereValues(:,1) = microsphereValues(:,1)*.0506;
brainValues(:,1) = brainValues(:,1)*.0506;
brainValues(end,:) = [];
microsphereDarkValues(:,1) = microsphereDarkValues(:,1)*.0506;

plot(brainValues(:,1),brainValues(:,2),'k')
hold on
plot(microsphereValues(:,1),microsphereValues(:,2),'k')
plot(microsphereDarkValues(:,1),microsphereDarkValues(:,2),'k')
xlim([148 398])

stationaryFind = microsphereValues(:,1) >= 148 & microsphereValues(:,1) <= 171;
stationaryMicrosphereMean = mean(microsphereValues(stationaryFind,2));
stationaryBrainMean = mean(brainValues(stationaryFind,2));
movementFind = microsphereValues(:,1) >= 180 & microsphereValues(:,1) <= 300;
movementMicrosphereMean = mean(microsphereValues(movementFind,2));
movementBrainMean = mean(brainValues(movementFind,2));

plot([148 398],[stationaryMicrosphereMean,stationaryMicrosphereMean],'m')
plot([148 398],[stationaryBrainMean,stationaryBrainMean],'g')
plot([148 398],[movementMicrosphereMean,movementMicrosphereMean],'m')
plot([148 398],[movementBrainMean,movementBrainMean],'g')

meanDiffMicrosphereMoveStationary = stationaryMicrosphereMean-movementMicrosphereMean;
meanDiffMicrosphereTotal = stationaryMicrosphereMean-mean(microsphereDarkValues(:,2));
meanDiffBrainMoveStationary = stationaryBrainMean-movementBrainMean;
meanDiffBrainTotal = stationaryBrainMean-mean(microsphereDarkValues(:,2));

percentChangeMicrospheres = meanDiffMicrosphereMoveStationary/meanDiffMicrosphereTotal;
percentChangeBrain = meanDiffBrainMoveStationary/meanDiffBrainTotal;

text(148,stationaryMicrosphereMean+250,[num2str(percentChangeMicrospheres*100) '% change in fluorescence - microspheres'])
text(148,stationaryBrainMean+250,[num2str(percentChangeBrain*100) '% change in fluorescence - microspheres'])

figure(2)
load('PSFData_FS.mat')
plot(PSFData.ZIntensity,PSFData.Z,'k')
hold on
scatter(PSFData.ZPointsIntensity,PSFData.ZPoints,'kx')
plot([1-percentChangeMicrospheres,1-percentChangeMicrospheres],[-6,6],'m')
plot([1-percentChangeBrain,1-percentChangeBrain],[-6,6],'g')
hold off
ylabel('\mum')
xlabel('Normalized Pixel Intensity')
title('Z')
xlim([0 1])
ylim([-7 7])
set(gca,'YDir','reverse')