%load('211216_001_processe_2layerBrainInSkullDataFinal.mat');
load('D:\21-12-16_MouseExp\211216_002_processe_2layerBrainInSkullDataFinal.mat');

figure(444);
subplot(311)
movement_time=movementData.secondsPerFrame*(1:length(movementData.targetPosition));
stackedplot(movement_time, movementData.targetPosition)
%%
subplot(312)
plot(movementData.emgData(:,1),movementData.emgData(:,2))
xlim([min(movementData.emgData(:,1)) max(movementData.emgData(:,1))])
set(gca, 'YScale', 'log')

movementData.emgDataInterp=zeros(size(movementData.targetPosition));
movementData.emgDataInterp(:,1)=movement_time;
movementData.emgDataInterp(:,2)=interp1(movementData.emgData(:,1),movementData.emgData(:,2),movement_time,'linear');

% subplot(313)
% % ballData = load('H:\21-12-16_MouseExp\211216_002.txt');
% emgDataOnly = [ballData(:,1) ballData(:,3)];
% emgDataOnly = downsample(emgDataOnly,200);
% plot(emgDataOnly(:,1),emgDataOnly(:,2))
% xlim([min(movementData.emgData(:,1)) max(movementData.emgData(:,1))])

%%
figure(433)
maxlag=500;
xc_1=xcorr(detrend(movementData.emgDataInterp(1:(end-100),2))', detrend(movementData.targetPosition(1:(end-100),1))',maxlag,'coeff');
xc_2=xcorr(detrend(movementData.emgDataInterp(1:(end-100),2))', detrend(movementData.targetPosition(1:(end-100),2))',maxlag,'coeff');
hold off
plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_1)
hold on
plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_2)

best_lag=round(1/movementData.secondsPerFrame)%about 1 second

%%
figure(432)
emg_bins=.5:.1:3.5;
pixel_bins=-3:.5:7;
subplot(221)
scatter(movementData.emgDataInterp(1:(end-best_lag),2),movementData.targetPosition(best_lag:end-1,1))
xlabel('emg')
ylabel('shift, pixels')
subplot(222)
scatter(movementData.emgDataInterp(1:(end-best_lag),2),movementData.targetPosition(best_lag:end-1,2))
xlabel('emg')
ylabel('shift, pixels')
subplot(223)
hh1=histogram2(movementData.emgDataInterp(1:(end-best_lag),2),movementData.targetPosition(best_lag:end-1,1), emg_bins, pixel_bins,...
    'DisplayStyle','tile','ShowEmptyBins','on');
imagesc(emg_bins, pixel_bins,log(hh1.Values'+1))
axis xy
xlabel('emg')
ylabel('shift, pixels')
subplot(224)
hh2=histogram2(movementData.emgDataInterp(1:(end-best_lag),2),movementData.targetPosition(best_lag:end-1,2), emg_bins, pixel_bins,...
    'DisplayStyle','tile','ShowEmptyBins','on');
imagesc(emg_bins, pixel_bins,log(hh2.Values'+1))
axis xy
colormap hot
xlabel('emg')
ylabel('shift, pixels')
%%

movementData.locDataInterp=zeros(size(movementData.targetPosition));
movementData.locDataInterp(:,1)=movement_time;
movementData.locDataInterp(:,2)=interp1(movementData.ballData(:,1),movementData.ballData(:,2),movement_time,'linear');

figure(700)
maxlag=500;
xc_1=xcorr(detrend(movementData.locDataInterp(1:(end-100),2))', detrend(movementData.targetPosition(1:(end-100),1))',maxlag,'coeff');
xc_2=xcorr(detrend(movementData.locDataInterp(1:(end-100),2))', detrend(movementData.targetPosition(1:(end-100),2))',maxlag,'coeff');
hold off
plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_1)
hold on
plot(movementData.secondsPerFrame*(-maxlag:maxlag),xc_2)

best_lag=round(1/movementData.secondsPerFrame)%about 1 second
