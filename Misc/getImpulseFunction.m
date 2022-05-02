close all
ballData = load('H:\21-12-16_MouseExp\211216_002.txt');
load('H:\21-12-16_MouseExp\211216_002_processed_Layer1_combined.mat');
% Fs = 19.775;
Fs = 20;
[~, ~, speedaveBinary] = speedProcess_2PLSM(convertBallVoltToMPS(ballData(:,2)), 10000, Fs, size(movementData.targetPosition,1));
[~, ~, speedaveBinary] = speedProcess_2PLSM(2*pi*0.06*(ballData(:,2))/10, 10000, Fs, size(movementData.targetPosition,1));

[zeroa, poleb, gain] = butter(5,(1/(Fs/2)),'low'); % 10 Hz low pass filter design, 5th order Butterworth
[sos,g] = zp2sos(zeroa,poleb,gain);
latMedFilt = filtfilt(sos,g,movementData.targetPosition(:,1));
rosCauFilt = filtfilt(sos,g,movementData.targetPosition(:,2));

out = OXY_HRF(speedaveBinary,latMedFilt,Fs);
figure(1)
[params, HRF_gamma, R2] = HRF_motion_gammaFit(out.HRF, Fs, 'unconstrained');
plot(out.HRF_time, out.HRF,'k'); hold on; plot(out.HRF_time, HRF_gamma,'r');
title('Medial/Lateral Impulse Response')
xlabel('Time (s)')
ylabel('Medial/Lateral Motion (\mum)')
ylim([-0.04 0.06])

out = OXY_HRF(speedaveBinary,-1*rosCauFilt,Fs);
figure(2)
[params, HRF_gamma, R2] = HRF_motion_gammaFit(out.HRF, Fs, 'unconstrained');
plot(out.HRF_time, out.HRF,'k'); hold on; plot(out.HRF_time, HRF_gamma,'r');
title('Rostral/Caudal Impulse Response')
xlabel('Time (s)')
ylabel('Rostral/Caudal Motion (\mum)')


% Gamma function fit


% Fs = 30.2767;
Fs = 31;
[~, ~, speedaveBinary] = speedProcess_2PLSM(convertBallVoltToMPS(ballData(:,2)), 10000, Fs, size(movementData.emgData,1));

[zeroa, poleb, gain] = butter(5,(10/(Fs/2)),'low'); % 10 Hz low pass filter design, 5th order Butterworth
[sos,g] = zp2sos(zeroa,poleb,gain);
emgFilt = filtfilt(sos,g,movementData.emgData(:,1));

out = OXY_HRF(speedaveBinary,emgFilt(:,2),Fs);
figure(3)
plot(out.HRF_time,out.HRF)