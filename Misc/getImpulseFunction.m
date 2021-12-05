close all
ballData = load('H:\21-05-24_MouseExp\210524_001.txt');
load('H:\21-05-24_MouseExp\210524_001_processed_Layer1_4.mat');
[speeddata, speedaccBinary, speedaveBinary] = speedProcess(20, convertBallVoltToMPS(ballData(:,2)), size(movementData.targetPosition,1));

out = OXY_HRF(speedaveBinary,movementData.targetPosition(:,1),20);

figure(1)
plot(out.HRF_time,out.HRF)

out = OXY_HRF(speedaveBinary,movementData.targetPosition(:,2),20);
figure(2)
plot(out.HRF_time,out.HRF)