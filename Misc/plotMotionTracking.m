%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:    plotMotionTracking
%
% FUNCTION:         plotMotionTracking(matFileName)
%
% DESCRIPTION:      Generates plots of motion tracking data from saved .mat file
%
% INPUT:
%
% VARIABLES:
%
% OUTPUT:
%
% FUNCTIONS USED:
%
% LIBARIES USED:
%
% NOTES:
%
% WRITTEN BY:       Spencer Garborg 1/28/19
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotMotionTracking(matFileName)

close all;

load(matFileName);

subplot(2,1,1)
plot(1:size(movementData.MoveDist,1),movementData.MoveDist(:,1),'r')
title('Object Movement Between Frames')
xlabel('Frame')
ylabel('X Movement (\mum)')
grid on
subplot(2,1,2)
plot(1:size(movementData.MoveDist,1),movementData.MoveDist(:,2),'b')
title('Object Movement Between Frames')
xlabel('Frame')
ylabel('Y Movement (\mum)')
grid on
% subplot(3,1,3)
% plot(1:size(ballData,1),ballData,'.k')

figure(2)
subplot(2,1,1)
plot(1:size(movementData.TargetPosition,1),movementData.TargetPosition(:,1),'r')
title('Object Position per Frame')
xlabel('Frame')
ylabel('X Position (\mum)')
grid on
subplot(2,1,2)
plot(1:size(movementData.TargetPosition,1),movementData.TargetPosition(:,2),'b')
title('Object Position per Frame')
xlabel('Frame')
ylabel('Y Position (\mum)')
grid on
% subplot(3,1,3)
% plot(1:size(ballData,1),ballData,'.k')

figure(3)
k = convhull(movementData.TargetPosition(:,1),movementData.TargetPosition(:,2));
plot(movementData.TargetPosition(k,1),movementData.TargetPosition(k,2),'b',movementData.TargetPosition(:,1),movementData.TargetPosition(:,2),'k');
maxVal = ceil(max(max(abs(movementData.TargetPosition)))/10)*10;
axis equal square
axis([-maxVal maxVal -maxVal maxVal])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title('Position of Target Object')
xlabel('\mum')
ylabel('\mum')
end