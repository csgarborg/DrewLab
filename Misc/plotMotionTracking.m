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

figure('Color','White')
subplot(3,1,1)
plot(1:size(movementData.moveDist,1),movementData.moveDist(:,1),'r')
title('Object Movement Between Frames')
xlabel('Frame')
ylabel('X Movement (\mum)')
grid on
axis([1 size(movementData.moveDist,1) floor(min(movementData.moveDist(:,1))/10)*10 ceil(max(movementData.moveDist(:,1))/10)*10])
subplot(3,1,2)
plot(1:size(movementData.moveDist,1),movementData.moveDist(:,2),'b')
title('Object Movement Between Frames')
xlabel('Frame')
ylabel('Y Movement (\mum)')
grid on
axis([1 size(movementData.moveDist,1) floor(min(movementData.moveDist(:,2))/10)*10 ceil(max(movementData.moveDist(:,2))/10)*10])
subplot(3,1,3)
plot(movementData.ballData(:,1),movementData.ballData(:,2),'k')
title('Ball Movement')
xlabel('Time (s)')
ylabel('Movement')
grid on
axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) -1 1])

figure('Color','White')
subplot(2,1,1)
plot(1:length(movementData.velocity),movementData.velocity,'r')
title('Object Velocity Between Frames')
xlabel('Frame')
ylabel('Velocity (micrometers/s)')
grid on
axis([1 size(movementData.velocity,1) floor(min(movementData.velocity(:,1))/10)*10 ceil(max(movementData.velocity(:,1))/10)*10])
subplot(2,1,2)
plot(movementData.ballData(:,1),movementData.ballData(:,2),'k')
title('Ball Movement')
xlabel('Time (s)')
ylabel('Movement')
grid on
axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) -1 1])

figure('Color','White')
subplot(3,1,1)
plot(1:size(movementData.targetPosition,1),movementData.targetPosition(:,1),'r')
title('Object Position per Frame')
xlabel('Frame')
ylabel('X Position (\mum)')
grid on
axis([1 ceil(size(movementData.targetPosition,1)/10)*10 floor(min(movementData.targetPosition(:,1))/10)*10 ceil(max(movementData.targetPosition(:,1))/10)*10])
subplot(3,1,2)
plot(1:size(movementData.targetPosition,1),movementData.targetPosition(:,2),'b')
title('Object Position per Frame')
xlabel('Frame')
ylabel('Y Position (\mum)')
grid on
axis([1 ceil(size(movementData.targetPosition,1)/10)*10 floor(min(movementData.targetPosition(:,2))/10)*10 ceil(max(movementData.targetPosition(:,2))/10)*10])
subplot(3,1,3)
plot(movementData.ballData(:,1),movementData.ballData(:,2),'k')
title('Ball Movement')
xlabel('Time (s)')
ylabel('Movement')
grid on
axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) -1 1])

figure('Color','White')
k = convhull(movementData.targetPosition(:,1),movementData.targetPosition(:,2));
plot(movementData.targetPosition(k,1),movementData.targetPosition(k,2),'b',movementData.targetPosition(:,1),movementData.targetPosition(:,2),'k');
maxVal = ceil(max(max(abs(movementData.targetPosition)))/10)*10;
axis equal square
axis([-maxVal maxVal -maxVal maxVal])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title('Position of Target Object')
xlabel('\mum')
ylabel('\mum')
end