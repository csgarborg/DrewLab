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

% Generate plots
subtitle = [num2str(1/movementData.secondsPerFrame) ' Frames/s, ' num2str(movementData.secondsPerFrame*(diff(movementData.frames)+1)) ' Seconds, ' num2str(movementData.objMag*movementData.digMag) 'x Magnification (' num2str(movementData.objMag) 'x Objective, ' num2str(movementData.digMag) 'x Digital), Turnabout = ' num2str(movementData.turnabout)];
h(1) = figure('Color','White');
subplot(3,1,1)
plot(1:size(movementData.moveDist,1),movementData.moveDist(:,1),'r')
title(['\fontsize{20pt}\bf{Object Movement Between Frames}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
xlabel('Frame')
ylabel('X Movement (\mum)')
grid on
axis([1 size(movementData.moveDist,1) floor(min([movementData.moveDist(:,1);movementData.moveDist(:,2)])/10)*10 ceil(max([movementData.moveDist(:,1);movementData.moveDist(:,2)])/10)*10])
subplot(3,1,2)
plot(1:size(movementData.moveDist,1),movementData.moveDist(:,2),'b')
xlabel('Frame')
ylabel('Y Movement (\mum)')
grid on
axis([1 size(movementData.moveDist,1) floor(min([movementData.moveDist(:,1);movementData.moveDist(:,2)])/10)*10 ceil(max([movementData.moveDist(:,1);movementData.moveDist(:,2)])/10)*10])
subplot(3,1,3)
plot(movementData.ballData(:,1),movementData.ballData(:,2),'k')
title('\fontsize{20pt}\bf{Ball Movement}')
xlabel('Time (s)')
ylabel('Movement')
grid on
axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) -1 ceil(max(movementData.ballData(:,2)))])

h(2) = figure('Color','White');
subplot(2,1,1)
plot(1:length(movementData.velocity),movementData.velocity,'r')
title(['\fontsize{20pt}\bf{Object Velocity Between Frames}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
xlabel('Frame')
ylabel('Velocity (\mum/s)')
grid on
axis([1 size(movementData.velocity,1) 0 ceil(max(movementData.velocity(:,1))/10)*10])
subplot(2,1,2)
plot(movementData.ballData(:,1),movementData.ballData(:,2),'k')
title('\fontsize{20pt}\bf{Ball Movement}')
xlabel('Time (s)')
ylabel('Movement')
grid on
axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) -1 ceil(max(movementData.ballData(:,2)))])

h(3) = figure('Color','White');
subplot(3,1,1)
plot(1:size(movementData.targetPosition,1),movementData.targetPosition(:,1),'r')
title(['\fontsize{20pt}\bf{Object Position per Frame}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
xlabel('Frame')
ylabel('X Position (\mum)')
grid on
axis([1 size(movementData.targetPosition,1) floor(min([movementData.targetPosition(:,1);movementData.targetPosition(:,2)])/10)*10 ceil(max([movementData.targetPosition(:,1);movementData.targetPosition(:,2)])/10)*10])
subplot(3,1,2)
plot(1:size(movementData.targetPosition,1),movementData.targetPosition(:,2),'b')
xlabel('Frame')
ylabel('Y Position (\mum)')
grid on
axis([1 size(movementData.targetPosition,1) floor(min([movementData.targetPosition(:,1);movementData.targetPosition(:,2)])/10)*10 ceil(max([movementData.targetPosition(:,1);movementData.targetPosition(:,2)])/10)*10])
subplot(3,1,3)
plot(movementData.ballData(:,1),movementData.ballData(:,2),'k')
title('\fontsize{20pt}\bf{Ball Movement}')
xlabel('Time (s)')
ylabel('Movement')
grid on
axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) -1 ceil(max(movementData.ballData(:,2)))])

h(4) = figure('Color','White');
k = convhull(movementData.targetPosition(:,1),movementData.targetPosition(:,2));
plot(movementData.targetPosition(k,1),movementData.targetPosition(k,2),'b',movementData.targetPosition(:,1),movementData.targetPosition(:,2),'k');
maxVal = ceil(max(max(abs(movementData.targetPosition)))/10)*10;
axis equal square
axis([-maxVal maxVal -maxVal maxVal])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title(['\fontsize{20pt}\bf{Position of Target Object}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
xlabel('\mum')
ylabel('\mum')

h(5) = figure('Color','White');
hist3(movementData.targetPosition,'CdataMode','auto','Nbins',[40 40]);
colorbar;
view(2);
axis equal square
axis([-maxVal maxVal -maxVal maxVal])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title(['\fontsize{20pt}\bf{Position of Target Object Histogram}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
xlabel('\mum')
ylabel('\mum')

h(6) = figure('Color','White');
hist(movementData.targetPosition(:,1),500);
title(['\fontsize{20pt}\bf{Position of Target Object Histogram (X)}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
xlabel('\mum')
ylabel('Number of Data Points')

h(7) = figure('Color','White');
hist(movementData.targetPosition(:,2),500);
title(['\fontsize{20pt}\bf{Position of Target Object Histogram (Y)}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
xlabel('\mum')
ylabel('Number of Data Points')

% h(5) = figure('Color','White');
% [uxy,~,idx] = unique([movementData.targetPosition(:,1),movementData.targetPosition(:,2)],'rows');
% szscale = histc(idx,unique(idx));
% [x1,y1] = meshgrid(-maxVal:.05:maxVal,-maxVal:.05:maxVal);
% z1 = griddata(uxy(:,1),uxy(:,2),szscale,x1,y1);
% pcolor(x1,y1,log(sqrt(z1)))
% shading interp
% axis equal square
% axis([-maxVal maxVal -maxVal maxVal])
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% title(['\fontsize{20pt}\bf{Position of Target Object}' 10 '\fontsize{10pt}\rm{' subtitle '}'])
% xlabel('\mum')
% ylabel('\mum')

% Save figures to single .fig file
savefig(h,[matFileName(1:end-14) '_outputPlots.fig']);
end