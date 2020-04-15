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
medFiltSize = 6;

% Generate plots
subtitle = [num2str(1/movementData.secondsPerFrame) ' Frames/s, ' num2str(movementData.secondsPerFrame*(diff(movementData.frames)+1)) ' Seconds, ' num2str(movementData.objMag*movementData.digMag) 'x Magnification (' num2str(movementData.objMag) 'x Objective, ' num2str(movementData.digMag) 'x Digital), Turnabout = ' num2str(movementData.turnabout)];
movementData.ballData(:,2) = abs(convertBallVoltToMPS(movementData.ballData(:,2)));
h(1) = figure('Color','White');
subplot(3,1,1)
plot(1:size(movementData.moveDist,1),medfilt1(movementData.moveDist(:,1),medFiltSize),'r')
title(['\fontsize{20pt}\bf{Object Movement Between Frames}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
xlabel('Frame')
ylabel('X Movement (\mum)')
grid on
axis([1 size(movementData.moveDist,1) floor(min(medfilt1([movementData.moveDist(:,1);movementData.moveDist(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.moveDist(:,1);movementData.moveDist(:,2)],medFiltSize)))])
subplot(3,1,2)
plot(1:size(movementData.moveDist,1),medfilt1(movementData.moveDist(:,2),medFiltSize),'b')
xlabel('Frame')
ylabel('Y Movement (\mum)')
grid on
axis([1 size(movementData.moveDist,1) floor(min(medfilt1([movementData.moveDist(:,1);movementData.moveDist(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.moveDist(:,1);movementData.moveDist(:,2)],medFiltSize)))])
subplot(3,1,3)
plot(movementData.ballData(:,1),movementData.ballData(:,2),'k')
title('\fontsize{20pt}\bf{Ball Movement}')
xlabel('Time (s)')
ylabel('m/s')
grid on
axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) 0 ceil(max(movementData.ballData(:,2))*10)/10])

h(2) = figure('Color','White');
subplot(2,1,1)
plot(1:length(movementData.velocity),medfilt1(movementData.velocity,medFiltSize),'r')
title(['\fontsize{20pt}\bf{Object Velocity Between Frames}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
xlabel('Frame')
ylabel('Velocity (\mum/s)')
grid on
axis([1 size(movementData.velocity,1) 0 ceil(max(medfilt1(movementData.velocity(:,1),medFiltSize))/10)*10])
subplot(2,1,2)
plot(movementData.ballData(:,1),movementData.ballData(:,2),'k')
title('\fontsize{20pt}\bf{Ball Movement}')
xlabel('Time (s)')
ylabel('m/s')
grid on
axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) 0 ceil(max(movementData.ballData(:,2))*10)/10])

h(3) = figure('Color','White');
x1 = subplot(3,1,1);
plot(1:size(movementData.targetPosition,1),medfilt1(movementData.targetPosition(:,1),medFiltSize),'r')
title(['\fontsize{20pt}\bf{Object Position per Frame}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
xlabel('Frame')
ylabel('X Position (\mum)')
grid on
axis([1 size(movementData.targetPosition,1) floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)))])
set(x1,'Position',[.05, .68, .9, .23])
x2 = subplot(3,1,2);
plot(1:size(movementData.targetPosition,1),medfilt1(movementData.targetPosition(:,2),medFiltSize),'b')
xlabel('Frame')
ylabel('Y Position (\mum)')
grid on
axis([1 size(movementData.targetPosition,1) floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)))])
set(x2,'Position',[.05, .39, .9, .23])
x3 = subplot(3,1,3);
plot(movementData.ballData(:,1),movementData.ballData(:,2),'k')
title('\fontsize{20pt}\bf{Ball Movement}')
xlabel('Time (s)')
ylabel('m/s')
grid on
axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) 0 ceil(max(movementData.ballData(:,2))*10)/10])
set(x3,'Position',[.05, .06, .9, .23])

h(4) = figure('Color','White');
medFiltData = [medfilt1(movementData.targetPosition(:,1),medFiltSize),medfilt1(movementData.targetPosition(:,2),medFiltSize)];
k = convhull(medFiltData(:,1),medFiltData(:,2));
plot(medFiltData(k,1),medFiltData(k,2),'b',medFiltData(:,1),medFiltData(:,2),'k');
maxVal = ceil(max(max(abs(movementData.targetPosition))));
axis equal square
axis([-maxVal maxVal -maxVal maxVal])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title(['\fontsize{20pt}\bf{Position of Target Object}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
xlabel('\mum')
ylabel('\mum')

h(5) = figure('Color','White');
hist3(medFiltData,'CdataMode','auto','Nbins',[15 15]);
colorbar;
view(2);
axis equal square
axis([-maxVal maxVal -maxVal maxVal])
ax = gca;
set(ax,'ColorScale','log')
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title(['\fontsize{20pt}\bf{Position of Target Object Histogram}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
xlabel('\mum')
ylabel('\mum')

h(6) = figure('Color','White');
histogram(medFiltData(:,1),500);
title(['\fontsize{20pt}\bf{Position of Target Object Histogram (X)}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
xlabel('\mum')
ylabel('Number of Data Points')
set(gca, 'YScale', 'log')


h(7) = figure('Color','White');
histogram(medFiltData(:,2),500);
title(['\fontsize{20pt}\bf{Position of Target Object Histogram (Y)}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
xlabel('\mum')
ylabel('Number of Data Points')
set(gca, 'YScale', 'log')

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