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
subplot(4,1,1)
plot((1:size(movementData.moveDist,1))*movementData.secondsPerFrame,medfilt1(movementData.moveDist(:,1),medFiltSize),'r')
title(['\fontsize{20pt}\bf{Object Movement Between Frames}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
xlabel('Time (s)')
ylabel('X Movement (\mum)')
grid on
axis([1 size(movementData.moveDist,1)*movementData.secondsPerFrame floor(min(medfilt1([movementData.moveDist(:,1);movementData.moveDist(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.moveDist(:,1);movementData.moveDist(:,2)],medFiltSize)))])
subplot(4,1,2)
plot((1:size(movementData.moveDist,1))*movementData.secondsPerFrame,medfilt1(movementData.moveDist(:,2),medFiltSize),'b')
xlabel('Time (s)')
ylabel('Y Movement (\mum)')
grid on
axis([1 size(movementData.moveDist,1)*movementData.secondsPerFrame floor(min(medfilt1([movementData.moveDist(:,1);movementData.moveDist(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.moveDist(:,1);movementData.moveDist(:,2)],medFiltSize)))])
subplot(4,1,3)
plot(movementData.ballData(:,1),movementData.ballData(:,2),'k')
title('\fontsize{20pt}\bf{Ball Movement}')
xlabel('Time (s)')
ylabel('m/s')
grid on
axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) 0 ceil(max(movementData.ballData(:,2))*10)/10])
subplot(4,1,4)
if all(movementData.emgData(:,2) == 0)
    title('\fontsize{20pt}\bf{No EMG Data}')
else
    semilogy(movementData.emgData(:,1),movementData.emgData(:,2),'k')
    title('\fontsize{20pt}\bf{EMG}')
    xlabel('Time (s)')
    ylabel('Amplitude (a.u.)')
    grid on
    axis([min(movementData.emgData(:,1)) max(movementData.emgData(:,1)) 0 ceil(max(movementData.emgData(:,2)))])
end

h(2) = figure('Color','White');
subplot(2,1,1)
plot((1:length(movementData.velocity))*movementData.secondsPerFrame,medfilt1(movementData.velocity,medFiltSize),'r')
title(['\fontsize{20pt}\bf{Object Velocity Between Frames}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
xlabel('Time (s)')
ylabel('Velocity (\mum/s)')
grid on
axis([1 size(movementData.velocity,1)*movementData.secondsPerFrame 0 ceil(max(medfilt1(movementData.velocity(:,1),medFiltSize))/10)*10])
subplot(2,1,2)
plot(movementData.ballData(:,1),movementData.ballData(:,2),'k')
title('\fontsize{20pt}\bf{Ball Movement}')
xlabel('Time (s)')
ylabel('m/s')
grid on
axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) 0 ceil(max(movementData.ballData(:,2))*10)/10])

h(3) = figure('Color','White');
x1 = subplot(4,1,1);
plot((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame,medfilt1(movementData.targetPosition(:,1),medFiltSize),'r')
if contains(matFileName,'combined')
    hold on;
    f = fill([(1:size(movementData.targetPosition,1))*movementData.secondsPerFrame flip((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame)],movementData.cIntFillPtsX,'r','Linestyle','none');
    set(f,'facea',[.2]);
%     plot((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame,movementData.targetPositionSGF(:,1),'k')
    hold off
end
title(['\fontsize{20pt}\bf{Object Position per Frame}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
xlabel('Time (s)')
ylabel('X Position (\mum)')
grid on
if movementData.hemisphere == 1
    text(0,ceil(max(medfilt1([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)],medFiltSize))),'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(0,floor(min(medfilt1([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)],medFiltSize))),'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
else
    text(0,ceil(max(medfilt1([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)],medFiltSize))),'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(0,floor(min(medfilt1([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)],medFiltSize))),'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
end
axis([1 size(movementData.targetPosition,1)*movementData.secondsPerFrame floor(min(medfilt1([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)],medFiltSize)))])
% set(x1,'Position',[.05, .68, .9, .23])
x2 = subplot(4,1,2);
plot((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame,-1*medfilt1(movementData.targetPosition(:,2),medFiltSize),'b')
if contains(matFileName,'combined')
    hold on;
    f = fill([(1:size(movementData.targetPosition,1))*movementData.secondsPerFrame flip((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame)],-1*movementData.cIntFillPtsY,'r','Linestyle','none');
    set(f,'facea',[.2]);
%     plot((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame,movementData.targetPositionSGF(:,2),'k')
    hold off
end
xlabel('Time (s)')
ylabel('Y Position (\mum)')
grid on
text(0,-floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))),'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
text(0,-ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))),'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
axis([1 size(movementData.targetPosition,1)*movementData.secondsPerFrame floor(min(medfilt1([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)],medFiltSize)))])
% set(x2,'Position',[.05, .39, .9, .23])
x3 = subplot(4,1,3);
plot(movementData.ballData(:,1),movementData.ballData(:,2),'k')
title('\fontsize{20pt}\bf{Ball Movement}')
xlabel('Time (s)')
ylabel('m/s')
grid on
axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) 0 ceil(max(movementData.ballData(:,2))*10)/10])
% set(x3,'Position',[.05, .06, .9, .23])
x4 = subplot(4,1,4);
% if all(movementData.emgData(:,2) == 0)
%     title('\fontsize{20pt}\bf{No EMG Data}')
% else
%     plot(movementData.emgData(:,1),movementData.emgData(:,2),'k')
%     title('\fontsize{20pt}\bf{EMG}')
%     xlabel('Time (s)')
%     ylabel('Power')
%     grid on
%     xlim([min(movementData.emgData(:,1)) max(movementData.emgData(:,1))])
%     % axis([min(movementData.emgData(:,1)) max(movementData.emgData(:,1)) 0.9 ceil(max(movementData.emgData(2:end,2))*10)/10])
% end
if all(movementData.emgData(:,2) == 0)
    title('\fontsize{20pt}\bf{No EMG Data}')
else
    semilogy(movementData.emgData(:,1),movementData.emgData(:,2),'k')
    title('\fontsize{20pt}\bf{EMG}')
    xlabel('Time (s)')
    ylabel('Amplitude (a.u.)')
    grid on
    axis([min(movementData.emgData(:,1)) max(movementData.emgData(:,1)) 0 ceil(max(movementData.emgData(:,2)))])
end

h(4) = figure('Color','White');
medFiltData = [medfilt1(movementData.targetPosition(:,1),medFiltSize),-1*medfilt1(movementData.targetPosition(:,2),medFiltSize)];
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
maxMoveData = ceil(max([medFiltData(:,1);medFiltData(:,2)]));
if movementData.hemisphere == 1
    text(maxMoveData,0,'Lateral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(-maxMoveData,0,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    text(0,maxMoveData,'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(0,-maxMoveData,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
else
    text(maxMoveData,0,'Medial','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(-maxMoveData,0,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    text(0,maxMoveData,'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(0,-maxMoveData,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
end

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
maxMoveData = ceil(max([medFiltData(:,1);medFiltData(:,2)]));
if movementData.hemisphere == 1
    text(maxMoveData,0,'Lateral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(-maxMoveData,0,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    text(0,maxMoveData,'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(0,-maxMoveData,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
else
    text(maxMoveData,0,'Medial','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(-maxMoveData,0,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    text(0,maxMoveData,'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(0,-maxMoveData,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
end

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

h(8) = figure('Color','White');
motionEventsLocationsX = [];
motionEventsLocationsY = [];
if size(movementData.motionEvents,1) == 0
    title('No Ball Motion Events To Plot')
else
    for n = 1:size(movementData.motionEvents,1)
        motionVectorX = movementData.targetPosition(movementData.motionEvents(n,4):movementData.motionEvents(n,6),1);
        motionVectorY = movementData.targetPosition(movementData.motionEvents(n,4):movementData.motionEvents(n,6),2);
        if n > 1
            if length(motionVectorX) > size(motionEventsLocationsX,2)
                motionVectorX = motionVectorX(1:size(motionEventsLocationsX,2));
            elseif length(motionVectorX) < size(motionEventsLocationsX,2)
                motionEventsLocationsX = motionEventsLocationsX(:,1:length(motionVectorX));
            end
            if length(motionVectorY) > size(motionEventsLocationsY,2)
                motionVectorY = motionVectorY(1:size(motionEventsLocationsY,2));
            elseif length(motionVectorY) < size(motionEventsLocationsY,2)
                motionEventsLocationsY = motionEventsLocationsY(:,1:length(motionVectorY));
            end
        end
        motionEventsLocationsX(end+1,:) = medfilt1(motionVectorX-motionVectorX(1),medFiltSize);
        motionEventsLocationsY(end+1,:) = medfilt1(motionVectorY-motionVectorY(1),medFiltSize);
    end
    [meanX,cIntFillPtsX] = getCIntMeanAndFillPts(motionEventsLocationsX,90);
    [meanY,cIntFillPtsY] = getCIntMeanAndFillPts(motionEventsLocationsY,90);
    timeVecX = linspace(round(movementData.motionEvents(1,1)-movementData.motionEvents(1,2)),round(movementData.motionEvents(1,3)-movementData.motionEvents(1,2)),length(meanX));
    timeVecY = linspace(round(movementData.motionEvents(1,1)-movementData.motionEvents(1,2)),round(movementData.motionEvents(1,3)-movementData.motionEvents(1,2)),length(meanY));
    
    subplot(3,1,1)
    plot((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame,medfilt1(movementData.targetPosition(:,1),medFiltSize),'b')
    moveMin = floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)));
    moveMax = ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)));
    hold on
    for n = 1:size(motionEventsLocationsX,1)
        plot([movementData.motionEvents(n,2) movementData.motionEvents(n,2)],[moveMin moveMax],'k--')
        f = fill([movementData.motionEvents(n,1) movementData.motionEvents(n,2) movementData.motionEvents(n,2) movementData.motionEvents(n,1)],[moveMin moveMin moveMax moveMax],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([movementData.motionEvents(n,2) movementData.motionEvents(n,3) movementData.motionEvents(n,3) movementData.motionEvents(n,2)],[moveMin moveMin moveMax moveMax],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
%     title(['\fontsize{20pt}\bf{X Position - Locomotion Events}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Time (s)')
    ylabel('X Position (\mum)')
    title('\fontsize{20pt}\bf{Motion During Locomotion Events}')
    grid on
    if movementData.hemisphere == 1
        text(0,ceil(max(medfilt1([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)],medFiltSize))),'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
        text(0,floor(min(medfilt1([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)],medFiltSize))),'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    else
        text(0,ceil(max(medfilt1([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)],medFiltSize))),'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
        text(0,floor(min(medfilt1([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)],medFiltSize))),'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    end
    axis([0 size(movementData.targetPosition,1)*movementData.secondsPerFrame floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)))])
    subplot(3,1,2)
    plot((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame,-1*medfilt1(movementData.targetPosition(:,2),medFiltSize),'b')
    moveMin = floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)));
    moveMax = ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)));
    hold on
    for n = 1:size(motionEventsLocationsX,1)
        plot([movementData.motionEvents(n,2) movementData.motionEvents(n,2)],[moveMin moveMax],'k--')
        f = fill([movementData.motionEvents(n,1) movementData.motionEvents(n,2) movementData.motionEvents(n,2) movementData.motionEvents(n,1)],[moveMin moveMin moveMax moveMax],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([movementData.motionEvents(n,2) movementData.motionEvents(n,3) movementData.motionEvents(n,3) movementData.motionEvents(n,2)],[moveMin moveMin moveMax moveMax],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
%     title(['\fontsize{20pt}\bf{Y Position - Locomotion Events}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Time (s)')
    ylabel('Y Position (\mum)')
    grid on
    text(0,-floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))),'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(0,-ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))),'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    axis([0 size(movementData.targetPosition,1)*movementData.secondsPerFrame floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)))])
    subplot(3,1,3)
    plot(movementData.ballData(:,1),movementData.ballData(:,2),'k')
    title('\fontsize{20pt}\bf{Ball Movement}')
    xlabel('Time (s)')
    ylabel('Movement')
    grid on
    axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) -1 ceil(max(movementData.ballData(:,2)))])
    hold on
    for n = 1:size(motionEventsLocationsX,1)
        plot([movementData.motionEvents(n,2) movementData.motionEvents(n,2)],[moveMin moveMax],'k--')
        f = fill([movementData.motionEvents(n,1) movementData.motionEvents(n,2) movementData.motionEvents(n,2) movementData.motionEvents(n,1)],[moveMin moveMin moveMax moveMax],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([movementData.motionEvents(n,2) movementData.motionEvents(n,3) movementData.motionEvents(n,3) movementData.motionEvents(n,2)],[moveMin moveMin moveMax moveMax],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
end

h(9) = figure('Color','White');
if size(movementData.motionEvents,1) == 0
    title('No Ball Motion Events To Plot')
else
    subplot(2,1,1)
    maxMeanVal = max(abs([meanX meanY]));
    plot(timeVecX,meanX,'k')
    hold on
    %     f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
    %     set(f,'facea',[.2]);
    plot([0 0],[-maxMeanVal maxMeanVal],'r')
    hold off
    title('\fontsize{20pt}\bf{Mean Motion During Locomotion Events}')
    xlabel('Time (s)')
    ylabel('X Position (\mum)')
    ylim([-maxMeanVal maxMeanVal])
    grid on
    subplot(2,1,2)
    plot(timeVecY,-1*meanY,'k')
    hold on
    %     f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
    %     set(f,'facea',[.2]);
    plot([0 0],[-maxMeanVal maxMeanVal],'r')
    hold off
    xlabel('Time (s)')
    ylabel('Y Position (\mum)')
    ylim([-maxMeanVal maxMeanVal])
    grid on
end

h(10) = figure('Color','White');
EMGEventsLocationsX = [];
EMGEventsLocationsY = [];
if size(movementData.EMGEvents,1) == 0
    title('No EMG Events To Plot')
else
    for n = 1:size(movementData.EMGEvents,1)
        motionVectorX = movementData.targetPosition(movementData.EMGEvents(n,4):movementData.EMGEvents(n,6),1);
        motionVectorY = movementData.targetPosition(movementData.EMGEvents(n,4):movementData.EMGEvents(n,6),2);
        if n > 1
            if length(motionVectorX) > size(EMGEventsLocationsX,2)
                motionVectorX = motionVectorX(1:size(EMGEventsLocationsX,2));
            elseif length(motionVectorX) < size(EMGEventsLocationsX,2)
                EMGEventsLocationsX = EMGEventsLocationsX(:,1:length(motionVectorX));
            end
            if length(motionVectorY) > size(EMGEventsLocationsY,2)
                motionVectorY = motionVectorY(1:size(EMGEventsLocationsY,2));
            elseif length(motionVectorY) < size(EMGEventsLocationsY,2)
                EMGEventsLocationsY = EMGEventsLocationsY(:,1:length(motionVectorY));
            end
        end
        EMGEventsLocationsX(end+1,:) = medfilt1(motionVectorX-motionVectorX(1),medFiltSize);
        EMGEventsLocationsY(end+1,:) = medfilt1(motionVectorY-motionVectorY(1),medFiltSize);
    end
    [meanX,cIntFillPtsX] = getCIntMeanAndFillPts(EMGEventsLocationsX,90);
    [meanY,cIntFillPtsY] = getCIntMeanAndFillPts(EMGEventsLocationsY,90);
    timeVecX = linspace(round(movementData.EMGEvents(1,1)-movementData.EMGEvents(1,2)),round(movementData.EMGEvents(1,3)-movementData.EMGEvents(1,2)),length(meanX));
    timeVecY = linspace(round(movementData.EMGEvents(1,1)-movementData.EMGEvents(1,2)),round(movementData.EMGEvents(1,3)-movementData.EMGEvents(1,2)),length(meanY));
    
    subplot(3,1,1)
    plot((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame,medfilt1(movementData.targetPosition(:,1),medFiltSize),'b')
    moveMin = floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)));
    moveMax = ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)));
    hold on
    for n = 1:size(EMGEventsLocationsY,1)
        plot([movementData.EMGEvents(n,2) movementData.EMGEvents(n,2)],[moveMin moveMax],'k--')
        f = fill([movementData.EMGEvents(n,1) movementData.EMGEvents(n,2) movementData.EMGEvents(n,2) movementData.EMGEvents(n,1)],[moveMin moveMin moveMax moveMax],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([movementData.EMGEvents(n,2) movementData.EMGEvents(n,3) movementData.EMGEvents(n,3) movementData.EMGEvents(n,2)],[moveMin moveMin moveMax moveMax],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
%     title(['\fontsize{20pt}\bf{X Position - EMG Events}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Time (s)')
    ylabel('X Position (\mum)')
    title('\fontsize{20pt}\bf{Motion During EMG Events}')
    grid on
    if movementData.hemisphere == 1
        text(0,ceil(max(medfilt1([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)],medFiltSize))),'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
        text(0,floor(min(medfilt1([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)],medFiltSize))),'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    else
        text(0,ceil(max(medfilt1([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)],medFiltSize))),'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
        text(0,floor(min(medfilt1([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)],medFiltSize))),'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    end
    axis([0 size(movementData.targetPosition,1)*movementData.secondsPerFrame floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)))])
    subplot(3,1,2)
    plot((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame,-1*medfilt1(movementData.targetPosition(:,2),medFiltSize),'b')
    moveMin = floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)));
    moveMax = ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)));
    hold on
    for n = 1:size(EMGEventsLocationsY,1)
        plot([movementData.EMGEvents(n,2) movementData.EMGEvents(n,2)],[moveMin moveMax],'k--')
        f = fill([movementData.EMGEvents(n,1) movementData.EMGEvents(n,2) movementData.EMGEvents(n,2) movementData.EMGEvents(n,1)],[moveMin moveMin moveMax moveMax],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([movementData.EMGEvents(n,2) movementData.EMGEvents(n,3) movementData.EMGEvents(n,3) movementData.EMGEvents(n,2)],[moveMin moveMin moveMax moveMax],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
    %     title(['\fontsize{20pt}\bf{Y Position - EMG Events}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Time (s)')
    ylabel('Y Position (\mum)')
    grid on
    text(0,-floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))),'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(0,-ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))),'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    axis([0 size(movementData.targetPosition,1)*movementData.secondsPerFrame floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)))])
    subplot(3,1,3)
    %     plot(movementData.emgData(:,1),movementData.emgData(:,2),'k')
    plot(movementData.emgData(:,1),movementData.emgData(:,2),'k')
    title('\fontsize{20pt}\bf{EMG}')
    xlabel('Time (s)')
    ylabel('Amplitude (a.u.)')
    grid on
    axis([min(movementData.emgData(:,1)) max(movementData.emgData(:,1)) .8 ceil(max(movementData.emgData(2:end,2)))])
    hold on
    for n = 1:size(EMGEventsLocationsY,1)
        plot([movementData.EMGEvents(n,2) movementData.EMGEvents(n,2)],[.8 ceil(max(movementData.emgData(:,2)))],'k--')
        f = fill([movementData.EMGEvents(n,1) movementData.EMGEvents(n,2) movementData.EMGEvents(n,2) movementData.EMGEvents(n,1)],[0 0 ceil(max(movementData.emgData(:,2))) ceil(max(movementData.emgData(:,2)))],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([movementData.EMGEvents(n,2) movementData.EMGEvents(n,3) movementData.EMGEvents(n,3) movementData.EMGEvents(n,2)],[0 0 ceil(max(movementData.emgData(:,2))) ceil(max(movementData.emgData(:,2)))],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
end

h(11) = figure('Color','White');
if size(movementData.EMGEvents,1) == 0
    title('No EMG Events To Plot')
else
    subplot(2,1,1)
    maxMeanVal = max(abs([meanX meanY]));
    plot(timeVecX,meanX,'k')
    hold on
    %     f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
    %     set(f,'facea',[.2]);
    plot([0 0],[-maxMeanVal maxMeanVal],'r')
    hold off
    title('\fontsize{20pt}\bf{Mean Motion During EMG Events}')
    xlabel('Time (s)')
    ylabel('X Position (\mum)')
    ylim([-maxMeanVal maxMeanVal])
    grid on
    subplot(2,1,2)
    plot(timeVecY,-1*meanY,'k')
    hold on
    %     f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
    %     set(f,'facea',[.2]);
    plot([0 0],[-maxMeanVal maxMeanVal],'r')
    hold off
    xlabel('Time (s)')
    ylabel('Y Position (\mum)')
    ylim([-maxMeanVal maxMeanVal])
    grid on
end

h(12) = figure('Color','White');
EMGNoMotionEventsLocationsX = [];
EMGNoMotionEventsLocationsY = [];
if size(movementData.EMGNoMotionEvents,1) == 0
    title('No EMG Without Ball Motion Events To Plot')
else
    for n = 1:size(movementData.EMGNoMotionEvents,1)
        motionVectorX = movementData.targetPosition(movementData.EMGNoMotionEvents(n,4):movementData.EMGNoMotionEvents(n,6),1);
        motionVectorY = movementData.targetPosition(movementData.EMGNoMotionEvents(n,4):movementData.EMGNoMotionEvents(n,6),2);
        if n > 1
            if length(motionVectorX) > size(EMGNoMotionEventsLocationsX,2)
                motionVectorX = motionVectorX(1:size(EMGNoMotionEventsLocationsX,2));
            elseif length(motionVectorX) < size(EMGNoMotionEventsLocationsX,2)
                EMGNoMotionEventsLocationsX = EMGNoMotionEventsLocationsX(:,1:length(motionVectorX));
            end
            if length(motionVectorY) > size(EMGNoMotionEventsLocationsY,2)
                motionVectorY = motionVectorY(1:size(EMGNoMotionEventsLocationsY,2));
            elseif length(motionVectorY) < size(EMGNoMotionEventsLocationsY,2)
                EMGNoMotionEventsLocationsY = EMGNoMotionEventsLocationsY(:,1:length(motionVectorY));
            end
        end
        EMGNoMotionEventsLocationsX(end+1,:) = medfilt1(motionVectorX-motionVectorX(1),medFiltSize);
        EMGNoMotionEventsLocationsY(end+1,:) = medfilt1(motionVectorY-motionVectorY(1),medFiltSize);
    end
    [meanX,cIntFillPtsX] = getCIntMeanAndFillPts(EMGNoMotionEventsLocationsX,90);
    [meanY,cIntFillPtsY] = getCIntMeanAndFillPts(EMGNoMotionEventsLocationsY,90);
    timeVecX = linspace(round(movementData.EMGNoMotionEvents(1,1)-movementData.EMGNoMotionEvents(1,2)),round(movementData.EMGNoMotionEvents(1,3)-movementData.EMGNoMotionEvents(1,2)),length(meanX));
    timeVecY = linspace(round(movementData.EMGNoMotionEvents(1,1)-movementData.EMGNoMotionEvents(1,2)),round(movementData.EMGNoMotionEvents(1,3)-movementData.EMGNoMotionEvents(1,2)),length(meanY));
    
    subplot(4,1,1)
    plot((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame,medfilt1(movementData.targetPosition(:,1),medFiltSize),'b')
    moveMin = floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)));
    moveMax = ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)));
    hold on
    for n = 1:size(EMGNoMotionEventsLocationsX,1)
        plot([movementData.EMGNoMotionEvents(n,2) movementData.EMGNoMotionEvents(n,2)],[moveMin moveMax],'k--')
        f = fill([movementData.EMGNoMotionEvents(n,1) movementData.EMGNoMotionEvents(n,2) movementData.EMGNoMotionEvents(n,2) movementData.EMGNoMotionEvents(n,1)],[moveMin moveMin moveMax moveMax],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([movementData.EMGNoMotionEvents(n,2) movementData.EMGNoMotionEvents(n,3) movementData.EMGNoMotionEvents(n,3) movementData.EMGNoMotionEvents(n,2)],[moveMin moveMin moveMax moveMax],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
%     title(['\fontsize{20pt}\bf{X Position - EMG (No Locomotion) Events}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Time (s)')
    ylabel('X Position (\mum)')
    title('\fontsize{20pt}\bf{Motion During EMG Events With No Locomotion}')
    grid on
    if movementData.hemisphere == 1
        text(0,ceil(max(medfilt1([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)],medFiltSize))),'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
        text(0,floor(min(medfilt1([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)],medFiltSize))),'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    else
        text(0,ceil(max(medfilt1([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)],medFiltSize))),'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
        text(0,floor(min(medfilt1([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)],medFiltSize))),'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    end
    %     axis([1 size(movementData.targetPosition,1) floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)))])
    subplot(4,1,2)
    plot((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame,-1*medfilt1(movementData.targetPosition(:,2),medFiltSize),'b')
    moveMin = floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)));
    moveMax = ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)));
    hold on
    for n = 1:size(EMGNoMotionEventsLocationsX,1)
        plot([movementData.EMGNoMotionEvents(n,2) movementData.EMGNoMotionEvents(n,2)],[moveMin moveMax],'k--')
        f = fill([movementData.EMGNoMotionEvents(n,1) movementData.EMGNoMotionEvents(n,2) movementData.EMGNoMotionEvents(n,2) movementData.EMGNoMotionEvents(n,1)],[moveMin moveMin moveMax moveMax],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([movementData.EMGNoMotionEvents(n,2) movementData.EMGNoMotionEvents(n,3) movementData.EMGNoMotionEvents(n,3) movementData.EMGNoMotionEvents(n,2)],[moveMin moveMin moveMax moveMax],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
%     title(['\fontsize{20pt}\bf{Y Position - EMG (No Locomotion) Events}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Time (s)')
    ylabel('Y Position (\mum)')
    grid on
    text(0,-floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))),'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(0,-ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))),'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
%     axis([1 size(movementData.targetPosition,1) floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)))])
    subplot(4,1,3)
    plot(movementData.ballData(:,1),movementData.ballData(:,2),'k')
    title('Ball Movement')
    xlabel('Time (s)')
    ylabel('Movement')
    grid on
    axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) -1 ceil(max(movementData.ballData(:,2)))])
    hold on
    for n = 1:size(EMGNoMotionEventsLocationsX,1)
        plot([movementData.EMGNoMotionEvents(n,2) movementData.EMGNoMotionEvents(n,2)],[moveMin moveMax],'k--')
        f = fill([movementData.EMGNoMotionEvents(n,1) movementData.EMGNoMotionEvents(n,2) movementData.EMGNoMotionEvents(n,2) movementData.EMGNoMotionEvents(n,1)],[moveMin moveMin moveMax moveMax],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([movementData.EMGNoMotionEvents(n,2) movementData.EMGNoMotionEvents(n,3) movementData.EMGNoMotionEvents(n,3) movementData.EMGNoMotionEvents(n,2)],[moveMin moveMin moveMax moveMax],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
    subplot(4,1,4)
    %     plot(movementData.emgData(:,1),movementData.emgData(:,2),'k')
    plot(movementData.emgData(:,1),movementData.emgData(:,2),'k')
    title('EMG')
    xlabel('Time (s)')
    ylabel('Amplitude (a.u.)')
    grid on
    axis([min(movementData.emgData(:,1)) max(movementData.emgData(:,1)) .8 ceil(max(movementData.emgData(2:end,2)))])
    hold on
    for n = 1:size(EMGNoMotionEventsLocationsX,1)
        plot([movementData.EMGNoMotionEvents(n,2) movementData.EMGNoMotionEvents(n,2)],[moveMin moveMax],'k--')
        f = fill([movementData.EMGNoMotionEvents(n,1) movementData.EMGNoMotionEvents(n,2) movementData.EMGNoMotionEvents(n,2) movementData.EMGNoMotionEvents(n,1)],[moveMin moveMin moveMax moveMax],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([movementData.EMGNoMotionEvents(n,2) movementData.EMGNoMotionEvents(n,3) movementData.EMGNoMotionEvents(n,3) movementData.EMGNoMotionEvents(n,2)],[moveMin moveMin moveMax moveMax],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
end

h(13) = figure('Color','White');
if size(movementData.EMGNoMotionEvents,1) == 0
    title('No EMG Without Ball Motion Events To Plot')
else
    subplot(2,1,1)
    maxMeanVal = max(abs([meanX meanY]));
    plot(timeVecX,meanX,'k')
    hold on
    %     f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
    %     set(f,'facea',[.2]);
    plot([0 0],[-maxMeanVal maxMeanVal],'r')
    hold off
    title('\fontsize{20pt}\bf{Mean Motion During EMG Events With No Locomotion}')
    xlabel('Time (s)')
    ylabel('X Position (\mum)')
    ylim([-maxMeanVal maxMeanVal])
    grid on
    subplot(2,1,2)
    plot(timeVecY,-1*meanY,'k')
    hold on
    %     f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
    %     set(f,'facea',[.2]);
    plot([0 0],[-maxMeanVal maxMeanVal],'r')
    hold off
    xlabel('Time (s)')
    ylabel('Y Position (\mum)')
    ylim([-maxMeanVal maxMeanVal])
    grid on
end

motionVec = pcaMotionAnalysis(medFiltData);
h(14) = figure('Color','White');
scatter(medFiltData(:,1),medFiltData(:,2),10)
hold on
drawArrow([0;0],[motionVec(1);motionVec(2)]);
hold off
axis equal square
axis([-maxMoveData maxMoveData -maxMoveData maxMoveData])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title(['\fontsize{20pt}\bf{Position of Brain in Skull}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
xlabel('\mum')
ylabel('\mum')
if movementData.hemisphere == 1
    text(maxMoveData,0,'Lateral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(-maxMoveData,0,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    text(0,maxMoveData,'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(0,-maxMoveData,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
else
    text(maxMoveData,0,'Medial','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(-maxMoveData,0,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    text(0,maxMoveData,'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    text(0,-maxMoveData,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
end
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
if contains(matFileName,'combined')
    matFileName = [matFileName(1:end-23) '_outputPlots_combined.fig'];
else
    fileNum = matFileName(end-4);
    matFileName = [matFileName(1:end-16) '_outputPlots_' fileNum '.fig'];
end
savefig(h,matFileName);
