%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:    plotMotionTracking2Layer
%
% FUNCTION:         plotMotionTracking2Layer(matFileName)
%
% DESCRIPTION:      Generates plots of motion tracking data from saved .mat file
%
% INPUT:            First mat file should be the stationary movement data
%                   to stabilize around (skull) and the second should be
%                   the movement data to be analyzed (brain)
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
% WRITTEN BY:       Spencer Garborg 1/10/20
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotMotionTracking2Layer(matFileName1,matFileName2)

close all;


if ~exist('matFileName2','var')
    load(matFileName1);
    medFiltSize = 6;
    
    % Generate plots
    subtitle = [num2str(1/movementData.secondsPerFrame) ' Frames/s, ' num2str(movementData.secondsPerFrame*(diff(movementData.frames)+1)) ' Seconds, ' num2str(movementData.objMag*movementData.digMag) 'x Magnification (' num2str(movementData.objMag) 'x Objective, ' num2str(movementData.digMag) 'x Digital), Turnabout = ' num2str(movementData.turnabout)];
    h(1) = figure('Color','White');
    subplot(4,1,1)
    plot(1:size(movementData.moveDist,1),medfilt1(movementData.moveDist(:,1),medFiltSize),'r')
    title(['\fontsize{20pt}\bf{Object Movement Between Frames}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Frame')
    ylabel('X Movement (\mum)')
    grid on
    axis([1 size(movementData.moveDist,1) floor(min(medfilt1([movementData.moveDist(:,1);movementData.moveDist(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.moveDist(:,1);movementData.moveDist(:,2)],medFiltSize)))])
    subplot(4,1,2)
    plot(1:size(movementData.moveDist,1),medfilt1(movementData.moveDist(:,2),medFiltSize),'b')
    xlabel('Frame')
    ylabel('Y Movement (\mum)')
    grid on
    axis([1 size(movementData.moveDist,1) floor(min(medfilt1([movementData.moveDist(:,1);movementData.moveDist(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.moveDist(:,1);movementData.moveDist(:,2)],medFiltSize)))])
    subplot(4,1,3)
    plot(movementData.ballData(:,1),movementData.ballData(:,2),'k')
    title('\fontsize{20pt}\bf{Ball Movement}')
    xlabel('Time (s)')
    ylabel('Movement')
    grid on
    axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) -1 ceil(max(movementData.ballData(:,2)))])
    subplot(4,1,4)
    if all(movementData.emgData(:,2) == 0)
        title('\fontsize{20pt}\bf{No EMG Data0}')
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
    ylabel('Movement')
    grid on
    axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) -1 ceil(max(movementData.ballData(:,2)))])
    
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
    ylabel('Movement')
    grid on
    axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) -1 ceil(max(movementData.ballData(:,2)))])
    set(x3,'Position',[.05, .06, .9, .23])
    
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
        timeVecX = linspace(round(movementData.motionEvents(1,1)-movementData.motionEvents(1,2)),round(movementData.motionEvents(1,3)-movementData.motionEvents(1,3)),length(meanX));
        timeVecY = linspace(round(movementData.motionEvents(1,1)-movementData.motionEvents(1,2)),round(movementData.motionEvents(1,3)-movementData.motionEvents(1,3)),length(meanY));
        
        subplot(4,1,1)
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
        title(['\fontsize{20pt}\bf{X Position - Locomotion Events}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
        xlabel('Time (s)')
        ylabel('X Position (\mum)')
        grid on
        axis([1 size(movementData.targetPosition,1) floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)))])
        subplot(4,1,2)
        plot(timeVecX,meanX,'k')
        hold on
        f = fill([meanX flip(meanX)],cIntFillPtsX,'r','Linestyle','none');
        set(f,'facea',[.2]);
        hold off
        xlabel('Time (s)')
        ylabel('X Position (\mum)')
        
        subplot(4,1,3)
        plot((1:size(movementData.targetPosition,2))*movementData.secondsPerFrame,medfilt1(movementData.targetPosition(:,2),medFiltSize),'b')
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
        title(['\fontsize{20pt}\bf{Y Position - Locomotion Events}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
        xlabel('Time (s)')
        ylabel('Y Position (\mum)')
        grid on
        axis([1 size(movementData.targetPosition,1) floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)))])
        subplot(4,1,4)
        plot(timeVecY,meanY,'k')
        hold on
        f = fill([meanY flip(meanY)],cIntFillPtsY,'r','Linestyle','none');
        set(f,'facea',[.2]);
        hold off
        xlabel('Time (s)')
        ylabel('Y Position (\mum)')
    end
    
    h(9) = figure('Color','White');
    EMGEventsLocationsY = [];
    EMGEventsLocationsY = [];
    if size(movementData.EMGEvents,1) == 0
        title('No EMG Events To Plot')
    else
        for n = 1:size(movementData.EMGEvents,1)
            motionVectorX = movementData.targetPosition(movementData.EMGEvents(n,4):movementData.EMGEvents(n,6),1);
            motionVectorY = movementData.targetPosition(movementData.EMGEvents(n,4):movementData.EMGEvents(n,6),2);
            if n > 1
                if length(motionVectorX) > size(EMGEventsLocationsY,2)
                    motionVectorX = motionVectorX(1:size(EMGEventsLocationsY,2));
                elseif length(motionVectorX) < size(EMGEventsLocationsY,2)
                    EMGEventsLocationsY = EMGEventsLocationsY(:,1:length(motionVectorX));
                end
                if length(motionVectorY) > size(EMGEventsLocationsY,2)
                    motionVectorY = motionVectorY(1:size(EMGEventsLocationsY,2));
                elseif length(motionVectorY) < size(EMGEventsLocationsY,2)
                    EMGEventsLocationsY = EMGEventsLocationsY(:,1:length(motionVectorY));
                end
            end
            EMGEventsLocationsY(end+1,:) = medfilt1(motionVectorX-motionVectorX(1),medFiltSize);
            EMGEventsLocationsY(end+1,:) = medfilt1(motionVectorY-motionVectorY(1),medFiltSize);
        end
        [meanX,cIntFillPtsX] = getCIntMeanAndFillPts(EMGEventsLocationsY,90);
        [meanY,cIntFillPtsY] = getCIntMeanAndFillPts(EMGEventsLocationsY,90);
        timeVecX = linspace(round(movementData.EMGEvents(1,1)-movementData.EMGEvents(1,2)),round(movementData.EMGEvents(1,3)-movementData.EMGEvents(1,3)),length(meanX));
        timeVecY = linspace(round(movementData.EMGEvents(1,1)-movementData.EMGEvents(1,2)),round(movementData.EMGEvents(1,3)-movementData.EMGEvents(1,3)),length(meanY));
        
        subplot(4,1,1)
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
        title(['\fontsize{20pt}\bf{X Position - EMG Events}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
        xlabel('Time (s)')
        ylabel('X Position (\mum)')
        grid on
        axis([1 size(movementData.targetPosition,1) floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)))])
        subplot(4,1,2)
        plot(timeVecX,meanX,'k')
        hold on
        f = fill([meanX flip(meanX)],cIntFillPtsX,'r','Linestyle','none');
        set(f,'facea',[.2]);
        hold off
        xlabel('Time (s)')
        ylabel('X Position (\mum)')
        
        subplot(4,1,3)
        plot((1:size(movementData.targetPosition,2))*movementData.secondsPerFrame,medfilt1(movementData.targetPosition(:,2),medFiltSize),'b')
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
        title(['\fontsize{20pt}\bf{Y Position - EMG Events}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
        xlabel('Time (s)')
        ylabel('Y Position (\mum)')
        grid on
        axis([1 size(movementData.targetPosition,1) floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)))])
        subplot(4,1,4)
        plot(timeVecY,meanY,'k')
        hold on
        f = fill([meanY flip(meanY)],cIntFillPtsY,'r','Linestyle','none');
        set(f,'facea',[.2]);
        hold off
        xlabel('Time (s)')
        ylabel('Y Position (\mum)')
    end
    
    h(10) = figure('Color','White');
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
        timeVecX = linspace(round(movementData.EMGNoMotionEvents(1,1)-movementData.EMGNoMotionEvents(1,2)),round(movementData.EMGNoMotionEvents(1,3)-movementData.EMGNoMotionEvents(1,3)),length(meanX));
        timeVecY = linspace(round(movementData.EMGNoMotionEvents(1,1)-movementData.EMGNoMotionEvents(1,2)),round(movementData.EMGNoMotionEvents(1,3)-movementData.EMGNoMotionEvents(1,3)),length(meanY));
        
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
        title(['\fontsize{20pt}\bf{X Position - EMG (No Locomotion) Events}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
        xlabel('Time (s)')
        ylabel('X Position (\mum)')
        grid on
        axis([1 size(movementData.targetPosition,1) floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)))])
        subplot(4,1,2)
        plot(timeVecX,meanX,'k')
        hold on
        f = fill([meanX flip(meanX)],cIntFillPtsX,'r','Linestyle','none');
        set(f,'facea',[.2]);
        hold off
        xlabel('Time (s)')
        ylabel('X Position (\mum)')
        
        subplot(4,1,3)
        plot((1:size(movementData.targetPosition,2))*movementData.secondsPerFrame,medfilt1(movementData.targetPosition(:,2),medFiltSize),'b')
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
        title(['\fontsize{20pt}\bf{Y Position - EMG (No Locomotion) Events}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
        xlabel('Time (s)')
        ylabel('Y Position (\mum)')
        grid on
        axis([1 size(movementData.targetPosition,1) floor(min(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize))) ceil(max(medfilt1([movementData.targetPosition(:,1);movementData.targetPosition(:,2)],medFiltSize)))])
        subplot(4,1,4)
        plot(timeVecY,meanY,'k')
        hold on
        f = fill([meanY flip(meanY)],cIntFillPtsY,'r','Linestyle','none');
        set(f,'facea',[.2]);
        hold off
        xlabel('Time (s)')
        ylabel('Y Position (\mum)')
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
    savefig(h,[matFileName1(1:end-21) '_outputPlots_Layer' num2str(movementData.layer) '.fig']);
else
    load(matFileName1);
    stationaryData = movementData;
    load(matFileName2);
    medFiltSize = 6;
    
    subtitle = [num2str(1/movementData.secondsPerFrame) ' Frames/s, ' num2str(movementData.secondsPerFrame*(diff(movementData.frames)+1)) ' Seconds, ' num2str(movementData.objMag*movementData.digMag) 'x Magnification (' num2str(movementData.objMag) 'x Objective, ' num2str(movementData.digMag) 'x Digital), Turnabout = ' num2str(movementData.turnabout)];
    ballDataNoAbs = convertBallVoltToMPS(movementData.ballData(:,2));
    movemeplontData.ballData(:,2) = abs(convertBallVoltToMPS(movementData.ballData(:,2)));
    
    h(1) = figure('Color','White');
    x1 = subplot(4,1,1);
    plot([1:size(movementData.targetPosition,1)]*2*movementData.secondsPerFrame,medfilt1(movementData.targetPosition(:,1)-movementData.targetPosition(1,1),medFiltSize),'r')
    hold on
    plot([1:size(stationaryData.targetPosition,1)]*2*movementData.secondsPerFrame,medfilt1(stationaryData.targetPosition(:,1)-stationaryData.targetPosition(1,1),medFiltSize),'b')
    legend('Brain','Skull')
    title(['\fontsize{20pt}\bf{Object Position per Frame}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Time (s)')
    ylabel('X Position (\mum)')
    grid on
    maxVal = ceil(max(abs(medfilt1([movementData.targetPosition(:,1)-movementData.targetPosition(1,1);movementData.targetPosition(:,2)-movementData.targetPosition(1,2);stationaryData.targetPosition(:,1)-stationaryData.targetPosition(1,1);stationaryData.targetPosition(:,2)-stationaryData.targetPosition(1,2)],medFiltSize))));
    axis([0 size(movementData.targetPosition,1)*2*movementData.secondsPerFrame -maxVal maxVal])
%     set(x1,'Position',[.05, .68, .9, .23])
    if movementData.hemisphere == 1
        text(0,maxVal,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
        text(0,-maxVal,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    else
        text(0,maxVal,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
        text(0,-maxVal,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    end
    x2 = subplot(4,1,2);
    plot([1:size(movementData.targetPosition,1)]*2*movementData.secondsPerFrame,-1*medfilt1(movementData.targetPosition(:,2)-movementData.targetPosition(1,2),medFiltSize),'r')
    hold on
    plot([1:size(stationaryData.targetPosition,1)]*2*movementData.secondsPerFrame,-1*medfilt1(stationaryData.targetPosition(:,2)-stationaryData.targetPosition(1,2),medFiltSize),'b')
    legend('Brain','Skull')
    xlabel('Time (s)')
    ylabel('Y Position (\mum)')
    grid on
    axis([0 size(movementData.targetPosition,1)*2*movementData.secondsPerFrame -maxVal maxVal])
%     set(x2,'Position',[.05, .39, .9, .23])
    text(0,maxVal,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(0,-maxVal,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    x3 = subplot(4,1,3);
    plot(movementData.ballData(:,1),movementData.ballData(:,2),'k')
    title('\fontsize{20pt}\bf{Ball Movement}')
    xlabel('Time (s)')
    ylabel('m/s')
    grid on
    axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) 0 ceil(max(movementData.ballData(:,2)*10))/10])
%     set(x3,'Position',[.05, .06, .9, .23])
    x4 = subplot(4,1,4);
    semilogy(movementData.emgData(:,1),movementData.emgData(:,2),'k')
    title('\fontsize{20pt}\bf{EMG}')
    xlabel('Time (s)')
    ylabel('Amplitude (a.u.)')
    grid on
    axis([min(movementData.emgData(:,1)) max(movementData.emgData(:,1)) 0 ceil(max(movementData.emgData(:,2)))])
    
    h(2) = figure('Color','White');
    posL1 = [medfilt1(movementData.targetPosition(:,1) - movementData.targetPosition(1,1),medFiltSize), medfilt1(movementData.targetPosition(:,2) - movementData.targetPosition(1,2),medFiltSize)];
    meanPosValL1 = [.5*(posL1(1:end-1,1) + posL1(2:end,1)),.5*(posL1(1:end-1,2) + posL1(2:end,2))];
    posL1 = zipperVecs(posL1,meanPosValL1);
%     posL1 = posL1(2:end,:);
    posL2 = [medfilt1(stationaryData.targetPosition(:,1) - stationaryData.targetPosition(1,1),medFiltSize), medfilt1(stationaryData.targetPosition(:,2) - stationaryData.targetPosition(1,2),medFiltSize)];
    meanPosValL2 = [.5*(posL2(1:end-1,1) + posL2(2:end,1)),.5*(posL2(1:end-1,2) + posL2(2:end,2))];
    posL2 = zipperVecs(posL2,meanPosValL2);
    if size(posL1,1) > size(posL2,1)
        posL1 = posL1(1:end-1,:);
        posL2 = [posL2(1,:); posL2(:,:)];
    elseif size(posL1,1) < size(posL2,1)
        posL2 = posL2(1:end-1,:);
        posL1 = [posL1(1,:); posL1(:,:)];
    end
%     posDiff = [posL2(:,1)-posL1(:,1),posL2(:,2)-posL1(:,2)];
    medFiltData = [posL1(:,1)-posL2(:,1),-1*(posL1(:,2)-posL2(:,2))];
    maxPosDiff = ceil(max(max(medFiltData)));
    maxMoveData = ceil(max([posL1(:,1);posL1(:,2);posL2(:,1);posL2(:,2);]));
    x1 = subplot(4,1,1);
    plot([1:size(medFiltData,1)]*movementData.secondsPerFrame,medFiltData(:,1),'k')
    title(['\fontsize{20pt}\bf{Position of Brain in Skull}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Time (s)')
    ylabel('X Position (\mum)')
    grid on
    axis([0 size(movementData.targetPosition,1)*2*movementData.secondsPerFrame -maxPosDiff maxPosDiff])
%     set(x1,'Position',[.05, .68, .9, .23])
    if movementData.hemisphere == 1
        text(0,maxVal,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
        text(0,-maxVal,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    else
        text(0,maxVal,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
        text(0,-maxVal,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    end
    x2 = subplot(4,1,2);
    plot([1:size(medFiltData,1)]*movementData.secondsPerFrame,medFiltData(:,2),'k')
    xlabel('Time (s)')
    ylabel('Y Position (\mum)')
    grid on
    axis([0 size(movementData.targetPosition,1)*2*movementData.secondsPerFrame -maxPosDiff maxPosDiff])
%     set(x2,'Position',[.05, .39, .9, .23])
    text(0,maxVal,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(0,-maxVal,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    x3 = subplot(4,1,3);
    plot(movementData.ballData(:,1),movementData.ballData(:,2),'k')
    title('\fontsize{20pt}\bf{Ball Movement}')
    xlabel('Time (s)')
    ylabel('m/s')
    grid on
    axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) 0 ceil(max(movementData.ballData(:,2)*10))/10])
%     set(x3,'Position',[.05, .06, .9, .23])
    x4 = subplot(4,1,4);
    semilogy(movementData.emgData(:,1),movementData.emgData(:,2),'k')
    title('\fontsize{20pt}\bf{EMG}')
    xlabel('Time (s)')
    ylabel('Amplitude (a.u.)')
    grid on
    axis([min(movementData.emgData(:,1)) max(movementData.emgData(:,1)) 0 ceil(max(movementData.emgData(:,2)))])
    
    h(3) = figure('Color','White');   
    k = convhull(medFiltData(:,1),medFiltData(:,2));
    plot(medFiltData(k,1),medFiltData(k,2),'b',medFiltData(:,1),medFiltData(:,2),'k');
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
    
    h(4) = figure('Color','White');
    hist3(medFiltData,'CdataMode','auto','Nbins',[15 15]);
    colorbar;
    view(2);
    axis equal square
    axis([-maxVal maxVal -maxVal maxVal])
    ax = gca;
    set(ax,'ColorScale','log')
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    title(['\fontsize{20pt}\bf{Position of Brain in Skull}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('\mum')
    ylabel('\mum')
    
    h(5) = figure('Color','White');
    histogram(medFiltData(:,1),50);
    title(['\fontsize{20pt}\bf{Position of Target Object Histogram (X)}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('\mum')
    ylabel('Number of Data Points')
%     set(gca, 'YScale', 'log')
    xlimVal = xlim;
    ylimVal = ylim;
    if movementData.hemisphere == 1
        text(xlimVal(1),ylimVal(2),'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
        text(xlimVal(2),ylimVal(2),'Lateral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    else
        text(xlimVal(1),ylimVal(2),'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
        text(xlimVal(2),ylimVal(2),'Medial','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    end
    
    
    h(6) = figure('Color','White');
    histogram(medFiltData(:,2),50);
    title(['\fontsize{20pt}\bf{Position of Target Object Histogram (Y)}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('\mum')
    ylabel('Number of Data Points')
%     set(gca, 'YScale', 'log')
    xlimVal = xlim;
    ylimVal = ylim;
    text(xlimVal(1),ylimVal(2),'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    text(xlimVal(2),ylimVal(2),'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
    
    motionVec = pcaMotionAnalysis(medFiltData);
    h(7) = figure('Color','White');
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
    
    [PowerX,HzX,ErrorX,PowerY,HzY,ErrorY] = motionSpectrumAnalysis(medFiltData(1:1850,:));
    h(8) = figure('Color','White');
    semilogy(HzX,PowerX,'k')
    hold on
    f = fill([HzX flip(HzX)],[ErrorX(1,:) flip(ErrorX(2,:))],'r','Linestyle','none');
    set(f,'facea',[.2]);
    hold off
    title(['\fontsize{20pt}\bf{X Position Frequency Domain}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Frequency (Hz)')
    ylabel('Power')
    
    h(9) = figure('Color','White');
    semilogy(HzY,PowerY,'k')
    hold on
    f = fill([HzY flip(HzY)],[ErrorY(1,:) flip(ErrorY(2,:))],'r','Linestyle','none');
    set(f,'facea',[.2]);
    hold off
    title(['\fontsize{20pt}\bf{Y Position Frequency Domain}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Frequency (Hz)')
    ylabel('Power')
    
    h(10) = figure('Color','White');
    movementFrameSec = [1:size(medFiltData,1)]*movementData.secondsPerFrame;
    secondsPerSample = movementData.ballData(2,1) - movementData.ballData(1,1);
    samplesPerSecond = 1/secondsPerSample;
    acc = diff(medfilt1(ballDataNoAbs,15))/secondsPerSample;
    accEventSec = [];
    xPosBefore = [];
    yPosBefore = [];
    xPosAfter = [];
    yPosAfter = [];
    xMoveBefore = [];
    yMoveBefore = [];
    xMoveAfter = [];
    yMoveAfter = [];
    accSampleNum = ceil(samplesPerSecond);
    while accSampleNum <= length(acc)-samplesPerSecond
        if acc(accSampleNum) >= 15
            accSec = accSampleNum*secondsPerSample;
            accEventSec(end+1,1) = movementFrameSec(abs(movementFrameSec - accSec) == min(abs(movementFrameSec - accSec)));
            secBeforeInd = ((accEventSec(end) - 1) <= movementFrameSec) & (movementFrameSec <= accEventSec(end));
            secAfterInd = (accEventSec(end) <= movementFrameSec) & (movementFrameSec <= (accEventSec(end) + 1));
            xPosBefore(end+1,1) = mean(medFiltData(secBeforeInd,1));
            yPosBefore(end+1,1) = mean(medFiltData(secBeforeInd,2));
            xPosAfter(end+1,1) = mean(medFiltData(secAfterInd,1));
            yPosAfter(end+1,1) = mean(medFiltData(secAfterInd,2));
            xMoveBefore(end+1,1) = max(medFiltData(secBeforeInd,1)) -  min(medFiltData(secBeforeInd,1));
            yMoveBefore(end+1,1) = max(medFiltData(secBeforeInd,2)) -  min(medFiltData(secBeforeInd,2));
            xMoveAfter(end+1,1) = max(medFiltData(secAfterInd,1)) -  min(medFiltData(secAfterInd,1));
            yMoveAfter(end+1,1) = max(medFiltData(secAfterInd,2)) -  min(medFiltData(secAfterInd,2));
            accSampleNum = accSampleNum + ceil(samplesPerSecond);
        else
            accSampleNum = accSampleNum + 1;
        end
    end
    x1 = subplot(3,1,1);
    plot([1:size(medFiltData,1)]*movementData.secondsPerFrame,medFiltData(:,1),'k')
    title(['\fontsize{20pt}\bf{Position of Brain in Skull}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Time (s)')
    ylabel('X Position (\mum)')
    grid on
    axis([0 size(movementData.targetPosition,1)*2*movementData.secondsPerFrame -maxPosDiff maxPosDiff])
    set(x1,'Position',[.05, .68, .9, .23])
    if movementData.hemisphere == 1
        text(0,maxVal,'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
        text(0,-maxVal,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    else
        text(0,maxVal,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
        text(0,-maxVal,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    end
    hold on
    for n = 1:length(accEventSec)
        plot([accEventSec(n) accEventSec(n)],[-maxPosDiff maxPosDiff],'k--')
        f = fill([[accEventSec(n)-1 accEventSec(n)] flip([accEventSec(n)-1 accEventSec(n)])],[[maxPosDiff maxPosDiff] flip([-maxPosDiff -maxPosDiff])],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([[accEventSec(n) accEventSec(n)+1] flip([accEventSec(n) accEventSec(n)+1])],[[maxPosDiff maxPosDiff] flip([-maxPosDiff -maxPosDiff])],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
    x2 = subplot(3,1,2);
    plot([1:size(medFiltData,1)]*movementData.secondsPerFrame,medFiltData(:,2),'k')
    xlabel('Time (s)')
    ylabel('Y Position (\mum)')
    grid on
    axis([0 size(movementData.targetPosition,1)*2*movementData.secondsPerFrame -maxPosDiff maxPosDiff])
    set(x2,'Position',[.05, .39, .9, .23])
    text(0,maxVal,'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(0,-maxVal,'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    hold on
    for n = 1:length(accEventSec)
        plot([accEventSec(n) accEventSec(n)],[-maxPosDiff maxPosDiff],'k--')
        f = fill([[accEventSec(n)-1 accEventSec(n)] flip([accEventSec(n)-1 accEventSec(n)])],[[maxPosDiff maxPosDiff] flip([-maxPosDiff -maxPosDiff])],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([[accEventSec(n) accEventSec(n)+1] flip([accEventSec(n) accEventSec(n)+1])],[[maxPosDiff maxPosDiff] flip([-maxPosDiff -maxPosDiff])],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
    x3 = subplot(3,1,3);
    plot([1:length(acc)]*secondsPerSample,acc,'k')
    title('\fontsize{20pt}\bf{Ball Movement}')
    xlabel('Time (s)')
    ylabel('m/s^2')
    grid on
    axis([0 length(acc)*secondsPerSample -ceil(max(abs(acc))/10)*10 ceil(max(abs(acc))/10)*10])
    set(x3,'Position',[.05, .06, .9, .23])
    
    h(11) = figure('Color','White');
    x1 = subplot(2,2,1);
    plotSpread([xPosBefore xPosAfter],'xNames',{'Before','After'})
    title('\fontsize{20pt}\bf{Brain Mean X Position 1s Before and After Acceleration Spike}')
    ylabel('Brain Position From Baseline (\mum)')
    if movementData.hemisphere == 1
        text(0,ceil(max([xPosBefore;xPosAfter])),'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
        text(0,floor(min([xPosBefore;xPosAfter])),'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    else
        text(0,ceil(max([xPosBefore;xPosAfter])),'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
        text(0,floor(min([xPosBefore;xPosAfter])),'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    end
    x2 = subplot(2,2,2);
    plotSpread([yPosBefore yPosAfter],'xNames',{'Before','After'})
    title('\fontsize{20pt}\bf{Brain Mean Y Position 1s Before and After Acceleration Spike}')
    ylabel('Brain Position From Baseline (\mum)')
    text(0,ceil(max([yPosBefore;yPosAfter])),'Rostral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    text(0,floor(min([yPosBefore;yPosAfter])),'Caudal','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    x3 = subplot(2,2,3);
    plotSpread([xMoveBefore xMoveAfter],'xNames',{'Before','After'})
    title('\fontsize{20pt}\bf{Brain X Movement 1s Before and After Acceleration Spike}')
    ylabel('Brain Displacement (\mum)')
    if movementData.hemisphere == 1
        text(0,ceil(max([xMoveBefore;xMoveAfter])),'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
        text(0,floor(min([xMoveBefore;xMoveAfter])),'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    else
        text(0,ceil(max([xMoveBefore;xMoveAfter])),'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
        text(0,floor(min([xMoveBefore;xMoveAfter])),'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    end
    x4 = subplot(2,2,4);
    plotSpread([yMoveBefore yMoveAfter],'xNames',{'Before','After'})
    title('\fontsize{20pt}\bf{Brain Y Movement 1s Before and After Acceleration Spike}')
    ylabel('Brain Displacement (\mum)')
    text(0,ceil(max([yMoveBefore;yMoveAfter])),'Rostral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    text(0,floor(min([yMoveBefore;yMoveAfter])),'Caudal','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    
%     h(7) = figure('Color','White');
%     title(['\fontsize{20pt}\bf{Position of Brain in Skull}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
%     xlabel('\mum')
%     ylabel('\mum')
%     if movementData.hemisphere == 1
%         text(maxMoveData,0,'Lateral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
%         text(-maxMoveData,0,'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
%         text(0,maxMoveData,'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
%         text(0,-maxMoveData,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
%     else
%         text(maxMoveData,0,'Medial','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
%         text(-maxMoveData,0,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
%         text(0,maxMoveData,'Rostral','VerticalAlignment','top','HorizontalAlignment','right','FontSize',15);
%         text(0,-maxMoveData,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15);
%     end
%     axis equal square
%     for n = 1:length(medFiltData(:,1))-1
% %         quiver(medFiltData(n,1),medFiltData(n,2),medFiltData(n+1,1)-medFiltData(n,1),medFiltData(n+1,2)-medFiltData(n,2),0)
%         drawArrow([medFiltData(n,1);medFiltData(n,2)],[medFiltData(n+1,1);medFiltData(n+1,2)]);
% %         hText = text(-maxMoveData+.2,-maxMoveData+.2,['Brain Displacement = ' num2str(sqrt((medFiltData(n+1,1))^2+(medFiltData(n+1,2))^2))]);
%         axis([-maxMoveData maxMoveData -maxMoveData maxMoveData])
%         ax = gca;
%         ax.XAxisLocation = 'origin';
%         ax.YAxisLocation = 'origin';
% %         pause(.001)
% %         delete(hText);
%     end
    
    
    
    % Save figures to single .fig file
    savefig(h,[matFileName1(1:end-21) '_2layerComparisonOutputPlots.fig']);
end
end

function zipPos = zipperVecs(rawPos,meanPos)
zipPos = [];
for n = 1:size(rawPos,1)-1
    zipPos(end+1,:) = [rawPos(n,1),rawPos(n,2)];
    zipPos(end+1,:) = [meanPos(n,1),meanPos(n,2)];
end
zipPos(end+1,:) = [rawPos(end,1),rawPos(end,2)];
end