dateCell = {'211021','211022','211101','211102','211105','211109','211112','211116','211117','211119','211203','211216','220203','220209','220210','220211','220214','220221','220223','220303','220308','220309','220314','220318','220404','220406','220407','220429','220509','220511','220712','220714','220718','220719','220808','220809','220813','220815','220816','220822','220823','221205','221207','221208','221213'};
% load('LTADataCell.mat')
% EMGDataCellSaveNew = EMGDataCell;
EMGDataCell = {};
close all
for n = 12:size(dateCell,2)
    folderName = ['I:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/'];
    disp(folderName)
    fileList = dir(folderName);
    fileNamesCell = struct2cell(fileList);
    fileNames = fileNamesCell(1,:);
    maxRun = 0;
    for i = 1:size(fileNames,2)
        if contains(fileNames{i},dateCell{n}) && str2double(fileNames{i}(8:10)) > maxRun
             maxRun = str2double(fileNames{i}(8:10));
        end
    end
    for i = 1:maxRun
        close all
        if i > 9
            runNumberStr = num2str(i);
        else
            runNumberStr = ['0' num2str(i)];
        end
        fileName = ['I:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processe_2layerBrainInSkullDataFinal.mat'];
        if exist(fileName,'file')
            load(fileName)
            [EMGEvents,stopEMGEvents,~] = detectEMGEvents(movementData.emgData,10000,movementData.secondsPerFrame);
            EMGDataCell = plotEMGEventsSub(EMGEvents,stopEMGEvents,movementData,fileName,EMGDataCell);
            goodRun = input('');
            while goodRun > 0
                if goodRun == 1
                    removePointsX = input('');
                    EMGEvents(removePointsX,:) = [];
                elseif goodRun == 2
                    removePointsY = input('');
                    stopEMGEvents(removePointsY,:) = [];
                elseif goodRun == 3
                    [EMGEvents,stopEMGEvents,~] = detectEMGEvents(movementData.emgData,10000,movementData.secondsPerFrame);
                elseif goodRun == 4
                    EMGEvents = [];
                    stopEMGEvents = [];
                end
                close all
                EMGDataCell = plotEMGEventsSub(EMGEvents,stopEMGEvents,movementData,fileName,EMGDataCell);
                goodRun = input('');
            end
        end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subfunctions
function [EMGEvents,stopEMGEvents,EMGVelocityThresh] = detectEMGEvents(procemgData,analogSampleRate,secondsPerFrame)

EMGEvents = [];
stopEMGEvents = [];
EMGVelocityThresh = [];
analogSampleRate = 1/(procemgData(2,1) - procemgData(1,1));
if size(procemgData,1) <= 7.1*analogSampleRate
    return
else
    procemgData = abs(procemgData);
    figure('units','normalized','outerposition',[0 0 1 1])
    plot(procemgData(:,1),procemgData(:,2));
    title('Select upper limit of ball EMG events')
    [~,EMGVelocityThresh] = ginput(1);
    close;
    startFrame = 2*analogSampleRate + 1;
    stopFrame = size(procemgData,1) - (6*analogSampleRate) - 1;
    i = find(procemgData(:,2)>EMGVelocityThresh);
    for n = 1:length(i)
        index = i(n);
        iStart = index - (2*analogSampleRate);
        iStop = index + (6*analogSampleRate);
        if index > startFrame && index < stopFrame && ~any(procemgData(iStart:index-1,2)>EMGVelocityThresh) && (sum(procemgData(index+1:iStop,2)>EMGVelocityThresh)/(iStop-index))>0.5
            EMGEvents(end+1,1:3) = [iStart/analogSampleRate index/analogSampleRate iStop/analogSampleRate];
            EMGEvents(end,4:6) = EMGEvents(end,1:3)./secondsPerFrame;
            if EMGEvents(end,5) - round(EMGEvents(end,5)) > 0
                EMGEvents(end,4:6) = floor(EMGEvents(end,4:6));
            elseif EMGEvents(end,5) - round(EMGEvents(end,5)) < 0
                EMGEvents(end,4:6) = ceil(EMGEvents(end,4:6));
            end
        end
    end
    i = find(procemgData(:,2)<EMGVelocityThresh);
    for n = 1:length(i)
        index = i(n);
        iStart = index - (2*analogSampleRate);
        iStop = index + (6*analogSampleRate);
        if index > startFrame && index < stopFrame && procemgData(index-1,2)>EMGVelocityThresh && all(procemgData(index+1:iStop,2)<EMGVelocityThresh) && (sum(procemgData(iStart:index-1,2)>EMGVelocityThresh)/(index-iStart))>0.3
            stopEMGEvents(end+1,1:3) = [iStart/analogSampleRate index/analogSampleRate iStop/analogSampleRate];
            stopEMGEvents(end,4:6) = stopEMGEvents(end,1:3)./secondsPerFrame;
            if stopEMGEvents(end,5) - round(stopEMGEvents(end,5)) > 0
                stopEMGEvents(end,4:6) = floor(stopEMGEvents(end,4:6));
            elseif stopEMGEvents(end,5) - round(stopEMGEvents(end,5)) < 0
                stopEMGEvents(end,4:6) = ceil(stopEMGEvents(end,4:6));
            end
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dataCell = plotEMGEventsSub(EMGEvents,stopEMGEvents,movementData,fileName,dataCell)
if isempty(dataCell) || ~any(strcmp(fileName,dataCell(:,1)))
    dataCell{end+1,1} = fileName;
    dataCellInd = size(dataCell,1);
else
    dataCellInd = find(strcmp(fileName,dataCell(:,1)));
end
h(1) = figure('Color','White');figure(1)
set(gcf,'units','normalized','outerposition',[0 0.5 .5 .5])
EMGEventsLocationsX = [];
EMGEventsLocationsY = [];
if size(EMGEvents,1) == 0
    title('No Ball EMG Events To Plot')
    dataCell(dataCellInd,2:4) = {EMGEvents,NaN,NaN};
else
    for n = 1:size(EMGEvents,1)
        EMGVectorX = movementData.targetPosition(EMGEvents(n,4):EMGEvents(n,6),1);
        EMGVectorY = movementData.targetPosition(EMGEvents(n,4):EMGEvents(n,6),2)*-1;
        if n > 1
            if length(EMGVectorX) > size(EMGEventsLocationsX,2)
                EMGVectorX = EMGVectorX(1:size(EMGEventsLocationsX,2));
            elseif length(EMGVectorX) < size(EMGEventsLocationsX,2)
                EMGEventsLocationsX = EMGEventsLocationsX(:,1:length(EMGVectorX));
            end
            if length(EMGVectorY) > size(EMGEventsLocationsY,2)
                EMGVectorY = EMGVectorY(1:size(EMGEventsLocationsY,2));
            elseif length(EMGVectorY) < size(EMGEventsLocationsY,2)
                EMGEventsLocationsY = EMGEventsLocationsY(:,1:length(EMGVectorY));
            end
        end
        EMGEventsLocationsX(end+1,:) = medfilt1(EMGVectorX-EMGVectorX(1),6);
        EMGEventsLocationsY(end+1,:) = medfilt1(EMGVectorY-EMGVectorY(1),6);
%         EMGEventsLocationsX(end+1,:) = medfilt1(EMGVectorX,6);
%         EMGEventsLocationsY(end+1,:) = medfilt1(EMGVectorY,6);
    end
    [meanX,cIntFillPtsX] = getCIntMeanAndFillPts(EMGEventsLocationsX,90);
    if movementData.hemisphere == 2
        meanX = meanX * -1;
    end
    [meanY,cIntFillPtsY] = getCIntMeanAndFillPts(EMGEventsLocationsY,90);
    timeVecX = linspace(round(EMGEvents(1,1)-EMGEvents(1,2)),round(EMGEvents(1,3)-EMGEvents(1,2)),length(meanX));
    timeVecY = linspace(round(EMGEvents(1,1)-EMGEvents(1,2)),round(EMGEvents(1,3)-EMGEvents(1,2)),length(meanY));
    subplot(3,1,1)
    plot((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame,movementData.targetPosition(:,1),'b')
    moveMin = floor(min([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1]));
    moveMax = ceil(max([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1]));
    hold on
    for n = 1:size(EMGEventsLocationsX,1)
        plot([EMGEvents(n,2) EMGEvents(n,2)],[moveMin moveMax],'k--')
        f = fill([EMGEvents(n,1) EMGEvents(n,2) EMGEvents(n,2) EMGEvents(n,1)],[moveMin moveMin moveMax moveMax],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([EMGEvents(n,2) EMGEvents(n,3) EMGEvents(n,3) EMGEvents(n,2)],[moveMin moveMin moveMax moveMax],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
    %     title(['\fontsize{20pt}\bf{X Position - LocoEMG Events}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Time (s)')
    ylabel('X Position (\mum)')
    title('\fontsize{20pt}\bf{EMG During LocoEMG Events}')
    grid on
    if movementData.hemisphere == 1
        text(0,ceil(max([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)*-1])),'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
        text(0,floor(min([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)*-1])),'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    else
        text(0,ceil(max([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)*-1])),'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
        text(0,floor(min([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)*-1])),'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    end
    axis([0 size(movementData.targetPosition,1)*movementData.secondsPerFrame floor(min([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1])) ceil(max([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1]))])
    subplot(3,1,2)
    title(fileName)
    plot((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame,-1*movementData.targetPosition(:,2)*-1,'b')
    moveMin = floor(min([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1]));
    moveMax = ceil(max([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1]));
    hold on
    for n = 1:size(EMGEventsLocationsX,1)
        plot([EMGEvents(n,2) EMGEvents(n,2)],[moveMin moveMax],'k--')
        f = fill([EMGEvents(n,1) EMGEvents(n,2) EMGEvents(n,2) EMGEvents(n,1)],[moveMin moveMin moveMax moveMax],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([EMGEvents(n,2) EMGEvents(n,3) EMGEvents(n,3) EMGEvents(n,2)],[moveMin moveMin moveMax moveMax],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
    %     title(['\fontsize{20pt}\bf{Y Position - LocoEMG Events}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Time (s)')
    ylabel('Y Position (\mum)')
    grid on
    text(0,-floor(min([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1])),'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(0,-ceil(max([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1])),'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    axis([0 size(movementData.targetPosition,1)*movementData.secondsPerFrame floor(min([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1])) ceil(max([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1]))])
    subplot(3,1,3)
    plot(movementData.emgData(:,1),movementData.emgData(:,2),'k')
    title('\fontsize{20pt}\bf{Ball Movement}')
    xlabel('Time (s)')
    ylabel('m/s')
    grid on
    axis([min(movementData.emgData(:,1)) max(movementData.emgData(:,1)) 0 ceil(max(movementData.emgData(:,2)*10))/10])
    hold on
    for n = 1:size(EMGEventsLocationsX,1)
        plot([EMGEvents(n,2) EMGEvents(n,2)],[moveMin moveMax],'k--')
        f = fill([EMGEvents(n,1) EMGEvents(n,2) EMGEvents(n,2) EMGEvents(n,1)],[moveMin moveMin moveMax moveMax],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([EMGEvents(n,2) EMGEvents(n,3) EMGEvents(n,3) EMGEvents(n,2)],[moveMin moveMin moveMax moveMax],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
    dataCell(dataCellInd,2:4) = {EMGEvents,[timeVecX;meanX],[timeVecY;meanY]};
end

h(2) = figure('Color','White');
set(gcf,'units','normalized','outerposition',[0 0 .5 .5])
if size(EMGEvents,1) == 0
    title('No Ball EMG Events To Plot')
else
    subplot(2,1,1)
    maxMeanVal = max(abs([meanX meanY]));
    plot(timeVecX,meanX,'k')
    hold on
    %     f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
    %     set(f,'facea',[.2]);
    plot([0 0],[-maxMeanVal maxMeanVal],'r')
    hold off
    title('\fontsize{20pt}\bf{Mean EMG During LocoEMG Events}')
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

h(3) = figure('Color','White');
set(gcf,'units','normalized','outerposition',[0.5 0.5 .5 .5])
stopEMGEventsLocationsX = [];
stopEMGEventsLocationsY = [];
if size(stopEMGEvents,1) == 0
    title('No Stop Ball EMG Events To Plot')
    dataCell(dataCellInd,5:7) = {stopEMGEvents,NaN,NaN};
else
    for n = 1:size(stopEMGEvents,1)
        EMGVectorX = movementData.targetPosition(stopEMGEvents(n,4):stopEMGEvents(n,6),1);
        EMGVectorY = movementData.targetPosition(stopEMGEvents(n,4):stopEMGEvents(n,6),2)*-1;
        if n > 1
            if length(EMGVectorX) > size(stopEMGEventsLocationsX,2)
                EMGVectorX = EMGVectorX(1:size(stopEMGEventsLocationsX,2));
            elseif length(EMGVectorX) < size(stopEMGEventsLocationsX,2)
                stopEMGEventsLocationsX = stopEMGEventsLocationsX(:,1:length(EMGVectorX));
            end
            if length(EMGVectorY) > size(stopEMGEventsLocationsY,2)
                EMGVectorY = EMGVectorY(1:size(stopEMGEventsLocationsY,2));
            elseif length(EMGVectorY) < size(stopEMGEventsLocationsY,2)
                stopEMGEventsLocationsY = stopEMGEventsLocationsY(:,1:length(EMGVectorY));
            end
        end
        stopEMGEventsLocationsX(end+1,:) = medfilt1(EMGVectorX-EMGVectorX(1),6);
        stopEMGEventsLocationsY(end+1,:) = medfilt1(EMGVectorY-EMGVectorY(1),6);
%         stopEMGEventsLocationsX(end+1,:) = medfilt1(EMGVectorX,6);
%         stopEMGEventsLocationsY(end+1,:) = medfilt1(EMGVectorY,6);
    end
    [meanX,cIntFillPtsX] = getCIntMeanAndFillPts(stopEMGEventsLocationsX,90);
    if movementData.hemisphere == 2
        meanX = meanX * -1;
    end
    [meanY,cIntFillPtsY] = getCIntMeanAndFillPts(stopEMGEventsLocationsY,90);
    timeVecX = linspace(round(stopEMGEvents(1,1)-stopEMGEvents(1,2)),round(stopEMGEvents(1,3)-stopEMGEvents(1,2)),length(meanX));
    timeVecY = linspace(round(stopEMGEvents(1,1)-stopEMGEvents(1,2)),round(stopEMGEvents(1,3)-stopEMGEvents(1,2)),length(meanY));
    
    subplot(3,1,1)
    plot((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame,movementData.targetPosition(:,1),'b')
    moveMin = floor(min([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1]));
    moveMax = ceil(max([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1]));
    hold on
    for n = 1:size(stopEMGEventsLocationsX,1)
        plot([stopEMGEvents(n,2) stopEMGEvents(n,2)],[moveMin moveMax],'k--')
        f = fill([stopEMGEvents(n,1) stopEMGEvents(n,2) stopEMGEvents(n,2) stopEMGEvents(n,1)],[moveMin moveMin moveMax moveMax],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([stopEMGEvents(n,2) stopEMGEvents(n,3) stopEMGEvents(n,3) stopEMGEvents(n,2)],[moveMin moveMin moveMax moveMax],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
    %     title(['\fontsize{20pt}\bf{X Position - LocoEMG Events}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Time (s)')
    ylabel('X Position (\mum)')
    title('\fontsize{20pt}\bf{EMG During Stop LocoEMG Events}')
    grid on
    if movementData.hemisphere == 1
        text(0,ceil(max([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)*-1])),'Lateral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
        text(0,floor(min([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)*-1])),'Medial','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    else
        text(0,ceil(max([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)*-1])),'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
        text(0,floor(min([movementData.targetPosition(:,1);-1*movementData.targetPosition(:,2)*-1])),'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    end
    axis([0 size(movementData.targetPosition,1)*movementData.secondsPerFrame floor(min([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1])) ceil(max([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1]))])
    subplot(3,1,2)
    plot((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame,-1*movementData.targetPosition(:,2)*-1,'b')
    moveMin = floor(min([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1]));
    moveMax = ceil(max([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1]));
    hold on
    for n = 1:size(stopEMGEventsLocationsX,1)
        plot([stopEMGEvents(n,2) stopEMGEvents(n,2)],[moveMin moveMax],'k--')
        f = fill([stopEMGEvents(n,1) stopEMGEvents(n,2) stopEMGEvents(n,2) stopEMGEvents(n,1)],[moveMin moveMin moveMax moveMax],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([stopEMGEvents(n,2) stopEMGEvents(n,3) stopEMGEvents(n,3) stopEMGEvents(n,2)],[moveMin moveMin moveMax moveMax],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
    %     title(['\fontsize{20pt}\bf{Y Position - LocoEMG Events}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Time (s)')
    ylabel('Y Position (\mum)')
    grid on
    text(0,-floor(min([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1])),'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(0,-ceil(max([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1])),'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    axis([0 size(movementData.targetPosition,1)*movementData.secondsPerFrame floor(min([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1])) ceil(max([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1]))])
    subplot(3,1,3)
    plot(movementData.emgData(:,1),movementData.emgData(:,2),'k')
    title('\fontsize{20pt}\bf{Ball Movement}')
    xlabel('Time (s)')
    ylabel('m/s')
    grid on
    axis([min(movementData.emgData(:,1)) max(movementData.emgData(:,1)) 0 ceil(max(movementData.emgData(:,2)*10))/10])
    hold on
    for n = 1:size(stopEMGEventsLocationsX,1)
        plot([stopEMGEvents(n,2) stopEMGEvents(n,2)],[moveMin moveMax],'k--')
        f = fill([stopEMGEvents(n,1) stopEMGEvents(n,2) stopEMGEvents(n,2) stopEMGEvents(n,1)],[moveMin moveMin moveMax moveMax],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([stopEMGEvents(n,2) stopEMGEvents(n,3) stopEMGEvents(n,3) stopEMGEvents(n,2)],[moveMin moveMin moveMax moveMax],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
    dataCell(dataCellInd,5:7) = {stopEMGEvents,[timeVecX;meanX],[timeVecY;meanY]};
end

h(4) = figure('Color','White');
set(gcf,'units','normalized','outerposition',[0.5 0 .5 .5])
if size(stopEMGEvents,1) == 0
    title('No Stop Ball EMG Events To Plot')
else
    subplot(2,1,1)
    maxMeanVal = max(abs([meanX meanY]));
    plot(timeVecX,meanX,'k')
    hold on
    %     f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
    %     set(f,'facea',[.2]);
    plot([0 0],[-maxMeanVal maxMeanVal],'r')
    hold off
    title('\fontsize{20pt}\bf{Mean EMG During Stopping LocoEMG Events}')
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
end