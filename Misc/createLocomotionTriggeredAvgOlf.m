close all
dateCell = {'220813','220815','220816'};
close all
locDataCellOlf = {};
locDataCellOlf = createLocomotionTriggeredAvgOlfSub(dateCell);
save('LTADataCellOlf.mat',locDataCellOlf)
% fileName = {};
% for n = 1:size(dateCell,2)
%     folderName = ['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/'];
%     fileList = dir(folderName);
%     fileNamesCell = struct2cell(fileList);
%     fileNames = fileNamesCell(1,:);
%     maxRun = 0;
%     for i = 1:size(fileNames,2)
%         if contains(fileNames{i},dateCell{n}) && str2double(fileNames{i}(8:10)) > maxRun
%              maxRun = str2double(fileNames{i}(8:10));
%         end
%     end
%     for i = 1:maxRun
%         if i > 9
%             runNumberStr = num2str(i);
%         else
%             runNumberStr = ['0' num2str(i)];
%         end
%         if exist(['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processe_2layerBrainInSkullDataFinal.mat'],'file')
%             fileName{end+1} = ['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processe_2layerBrainInSkullDataFinal.mat'];
%         end
%     end
% end
% medFiltSize = 6;
% timeBeforeSqueeze = 2;
% timeAfterSqueeze = 6;
% motionEventsLocationsX = [];
% motionEventsLocationsY = [];
% squeezeDataCell = {};
% for i = 1:numel(fileName)
%     load(fileName{i});
%     [squeezeMotionX,squeezeMotionY,squeezeIdxMat] = findTriggeredLocomotions(movementData,timeBeforeSqueeze,timeAfterSqueeze);
%     motionEventsLocationsX(end+1,:) = squeezeMotionX;
%     motionEventsLocationsY(end+1,:) = squeezeMotionY;
%     timeVec = linspace(timeBeforeSqueeze*-1,timeAfterSqueeze,length(squeezeMotionX));
%     squeezeDataCell(i,:) = {fileName,squeezeIdxMat,[timeVec;squeezeMotionX],[timeVec;squeezeMotionY*-1]};
%     clear movementData
% end
% save('squeezeDataCell_FS.mat','squeezeDataCell')
% load(fileName{1});
% [meanX,cIntFillPtsX] = getCIntMeanAndFillPts(motionEventsLocationsX,90);
% [meanY,cIntFillPtsY] = getCIntMeanAndFillPts(motionEventsLocationsY,90);
% timeVecX = linspace(timeBeforeSqueeze*-1,timeAfterSqueeze,length(meanX));
% timeVecY = linspace(timeBeforeSqueeze*-1,timeAfterSqueeze,length(meanY));
% 
% h(1) = figure('Color','White');
% subplot(2,1,1)
% maxMeanVal = max(abs([meanX meanY cIntFillPtsX cIntFillPtsY]));
% plot(timeVecX,meanX,'k')
% hold on
% f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
% set(f,'facea',[.2]);
% % plot([0 0],[-1 5],'r')
% % plot([2 2],[-1 5],'r')
% f = fill([0 2 2 0],[4.9 4.9 -.9 -.9],'g','Linestyle','none','FaceAlpha',0.1);
% for n = 1:size(motionEventsLocationsX,1)
%     plot(timeVecX,motionEventsLocationsX(n,:),'Color',[0,0,1,0.1])
% end
% % plot(timeVecX(brainMotionStart),meanX(brainMotionStart),'rx')
% hold off
% text(5,-1,'Medial','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(5,5,'Lateral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% title(['Figure 3f(1)' 10 '\fontsize{20pt}\bf{Mean Motion During Locomotion Events, n = ' num2str(size(motionEventsLocationsX,1)) '}'])
% xlabel('Time (s)')
% ylabel('\Delta Brian Shift (\mum)')
% ylim([-1 5])
% xlim([-2 6])
% grid on
% 
% subplot(2,1,2)
% plot(timeVecY,meanY,'k')
% hold on
% f = fill([timeVecY flip(timeVecY)],cIntFillPtsY,'r','Linestyle','none');
% set(f,'facea',[.2]);
% % plot([0 0],[-1 5],'r')
% % plot([2 2],[-1 5],'r')
% f = fill([0 2 2 0],[4.9 4.9 -.9 -.9],'g','Linestyle','none','FaceAlpha',0.1);
% for n = 1:size(motionEventsLocationsY,1)
%     plot(timeVecY,motionEventsLocationsY(n,:),'Color',[0,0,1,0.1])
% end
% % plot(timeVecY(brainMotionStart),meanY(brainMotionStart),'rx')
% hold off
% text(5,-1,'Caudal','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
% text(5,5,'Rostral','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
% xlabel('Time (s)')
% ylabel('\Delta Brian Shift (\mum)')
% ylim([-1 5])
% xlim([-2 6])
% grid on
% clear movementData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dataCell = createLocomotionTriggeredAvgOlfSub(dateCell)
% dateCell = {'211021','211022','211101','211102','211105','211109','211112','211116','211117','211119','211203','211216','220203','220209','220210','220211','220214','220221','220223','220303','220308','220309','220314','220318','220404','220406','220407','220429','220509','220511','220712','220714','220718','220719','220808','220809','220813','220815','220816','220822','220823','221205','221207','221208','221213'};
% load('LTADataCell.mat')
% dataCellSaveNew = dataCell;
dataCell = {};
close all
for n = 1:size(dateCell,2)
    folderName = ['D:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/'];
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
        fileName = ['D:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processe_2layerBrainInSkullDataFinal.mat'];
        if exist(fileName,'file')
            load(fileName)
            [motionEvents,stopMotionEvents,~] = detectMotionEvents(movementData.ballData,10000,movementData.secondsPerFrame);
            dataCell = plotMotionEventsSub(motionEvents,stopMotionEvents,movementData,fileName,dataCell);
            goodRun = input('');
            while goodRun > 0
                if goodRun == 1
                    removePointsX = input('');
                    motionEvents(removePointsX,:) = [];
                elseif goodRun == 2
                    removePointsY = input('');
                    stopMotionEvents(removePointsY,:) = [];
                elseif goodRun == 3
                    [motionEvents,stopMotionEvents,~] = detectMotionEvents(movementData.ballData,10000,movementData.secondsPerFrame);
                elseif goodRun == 4
                    motionEvents = [];
                    stopMotionEvents = [];
                end
                close all
                dataCell = plotMotionEventsSub(motionEvents,stopMotionEvents,movementData,fileName,dataCell);
                goodRun = input('');
            end
        end
    end
end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subfunctions
function [motionEvents,stopMotionEvents,motionVelocityThresh] = detectMotionEvents(procBallData,analogSampleRate,secondsPerFrame)

motionEvents = [];
stopMotionEvents = [];
motionVelocityThresh = [];
analogSampleRate = 1/(procBallData(2,1) - procBallData(1,1));
if size(procBallData,1) <= 7.1*analogSampleRate
    return
else
    procBallData = abs(procBallData);
    figure('units','normalized','outerposition',[0 0 1 1])
    plot(procBallData(:,1),procBallData(:,2));
    title('Select upper limit of ball motion events')
    [~,motionVelocityThresh] = ginput(1);
    close;
    startFrame = 2*analogSampleRate + 1;
    stopFrame = size(procBallData,1) - (3*analogSampleRate) - 1;
    i = find(procBallData(:,2)>motionVelocityThresh);
    for n = 1:length(i)
        index = i(n);
        iStart = index - (2*analogSampleRate);
        iStop = index + (3*analogSampleRate);
        if index > startFrame && index < stopFrame && ~any(procBallData(iStart:index-1,2)>motionVelocityThresh)
            motionEvents(end+1,1:3) = [iStart/analogSampleRate index/analogSampleRate iStop/analogSampleRate];
            motionEvents(end,4:6) = motionEvents(end,1:3)./secondsPerFrame;
            if motionEvents(end,5) - round(motionEvents(end,5)) > 0
                motionEvents(end,4:6) = floor(motionEvents(end,4:6));
            elseif motionEvents(end,5) - round(motionEvents(end,5)) < 0
                motionEvents(end,4:6) = ceil(motionEvents(end,4:6));
            end
        end
    end
    i = find(procBallData(:,2)<motionVelocityThresh);
    for n = 1:length(i)
        index = i(n);
        iStart = index - (2*analogSampleRate);
        iStop = index + (3*analogSampleRate);
        if index > startFrame && index < stopFrame && procBallData(index-1,2)>motionVelocityThresh && all(procBallData(index+1:iStop,2)<motionVelocityThresh)
            stopMotionEvents(end+1,1:3) = [iStart/analogSampleRate index/analogSampleRate iStop/analogSampleRate];
            stopMotionEvents(end,4:6) = stopMotionEvents(end,1:3)./secondsPerFrame;
            if stopMotionEvents(end,5) - round(stopMotionEvents(end,5)) > 0
                stopMotionEvents(end,4:6) = floor(stopMotionEvents(end,4:6));
            elseif stopMotionEvents(end,5) - round(stopMotionEvents(end,5)) < 0
                stopMotionEvents(end,4:6) = ceil(stopMotionEvents(end,4:6));
            end
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dataCell = plotMotionEventsSub(motionEvents,stopMotionEvents,movementData,fileName,dataCell)
if isempty(dataCell) || ~any(strcmp(fileName,dataCell(:,1)))
    dataCell{end+1,1} = fileName;
    dataCellInd = size(dataCell,1);
else
    dataCellInd = find(strcmp(fileName,dataCell(:,1)));
end
h(1) = figure('Color','White');figure(1)
set(gcf,'units','normalized','outerposition',[0 0.5 .5 .5])
motionEventsLocationsX = [];
motionEventsLocationsY = [];
if size(motionEvents,1) == 0
    title('No Ball Motion Events To Plot')
    dataCell(dataCellInd,2:4) = {motionEvents,NaN,NaN};
else
    for n = 1:size(motionEvents,1)
        motionVectorX = movementData.targetPosition(motionEvents(n,4):motionEvents(n,6),1);
        motionVectorY = movementData.targetPosition(motionEvents(n,4):motionEvents(n,6),2)*-1;
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
        motionEventsLocationsX(end+1,:) = medfilt1(motionVectorX-motionVectorX(1),6);
        motionEventsLocationsY(end+1,:) = medfilt1(motionVectorY-motionVectorY(1),6);
%         motionEventsLocationsX(end+1,:) = medfilt1(motionVectorX,6);
%         motionEventsLocationsY(end+1,:) = medfilt1(motionVectorY,6);
    end
    [meanX,cIntFillPtsX] = getCIntMeanAndFillPts(motionEventsLocationsX,90);
    if movementData.hemisphere == 2
        meanX = meanX * -1;
    end
    [meanY,cIntFillPtsY] = getCIntMeanAndFillPts(motionEventsLocationsY,90);
    timeVecX = linspace(round(motionEvents(1,1)-motionEvents(1,2)),round(motionEvents(1,3)-motionEvents(1,2)),length(meanX));
    timeVecY = linspace(round(motionEvents(1,1)-motionEvents(1,2)),round(motionEvents(1,3)-motionEvents(1,2)),length(meanY));
    subplot(3,1,1)
    plot((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame,movementData.targetPosition(:,1),'b')
    moveMin = floor(min([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1]));
    moveMax = ceil(max([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1]));
    hold on
    for n = 1:size(motionEventsLocationsX,1)
        plot([motionEvents(n,2) motionEvents(n,2)],[moveMin moveMax],'k--')
        f = fill([motionEvents(n,1) motionEvents(n,2) motionEvents(n,2) motionEvents(n,1)],[moveMin moveMin moveMax moveMax],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([motionEvents(n,2) motionEvents(n,3) motionEvents(n,3) motionEvents(n,2)],[moveMin moveMin moveMax moveMax],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
    %     title(['\fontsize{20pt}\bf{X Position - Locomotion Events}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Time (s)')
    ylabel('X Position (\mum)')
    title('\fontsize{20pt}\bf{Motion During Locomotion Events}')
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
    for n = 1:size(motionEventsLocationsX,1)
        plot([motionEvents(n,2) motionEvents(n,2)],[moveMin moveMax],'k--')
        f = fill([motionEvents(n,1) motionEvents(n,2) motionEvents(n,2) motionEvents(n,1)],[moveMin moveMin moveMax moveMax],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([motionEvents(n,2) motionEvents(n,3) motionEvents(n,3) motionEvents(n,2)],[moveMin moveMin moveMax moveMax],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
    %     title(['\fontsize{20pt}\bf{Y Position - Locomotion Events}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Time (s)')
    ylabel('Y Position (\mum)')
    grid on
    text(0,-floor(min([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1])),'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(0,-ceil(max([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1])),'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    axis([0 size(movementData.targetPosition,1)*movementData.secondsPerFrame floor(min([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1])) ceil(max([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1]))])
    subplot(3,1,3)
    plot(movementData.ballData(:,1),movementData.ballData(:,2),'k')
    title('\fontsize{20pt}\bf{Ball Movement}')
    xlabel('Time (s)')
    ylabel('m/s')
    grid on
    axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) 0 ceil(max(movementData.ballData(:,2)*10))/10])
    hold on
    for n = 1:size(motionEventsLocationsX,1)
        plot([motionEvents(n,2) motionEvents(n,2)],[moveMin moveMax],'k--')
        f = fill([motionEvents(n,1) motionEvents(n,2) motionEvents(n,2) motionEvents(n,1)],[moveMin moveMin moveMax moveMax],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([motionEvents(n,2) motionEvents(n,3) motionEvents(n,3) motionEvents(n,2)],[moveMin moveMin moveMax moveMax],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
    dataCell(dataCellInd,2:4) = {motionEvents,[timeVecX;meanX],[timeVecY;meanY]};
end

h(2) = figure('Color','White');
set(gcf,'units','normalized','outerposition',[0 0 .5 .5])
if size(motionEvents,1) == 0
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

h(3) = figure('Color','White');
set(gcf,'units','normalized','outerposition',[0.5 0.5 .5 .5])
stopMotionEventsLocationsX = [];
stopMotionEventsLocationsY = [];
if size(stopMotionEvents,1) == 0
    title('No Stop Ball Motion Events To Plot')
    dataCell(dataCellInd,5:7) = {stopMotionEvents,NaN,NaN};
else
    for n = 1:size(stopMotionEvents,1)
        motionVectorX = movementData.targetPosition(stopMotionEvents(n,4):stopMotionEvents(n,6),1);
        motionVectorY = movementData.targetPosition(stopMotionEvents(n,4):stopMotionEvents(n,6),2)*-1;
        if n > 1
            if length(motionVectorX) > size(stopMotionEventsLocationsX,2)
                motionVectorX = motionVectorX(1:size(stopMotionEventsLocationsX,2));
            elseif length(motionVectorX) < size(stopMotionEventsLocationsX,2)
                stopMotionEventsLocationsX = stopMotionEventsLocationsX(:,1:length(motionVectorX));
            end
            if length(motionVectorY) > size(stopMotionEventsLocationsY,2)
                motionVectorY = motionVectorY(1:size(stopMotionEventsLocationsY,2));
            elseif length(motionVectorY) < size(stopMotionEventsLocationsY,2)
                stopMotionEventsLocationsY = stopMotionEventsLocationsY(:,1:length(motionVectorY));
            end
        end
        stopMotionEventsLocationsX(end+1,:) = medfilt1(motionVectorX-motionVectorX(1),6);
        stopMotionEventsLocationsY(end+1,:) = medfilt1(motionVectorY-motionVectorY(1),6);
%         stopMotionEventsLocationsX(end+1,:) = medfilt1(motionVectorX,6);
%         stopMotionEventsLocationsY(end+1,:) = medfilt1(motionVectorY,6);
    end
    [meanX,cIntFillPtsX] = getCIntMeanAndFillPts(stopMotionEventsLocationsX,90);
    if movementData.hemisphere == 2
        meanX = meanX * -1;
    end
    [meanY,cIntFillPtsY] = getCIntMeanAndFillPts(stopMotionEventsLocationsY,90);
    timeVecX = linspace(round(stopMotionEvents(1,1)-stopMotionEvents(1,2)),round(stopMotionEvents(1,3)-stopMotionEvents(1,2)),length(meanX));
    timeVecY = linspace(round(stopMotionEvents(1,1)-stopMotionEvents(1,2)),round(stopMotionEvents(1,3)-stopMotionEvents(1,2)),length(meanY));
    
    subplot(3,1,1)
    plot((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame,movementData.targetPosition(:,1),'b')
    moveMin = floor(min([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1]));
    moveMax = ceil(max([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1]));
    hold on
    for n = 1:size(stopMotionEventsLocationsX,1)
        plot([stopMotionEvents(n,2) stopMotionEvents(n,2)],[moveMin moveMax],'k--')
        f = fill([stopMotionEvents(n,1) stopMotionEvents(n,2) stopMotionEvents(n,2) stopMotionEvents(n,1)],[moveMin moveMin moveMax moveMax],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([stopMotionEvents(n,2) stopMotionEvents(n,3) stopMotionEvents(n,3) stopMotionEvents(n,2)],[moveMin moveMin moveMax moveMax],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
    %     title(['\fontsize{20pt}\bf{X Position - Locomotion Events}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Time (s)')
    ylabel('X Position (\mum)')
    title('\fontsize{20pt}\bf{Motion During Stop Locomotion Events}')
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
    for n = 1:size(stopMotionEventsLocationsX,1)
        plot([stopMotionEvents(n,2) stopMotionEvents(n,2)],[moveMin moveMax],'k--')
        f = fill([stopMotionEvents(n,1) stopMotionEvents(n,2) stopMotionEvents(n,2) stopMotionEvents(n,1)],[moveMin moveMin moveMax moveMax],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([stopMotionEvents(n,2) stopMotionEvents(n,3) stopMotionEvents(n,3) stopMotionEvents(n,2)],[moveMin moveMin moveMax moveMax],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
    %     title(['\fontsize{20pt}\bf{Y Position - Locomotion Events}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Time (s)')
    ylabel('Y Position (\mum)')
    grid on
    text(0,-floor(min([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1])),'Rostral','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',15);
    text(0,-ceil(max([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1])),'Caudal','VerticalAlignment','top','HorizontalAlignment','left','FontSize',15);
    axis([0 size(movementData.targetPosition,1)*movementData.secondsPerFrame floor(min([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1])) ceil(max([movementData.targetPosition(:,1);movementData.targetPosition(:,2)*-1]))])
    subplot(3,1,3)
    plot(movementData.ballData(:,1),movementData.ballData(:,2),'k')
    title('\fontsize{20pt}\bf{Ball Movement}')
    xlabel('Time (s)')
    ylabel('m/s')
    grid on
    axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) 0 ceil(max(movementData.ballData(:,2)*10))/10])
    hold on
    for n = 1:size(stopMotionEventsLocationsX,1)
        plot([stopMotionEvents(n,2) stopMotionEvents(n,2)],[moveMin moveMax],'k--')
        f = fill([stopMotionEvents(n,1) stopMotionEvents(n,2) stopMotionEvents(n,2) stopMotionEvents(n,1)],[moveMin moveMin moveMax moveMax],'r','Linestyle','none');
        set(f,'facea',[.2]);
        f = fill([stopMotionEvents(n,2) stopMotionEvents(n,3) stopMotionEvents(n,3) stopMotionEvents(n,2)],[moveMin moveMin moveMax moveMax],'g','Linestyle','none');
        set(f,'facea',[.2]);
    end
    hold off
    dataCell(dataCellInd,5:7) = {stopMotionEvents,[timeVecX;meanX],[timeVecY;meanY]};
end

h(4) = figure('Color','White');
set(gcf,'units','normalized','outerposition',[0.5 0 .5 .5])
if size(stopMotionEvents,1) == 0
    title('No Stop Ball Motion Events To Plot')
else
    subplot(2,1,1)
    maxMeanVal = max(abs([meanX meanY]));
    plot(timeVecX,meanX,'k')
    hold on
    %     f = fill([timeVecX flip(timeVecX)],cIntFillPtsX,'r','Linestyle','none');
    %     set(f,'facea',[.2]);
    plot([0 0],[-maxMeanVal maxMeanVal],'r')
    hold off
    title('\fontsize{20pt}\bf{Mean Motion During Stopping Locomotion Events}')
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