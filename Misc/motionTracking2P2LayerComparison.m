function motionTracking2P2LayerComparison(dataFileStr,playMovieTF)
close all

load([dataFileStr '_processed_Layer1.mat'])
posL1 = [movementData.targetPosition(:,1) - movementData.targetPosition(1,1), movementData.targetPosition(:,2) - movementData.targetPosition(1,2)];
meanPosValL1 = [.5*(posL1(1:end-1,1) + posL1(2:end,1)),.5*(posL1(1:end-1,2) + posL1(2:end,2))];
posL1 = zipperVecs(posL1,meanPosValL1);
posL1 = posL1(2:end,:);

load([dataFileStr '_processed_Layer2.mat'])
posL2 = [movementData.targetPosition(:,1) - movementData.targetPosition(1,1), movementData.targetPosition(:,2) - movementData.targetPosition(1,2)];
meanPosValL2 = [.5*(posL2(1:end-1,1) + posL2(2:end,1)),.5*(posL2(1:end-1,2) + posL2(2:end,2))];
posL2 = zipperVecs(posL2,meanPosValL2);
if size(posL1,1) > size(posL2,1)
    posL1 = posL1(1:end-1,:);
else
    posL2 = posL2(1:end-1,:);
end

posDiff = [posL2(:,1)-posL1(:,1),posL2(:,2)-posL1(:,2)];
maxPosDiff = ceil(max(max(posDiff)));

maxMoveData = ceil(max([posL1(:,1);posL1(:,2);posL2(:,1);posL2(:,2);]));

movementData.ballData = [movementData.ballData(:,1) - movementData.ballData(1,1), movementData.ballData(:,2)];
maxBallData = ceil(max(movementData.ballData(:,2)*10))/10;
layer1Index = movementData.frames(1);
layer2Index = movementData.frames(1)+1;

if playMovieTF
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(3,4,5:8)
    x1 = [0 0];
    y1 = [-maxPosDiff maxPosDiff];
    plot((1:size(posDiff,1))*movementData.secondsPerFrame,medfilt1(posDiff(:,1),10),'r',(1:size(posDiff,1))*movementData.secondsPerFrame,medfilt1(posDiff(:,2),10),'b');
    hold on
    j = plot(x1,y1,'k--');
    j.XDataSource = 'x';
    j.YDataSource = 'y';
    hold off
    legend('X','Y','Time Marker')
    title('Position of Layer 2 vs Position of Layer 1')
    xlabel('Time (s)')
    ylabel('Distance (\mum)')
    grid on
    subplot(3,4,9:12)
    x = [0 0];
    y = [-maxBallData maxBallData];
    plot(movementData.ballData(:,1),movementData.ballData(:,2),'b-');
    hold on
    h = plot(x,y,'k--');
    h.XDataSource = 'x';
    h.YDataSource = 'y';
    hold off
    title('Ball Movement')
    xlabel('Time (s)')
    ylabel('Ball Rotation Speed')
    for n = 1:size(posL1)
%         tic
        subplot(3,4,1)
        plot(posL1(n,1),posL1(n,2),'rx',posL2(n,1),posL2(n,2),'bo','MarkerSize',20);
        axis equal square
        axis([-maxMoveData maxMoveData -maxMoveData maxMoveData])
        title('Layer Motion')
        xlabel('\mum')
        ylabel('\mum')
        legend('Layer 1','Layer 2','FontSize',14,'Location','southeast');
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        subplot(3,4,2)
        plot(posL1(n,1)-posL1(n,1),posL1(n,2)-posL1(n,2),'rx',posL2(n,1)-posL1(n,1),posL2(n,2)-posL1(n,2),'bo','MarkerSize',20);
        axis equal square
        axis([-maxMoveData maxMoveData -maxMoveData maxMoveData])
        title('Stabilized Layer Motion')
        xlabel('\mum')
        ylabel('\mum')
        legend('Layer 1','Layer 2','FontSize',14,'Location','southeast');
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        secondsFrame = n*movementData.secondsPerFrame;
        if n == 1
            subplot(3,4,3)
            imshow(im2double(imread([dataFileStr '.tif'],layer1Index)));
            title(['Layer 1 - Frame' num2str(layer1Index)])
            subplot(3,4,4)
            imshow(im2double(imread([dataFileStr '.tif'],layer2Index)));
            title(['Layer 2 - Frame' num2str(layer2Index)])
            layer2Index = layer2Index - 2;
        else
            if rem(n,2) == 1
                layer1Index = layer1Index + 2;
            else
                layer2Index = layer2Index + 2;
            end
            subplot(3,4,3)
            imshow(im2double(imread([dataFileStr '.tif'],layer1Index)));
            title(['Layer 1 - Frame ' num2str(layer1Index)])
            subplot(3,4,4)
            imshow(im2double(imread([dataFileStr '.tif'],layer2Index)));
            title(['Layer 2 - Frame ' num2str(layer2Index)])
        end
%         [~,idx] = min(abs(movementData.ballData(:,1)-secondsFrame));
%         if idx-450 < 1
%             minX = 1;
%         else
%             minX = idx-450;
%         end
%         if idx+450 > size(movementData.ballData,1)
%             maxX = size(movementData.ballData,1);
%         else
%             maxX = idx+450;
%         end
%         x = [movementData.ballData(idx,1),movementData.ballData(idx,1)];
x = [secondsFrame, secondsFrame];
        set(h,'XData',x);
        set(j,'XData',x)
%         toc
%         refreshdata
%         drawnow
% plot(movementData.ballData(:,1),movementData.ballData(:,2),'b-',[movementData.ballData(idx,1),movementData.ballData(idx,1)],[-maxBallData maxBallData],'k--');
% ylim([-maxBallData maxBallData]);
%         plot(movementData.ballData(minX:maxX,1),movementData.ballData(minX:maxX,2),'b-',[movementData.ballData(idx,1),movementData.ballData(idx,1)],[-maxBallData maxBallData],'k--');
%         axis([movementData.ballData(minX,1) movementData.ballData(maxX,1) -maxBallData maxBallData]);
        pause(.05)
        F = getframe(gcf);
        [image, ~] = frame2im(F);
        if n == 1
            imwrite(image,[dataFileStr '_comparisonOutput.tif']);
        else
            imwrite(image,[dataFileStr '_comparisonOutput.tif'],'WriteMode','append');
        end
    end
    close all
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