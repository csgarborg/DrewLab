close all
fileName = 'C:\Workspace\Code\DrewLab\movementDataLog.mat';
load(fileName)
if contains(fileName,'Skull') || contains(fileName,'Brain')
    xLoc = [];
    yLoc = [];
    UMic = [];
    VMic = [];
    a = {[] [] []};
    b = {[] [] []};
    c = {[] [] []};
    d = {[] [] []};
    e = {[] [] []};
    for i = 1:2:5
        skullMeanX = {[],[],[],[],[]};
        skullMeanY = {[],[],[],[],[]};
        skullStdX = {[],[],[],[],[]};
        skullStdY = {[],[],[],[],[]};
        brainMeanX = {[],[],[],[],[]};
        brainMeanY = {[],[],[],[],[]};
        brainStdX = {[],[],[],[],[]};
        brainStdY = {[],[],[],[],[]};
        figure(ceil(i/2))
        hold on
        for n = 1:size(moveDataMat,1)
            if strcmp(moveDataMat{n,3},'21')
                arrowColor = 'r';
                a{ceil(i/2)}(end+1) = sqrt(moveDataMat{n,7}(i)^2 + moveDataMat{n,7}(i+1)^2);
                catIdx = 0;
            elseif strcmp(moveDataMat{n,3},'27')
                b{ceil(i/2)}(end+1) = sqrt(moveDataMat{n,7}(i)^2 + moveDataMat{n,7}(i+1)^2);
                arrowColor = 'k';
                catIdx = 1;
            elseif strcmp(moveDataMat{n,3},'28')
                arrowColor = 'g';
                c{ceil(i/2)}(end+1) = sqrt(moveDataMat{n,7}(i)^2 + moveDataMat{n,7}(i+1)^2);
                catIdx = 2;
            elseif strcmp(moveDataMat{n,3},'39')
                arrowColor = 'b';
                d{ceil(i/2)}(end+1) = sqrt(moveDataMat{n,7}(i)^2 + moveDataMat{n,7}(i+1)^2);
                catIdx = 3;
            else
                arrowColor = 'm';
                e{ceil(i/2)}(end+1) = sqrt(moveDataMat{n,7}(i)^2 + moveDataMat{n,7}(i+1)^2);
                catIdx = 4;
            end
            plot(moveDataMat{n,4},-moveDataMat{n,5},[arrowColor 'o'],'MarkerSize',5)
            quiver(moveDataMat{n,4},-moveDataMat{n,5},moveDataMat{n,7}(i),moveDataMat{n,7}(i+1),150,'Color',arrowColor,'LineWidth',2.5,'MaxHeadSize',40)
        end
        xlim([-3000 3000])
        ylim([-3000 3000])
        plot(0,0,'kx','MarkerSize',12)
        plot(0,-2600,'kx','MarkerSize',12)
        quiver(-2500,2500,1,0,150,'LineWidth',2.5,'MaxHeadSize',40)
        title('Rostral, \mum')
        ylabel('Left, \mum')
        xlabel('Caudal, \mum')
        text(55,55,'Bregma')
        text(55,-2500,'Lambda')
        text(-2500,2400,'1 \mum')
        rectangle('Position',[-2650 2200 500 600])
    end
    
    catIdx = [zeros(length(skullMeanX{1}),1);ones(length(skullMeanX{2}),1);ones(length(skullMeanX{3}),1)+1;ones(length(skullMeanX{4}),1)+2;ones(length(skullMeanX{5}),1)+3];
    
    figure(4)
    plotSpread({a{1},b{1},c{1},d{1},e{1}},'categoryIdx',catIdx,'categoryColors',{'r','k','g','b','m'},'categoryLabels',{'SW21','SW27','SW28','SW39','SW25'});
    ylabel('Movement Magnitude \mum')
    xlabel('Mouse')
    title('Mean of top 20%')
    ylim([0 15])
    
    figure(5)
    plotSpread({a{2},b{2},c{2},d{2},e{2}},'categoryIdx',catIdx,'categoryColors',{'r','k','g','b','m'},'categoryLabels',{'SW21','SW27','SW28','SW39','SW25'});
    ylabel('Movement Magnitude \mum')
    xlabel('Mouse')
    title('75 Percentile')
    ylim([0 15])
    
    figure(6)
    plotSpread({a{3},b{3},c{3},d{3},e{3}},'categoryIdx',catIdx,'categoryColors',{'r','k','g','b','m'},'categoryLabels',{'SW21','SW27','SW28','SW39','SW25'});
    ylabel('Movement Magnitude \mum')
    xlabel('Mouse')
    title('95 Percentile')
    ylim([0 15])
else
    xLoc = [];
    yLoc = [];
    UMic = [];
    VMic = [];
    a = {[] [] []};
    b = {[] [] []};
    c = {[] [] []};
    d = {[] [] []};
    e = {[] [] []};
    for i = 1:2:5
        skullMeanX = {[],[],[],[],[]};
        skullMeanY = {[],[],[],[],[]};
        skullStdX = {[],[],[],[],[]};
        skullStdY = {[],[],[],[],[]};
        brainMeanX = {[],[],[],[],[]};
        brainMeanY = {[],[],[],[],[]};
        brainStdX = {[],[],[],[],[]};
        brainStdY = {[],[],[],[],[]};
        figure(ceil(i/2))
        hold on
        for n = 1:size(moveDataMat,1)
            if strcmp(moveDataMat{n,3},'21a')
                arrowColor = 'r';
                a{ceil(i/2)}(end+1) = sqrt(moveDataMat{n,7}(i)^2 + moveDataMat{n,7}(i+1)^2);
                catIdx = 0;
            elseif strcmp(moveDataMat{n,3},'27a')
                b{ceil(i/2)}(end+1) = sqrt(moveDataMat{n,7}(i)^2 + moveDataMat{n,7}(i+1)^2);
                arrowColor = 'm';
                catIdx = 1;
            elseif strcmp(moveDataMat{n,3},'28a')
                arrowColor = 'g';
                c{ceil(i/2)}(end+1) = sqrt(moveDataMat{n,7}(i)^2 + moveDataMat{n,7}(i+1)^2);
                catIdx = 2;
            elseif strcmp(moveDataMat{n,3},'39a')
                arrowColor = 'b';
                d{ceil(i/2)}(end+1) = sqrt(moveDataMat{n,7}(i)^2 + moveDataMat{n,7}(i+1)^2);
                catIdx = 3;
            else
                arrowColor = 'k';
                e{ceil(i/2)}(end+1) = sqrt(moveDataMat{n,7}(i)^2 + moveDataMat{n,7}(i+1)^2);
                catIdx = 4;
            end
            plot(moveDataMat{n,4},-moveDataMat{n,5},[arrowColor 'o'],'MarkerSize',5)
            quiver(moveDataMat{n,4},-moveDataMat{n,5},moveDataMat{n,7}(i),moveDataMat{n,7}(i+1),150,'Color',arrowColor,'LineWidth',2.5,'MaxHeadSize',40)
            skullMeanX{catIdx+1}(end+1) = moveDataMat{n,8};
            skullMeanY{catIdx+1}(end+1) = moveDataMat{n,9};
            skullStdX{catIdx+1}(end+1) = moveDataMat{n,10};
            skullStdY{catIdx+1}(end+1) = moveDataMat{n,11};
            brainMeanX{catIdx+1}(end+1) = moveDataMat{n,12};
            brainMeanY{catIdx+1}(end+1) = moveDataMat{n,13};
            brainStdX{catIdx+1}(end+1) = moveDataMat{n,14};
            brainStdY{catIdx+1}(end+1) = moveDataMat{n,15};
        end
        xlim([-8000 8000])
        ylim([-8000 8000])
        plot(0,0,'kx','MarkerSize',12)
        plot(0,-2600,'kx','MarkerSize',12)
        quiver(-2500,2500,1,0,150,'LineWidth',2.5,'MaxHeadSize',40)
        title('Rostral, \mum')
        ylabel('Left, \mum')
        xlabel('Caudal, \mum')
        text(55,55,'Bregma')
        text(55,-2500,'Lambda')
        text(-2500,2400,'1 \mum')
        rectangle('Position',[-2650 2200 500 600])
    end
    
    catIdx = [zeros(length(skullMeanX{1}),1);ones(length(skullMeanX{2}),1);ones(length(skullMeanX{3}),1)+1;ones(length(skullMeanX{4}),1)+2;ones(length(skullMeanX{5}),1)+3];
    
    figure(4)
    plotSpread({a{1},b{1},c{1},d{1},e{1}},'categoryIdx',catIdx,'categoryColors',{'r','k','g','b','m'},'categoryLabels',{'SW21','SW27','SW28','SW39','SW25'});
    ylabel('Movement Magnitude \mum')
    xlabel('Mouse')
    title('Mean of top 20%')
    ylim([0 15])
    
    figure(5)
    plotSpread({a{2},b{2},c{2},d{2},e{2}},'categoryIdx',catIdx,'categoryColors',{'r','k','g','b','m'},'categoryLabels',{'SW21','SW27','SW28','SW39','SW25'});
    ylabel('Movement Magnitude \mum')
    xlabel('Mouse')
    title('75 Percentile')
    ylim([0 15])
    
    figure(6)
    plotSpread({a{3},b{3},c{3},d{3},e{3}},'categoryIdx',catIdx,'categoryColors',{'r','k','g','b','m'},'categoryLabels',{'SW21','SW27','SW28','SW39','SW25'});
    ylabel('Movement Magnitude \mum')
    xlabel('Mouse')
    title('95 Percentile')
    ylim([0 15])
    
    
    % plot(ones(length(a)),a,'ro')
    % hold on
    % plot(ones(length(b))*2,b,'ko')
    % plot(ones(length(c))*3,c,'go')
    % plot(ones(length(d))*4,d,'bo')
    % plot(ones(length(e))*5,e,'mo')
    % hold off
    % xlim([0 6])
    % ylabel('Movement Magnitude \mum')
    % xlabel('Mouse')
    
    figure(7)
    plotSpread(skullMeanX,'categoryIdx',catIdx,'categoryColors',{'r','k','g','b','m'},'categoryLabels',{'SW21','SW27','SW28','SW39','SW25'});
    title('Skull Mean CI95 Lateral/Medial')
    ylabel('\mum')
    xlabel('Mouse')
    ylim([0 5])
    
    figure(8)
    plotSpread(skullMeanY,'categoryIdx',catIdx,'categoryColors',{'r','k','g','b','m'},'categoryLabels',{'SW21','SW27','SW28','SW39','SW25'});
    title('Skull Mean CI95 Rostral/Caudal')
    ylabel('\mum')
    xlabel('Mouse')
    ylim([0 5])
    
    figure(9)
    plotSpread(skullStdX,'categoryIdx',catIdx,'categoryColors',{'r','k','g','b','m'},'categoryLabels',{'SW21','SW27','SW28','SW39','SW25'});
    title('Skull Std CI95 Lateral/Medial')
    ylabel('\mum')
    xlabel('Mouse')
    ylim([0 5])
    
    figure(10)
    plotSpread(skullStdY,'categoryIdx',catIdx,'categoryColors',{'r','k','g','b','m'},'categoryLabels',{'SW21','SW27','SW28','SW39','SW25'});
    title('Skull Std CI95 Rostral/Caudal')
    ylabel('\mum')
    xlabel('Mouse')
    ylim([0 5])
    
    figure(11)
    plotSpread(brainMeanX,'categoryIdx',catIdx,'categoryColors',{'r','k','g','b','m'},'categoryLabels',{'SW21','SW27','SW28','SW39','SW25'});
    title('Brain Mean CI95 Lateral/Medial')
    ylabel('\mum')
    xlabel('Mouse')
    ylim([0 5])
    
    figure(12)
    plotSpread(brainMeanY,'categoryIdx',catIdx,'categoryColors',{'r','k','g','b','m'},'categoryLabels',{'SW21','SW27','SW28','SW39','SW25'});
    title('Brain Mean CI95 Rostral/Caudal')
    ylabel('\mum')
    xlabel('Mouse')
    ylim([0 5])
    
    figure(13)
    plotSpread(brainStdX,'categoryIdx',catIdx,'categoryColors',{'r','k','g','b','m'},'categoryLabels',{'SW21','SW27','SW28','SW39','SW25'});
    title('Brain Std CI95 Lateral/Medial')
    ylabel('\mum')
    xlabel('Mouse')
    ylim([0 5])
    
    figure(14)
    plotSpread(brainStdY,'categoryIdx',catIdx,'categoryColors',{'r','k','g','b','m'},'categoryLabels',{'SW21','SW27','SW28','SW39','SW25'});
    title('Brain Std CI95 Rostral/Caudal')
    ylabel('\mum')
    xlabel('Mouse')
    ylim([0 5])
end