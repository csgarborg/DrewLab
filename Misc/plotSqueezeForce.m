close all
circumMeters = 0.115;
widthMeters = 0.015;
surfArea = circumMeters*widthMeters;

forceCalArb = [1800;4400;7300;10000;12000;13700];
forceCalGrams = [50;100;150;200;250;300];
forceCalNewtons = (forceCalGrams./1000)*9.81;
forceCalPascals = forceCalNewtons./surfArea;
forceCalMMHg = forceCalPascals/133.322;
m = forceCalArb\forceCalMMHg;
plot(forceCalArb,forceCalMMHg)
hold on
plot([0;forceCalArb],[0;forceCalArb]*m)

figure
trial4 = squeezeTrial;
trial4(:,1) = (trial4(:,1)-trial4(1,1))/1000;
xInd = find(trial4(:,1)>150 & trial4(:,1)<300);
xStart = xInd(1);
plotLineX = trial4(xInd,1)-trial4(xStart,1);
plotLineY = (trial4(xInd,2)-trial4(xStart,2)).*m*-1;
plot(plotLineX,plotLineY)
xlim([0 150])
ylim([-1 20])
xlabel('Time (s)')
ylabel('mmHg')