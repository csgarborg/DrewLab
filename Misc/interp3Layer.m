close all
brainx = zeros(size(movementData.targetPosition,1),1);
brainy = zeros(size(movementData.targetPosition,1),1);
for n = 1:3
    load(['H:\22-12-05_MouseExp\221205_004_1_processed_' num2str(n) '.mat'])
    figure(1);
    set(gcf, 'Position', get(0, 'Screensize'));
    plot(movementData.targetPosition(:,1))
    offset = ginput(1);
    if offset(1) < 0
        offset(2) = 0;
    end
    brainx = brainx + movementData.targetPosition(:,1) - offset(2);
    close
    figure(1);
    set(gcf, 'Position', get(0, 'Screensize'));
    plot(movementData.targetPosition(:,2))
    offset = ginput(1);
    if offset(1) < 0
        offset(2) = 0;
    end
    brainy = brainy + movementData.targetPosition(:,2) - offset(2);
    close
end
avgBrainx = brainx ./ 3;
avgBrainy = brainy ./ 3;
timestamp1 = movementData.secondsPerFrame:movementData.secondsPerFrame:size(movementData.targetPosition,1)*movementData.secondsPerFrame;
newTimestamp = movementData.secondsPerFrame:movementData.secondsPerFrame/3:(size(movementData.targetPosition,1)*movementData.secondsPerFrame)+(movementData.secondsPerFrame*2/3);
moveDatax1 = medfilt1(interp1(timestamp1,avgBrainx,newTimestamp),8);
moveDatay1 = medfilt1(interp1(timestamp1,avgBrainy,newTimestamp),8);
durax = zeros(size(movementData.targetPosition,1),1);
duray = zeros(size(movementData.targetPosition,1),1);
for n = 1:3
    load(['H:\22-12-05_MouseExp\221205_004_2_processed_' num2str(n) '.mat'])
    figure(1);
    set(gcf, 'Position', get(0, 'Screensize'));
    plot(movementData.targetPosition(:,1))
    offset = ginput(1);
    if offset(1) < 0
        offset(2) = 0;
    end
    durax = durax + movementData.targetPosition(:,1) - offset(2);
    close
    figure(1);
    set(gcf, 'Position', get(0, 'Screensize'));
    plot(movementData.targetPosition(:,2))
    offset = ginput(1);
    if offset(1) < 0
        offset(2) = 0;
    end
    duray = duray + movementData.targetPosition(:,2) - offset(2);
    close
end
avgDurax = durax ./ 3;
avgDuray = duray ./ 3;
timestamp2 = timestamp1 + movementData.secondsPerFrame;
moveDatax2 = medfilt1(interp1(timestamp2,avgDurax,newTimestamp),8);
moveDatay2 = medfilt1(interp1(timestamp2,avgDuray,newTimestamp),8);
skullx = zeros(size(movementData.targetPosition,1),1);
skully = zeros(size(movementData.targetPosition,1),1);
for n = 1:3
    load(['H:\22-12-05_MouseExp\221205_004_3_processed_' num2str(n) '.mat'])
    figure(1);
    set(gcf, 'Position', get(0, 'Screensize'));
    plot(movementData.targetPosition(:,1))
    offset = ginput(1);
    if offset(1) < 0
        offset(2) = 0;
    end
    skullx = skullx + movementData.targetPosition(:,1) - offset(2);
    close
    figure(1);
    set(gcf, 'Position', get(0, 'Screensize'));
    plot(movementData.targetPosition(:,2))
    offset = ginput(1);
    if offset(1) < 0
        offset(2) = 0;
    end
    skully = skully + movementData.targetPosition(:,2) - offset(2);
    close
end
avgSkullx = skullx ./ 3;
avgSkully = skully ./ 3;
timestamp3 = timestamp2 + movementData.secondsPerFrame;
moveDatax3 = medfilt1(interp1(timestamp3,avgSkullx,newTimestamp),8);
moveDatay3 = medfilt1(interp1(timestamp3,avgSkully,newTimestamp),8);

% ballData = load('H:\22-12-05_MouseExp\221205_004.txt');
emgDataOnly = [ballData(:,1) ballData(:,4)];
procEMGData = filterEMGData(emgDataOnly,10000);

figure(1)
subplot(3,1,1)
plot(newTimestamp,moveDatax1*-1)
hold on
plot(newTimestamp,moveDatax2*-1)
plot(newTimestamp,moveDatax3*-1)
hold off
legend('Brain','Dural Vessel','Skull')
subplot(3,1,2)
plot(newTimestamp,moveDatay1*-1)
hold on
plot(newTimestamp,moveDatay2*-1)
plot(newTimestamp,moveDatay3*-1)
hold off
legend('Brain','Dural Vessel','Skull')
subplot(3,1,3)
plot(procEMGData(:,1),procEMGData(:,2))

for n = 1:3
    subplot(3,1,n)
    xlim([75 165])
end