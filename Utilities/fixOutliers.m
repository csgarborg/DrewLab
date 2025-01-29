% filename = 'H:\21-12-16_MouseExp\211216_002_processed_Layer2_3.mat';
% close all
% load(filename)
% for m = 1:2
%     for n = 1:size(movementData.targetPosition,1)
%         if abs(movementData.targetPosition(n,m)) > 2
%             movementData.targetPosition(n,m) = movementData.targetPosition(n-1,m);
%         end
%     end
% end
% save(filename,'movementData')
% plotMotionTracking2Layer(filename);figure(3)
x = load('H:\21-12-16_MouseExp\211216_001.txt');
procData = filterEMGData(x(:,[1 4]),10000);
for n = 1:4
    load(['H:\21-12-16_MouseExp\211216_001_processed_Layer2_' num2str(n) '.mat'])
    movementData.emgData = procData;
    save(['H:\21-12-16_MouseExp\211216_001_processed_Layer2_' num2str(n) '.mat'],'movementData')
    load(['H:\21-12-16_MouseExp\211216_001_processed_Layer1_' num2str(n) '.mat'])
    movementData.emgData = procData;
    save(['H:\21-12-16_MouseExp\211216_001_processed_Layer1_' num2str(n) '.mat'],'movementData')
end