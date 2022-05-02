% load('H:\21-11-16_MouseExp\211116_008_processed_Layer2_1.mat')
% movementData.targetPosition(movementData.targetPosition(:,1)>2,1) = 0;
% movementData.targetPosition(movementData.targetPosition(:,2)>2,2) = 0;
% movementData.targetPosition(movementData.targetPosition(:,1)<-2,1) = 0;
% movementData.targetPosition(movementData.targetPosition(:,2)<-2,2) = 0;
% save('H:\21-11-16_MouseExp\211116_008_processed_Layer2_1.mat','movementData')
% save('H:\21-11-16_MouseExp\211116_008_processed_Layer2_2.mat','movementData')
% save('H:\21-11-16_MouseExp\211116_008_processed_Layer2_3.mat','movementData')
% 
% load('H:\21-11-17_MouseExp\211117_018_processed_Layer1_1.mat')
% movementData.targetPosition(movementData.targetPosition(:,1)>5,1) = 5;
% movementData.targetPosition(movementData.targetPosition(:,2)>5,2) = 5;
% movementData.targetPosition(movementData.targetPosition(:,1)<-5,1) = -5;
% movementData.targetPosition(movementData.targetPosition(:,2)<-5,2) = -5;
% save('H:\21-11-17_MouseExp\211117_018_processed_Layer1_1.mat','movementData')
% save('H:\21-11-17_MouseExp\211117_018_processed_Layer1_2.mat','movementData')
% save('H:\21-11-17_MouseExp\211117_018_processed_Layer1_3.mat','movementData')

fileNames = {'H:\21-10-21_MouseExp\211021_00','H:\21-10-22_MouseExp\211022_00','H:\21-11-01_MouseExp\211101_00','H:\21-11-09_MouseExp\211109_00','H:\21-11-12_MouseExp\211112_00','H:\21-11-16_MouseExp\211116_00','H:\21-11-17_MouseExp\211117_00','H:\21-11-19_MouseExp\211119_00','H:\21-12-03_MouseExp\211203_00','H:\21-12-16_MouseExp\211216_00'};
vecs = {1:9,1:12,1:8,[1 2 3 5 6 7 8 9],[1:7 9:18],1:9,[1 2 3 10 11 13 14 15 16 17 18],[1:18 20:24],1:18,1:2};
for j = 1:numel(fileNames)
    for k = vecs{j}
        for n = 1:2
            if k <= 9
                fileName = [fileNames{j} num2str(k) '_processed_Layer' num2str(n) '_'];
            else
                fileName = [fileNames{j}(1:end-1) num2str(k) '_processed_Layer' num2str(n) '_'];
            end

            versionVec = [1 2 3];
            mode = 1;
            if n == 1
                brainTF = true;
            else
                brainTF = false;
            end
            combineMotionTracking(fileName,versionVec,mode,brainTF)
        end
    end
end

% vecs = {[],[],[],[],[5 6],[],[],[],[]};
swapTF = false;
for j = 1:numel(fileNames)
    for k = vecs{j}
        if k <= 9
            fileName = [fileNames{j} num2str(k) '_processed_Layer' num2str(n) '_'];
            plotMotionTracking2Layer([fileNames{j} num2str(k) '_processed_Layer2_combined.mat'],[fileNames{j} num2str(k) '_processed_Layer1_combined.mat'],swapTF)
        else
            fileName = [fileNames{j}(1:end-1) num2str(k) '_processed_Layer' num2str(n) '_'];
            plotMotionTracking2Layer([fileNames{j}(1:end-1) num2str(k) '_processed_Layer2_combined.mat'],[fileNames{j}(1:end-1) num2str(k) '_processed_Layer1_combined.mat'],swapTF)
        end
        figure(7)
        disp('stopping to check data')
    end
end