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

% fileNames = {'H:\21-10-21_MouseExp\211021_00','H:\21-10-22_MouseExp\211022_00','H:\21-11-01_MouseExp\211101_00','H:\21-11-09_MouseExp\211109_00','H:\21-11-12_MouseExp\211112_00','H:\21-11-16_MouseExp\211116_00','H:\21-11-17_MouseExp\211117_00','H:\21-11-19_MouseExp\211119_00','H:\21-12-03_MouseExp\211203_00','H:\21-12-16_MouseExp\211216_00'};
% vecs = {1:9,1:12,1:8,[1 2 3 5 6 7 8 9],[1:7 9:18],1:9,[1 2 3 10 11 13 14 15 16 17 18],[1:18 20:24],1:18,1:2};
% fileNames = {'H:\22-08-08_MouseExp\220808_00','H:\22-08-13_MouseExp\220813_00','H:\22-08-15_MouseExp\220815_00'};
% vecs = {1:2,3,[2 5 7]};
% for j = 1:numel(fileNames)
%     for k = vecs{j}
%         for n = 1:2
%             if k <= 9
%                 fileName = [fileNames{j} num2str(k) '_processed_Layer' num2str(n) '_'];
%             else
%                 fileName = [fileNames{j}(1:end-1) num2str(k) '_processed_Layer' num2str(n) '_'];
%             end
% 
%             versionVec = [1 2 3];
%             mode = 1;
%             if n == 1
%                 brainTF = true;
%             else
%                 brainTF = false;
%             end
%             combineMotionTracking(fileName,versionVec,mode,brainTF)
%         end
%     end
% end

% vecs = {[],[],[],[],[5 6],[],[],[],[]};
% swapTF = false;
% for j = 1:numel(fileNames)
%     for k = vecs{j}
%         if k <= 9
%             fileName = [fileNames{j} num2str(k) '_processed_Layer' num2str(n) '_'];
%             plotMotionTracking2Layer([fileNames{j} num2str(k) '_processed_Layer2_combined.mat'],[fileNames{j} num2str(k) '_processed_Layer1_combined.mat'],swapTF)
%         else
%             fileName = [fileNames{j}(1:end-1) num2str(k) '_processed_Layer' num2str(n) '_'];
%             plotMotionTracking2Layer([fileNames{j}(1:end-1) num2str(k) '_processed_Layer2_combined.mat'],[fileNames{j}(1:end-1) num2str(k) '_processed_Layer1_combined.mat'],swapTF)
%         end
%         figure(7)
%         disp('stopping to check data')
%     end
% end

% dateCell = {'211021','211022','211101','211102','211105','211109','211112','211116','211117','211119','211203','211216','220203','220209','220210','220211','220214','220221','220223','220303','220308','220309','220314','220318','220404','220406','220407','220429','220509','220511','220712','220714','220718','220719','220808','220809','220813','220815','220816','220822','220823'};
% dateCell = {'221205','221207','221208','221213'};
% dateCell = {'230404','230405','230410','230414','230421','230524','230527','230620','230622','230623','230624','230627','230630','230802'};
% dateCell = {'230828','230830'};
dateCell = {'240126'};
close all
for n = 1:size(dateCell,2)
    folderName = ['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/'];
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
        if i > 9
            runNumberStr = num2str(i);
        else
            runNumberStr = ['0' num2str(i)];
        end
        if exist(['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer1_3.mat'],'file') && exist(['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer2_3.mat'],'file')
            for k = 1:2
                fileName = ['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer' num2str(k) '_'];
                    
                versionVec = [1 2 3];
                mode = 1;
                if k == 1
                    brainTF = true;
                else
                    brainTF = false;
                end
                combineMotionTracking(fileName,versionVec,mode,brainTF)
            end
        end
    end
end




% dateCell = {'211021','211022','211101','211102','211105','211109','211112','211116','211117','211119','211203','211216','220203','220209','220210','220211','220214','220221','220223','220303','220308','220309','220314','220318','220404','220406','220407','220429','220509','220511','220712','220714','220718','220719','220808','220809','220813','220815','220816','220822','220823'};
% dateCell = {'221205','221207','221208','221213'};
% dateCell = {'230624','230627','230630','230802'};
% dateCell = {'230828','230830'};
dateCell = {'240126'};
close all
rerunCell = {};
for n = 1:size(dateCell,2)
    folderName = ['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/'];
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
        if i > 9
            runNumberStr = num2str(i);
        else
            runNumberStr = ['0' num2str(i)];
        end
        if exist(['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer1_3.mat'],'file') && exist(['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer2_3.mat'],'file')
            plotMotionTracking2Layer(['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer2_combined.mat'],['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer1_combined.mat'],false)
            figure(1)
            set(gcf,'units','normalized','outerposition',[0 0.5 .5 .5])
            figure(2)
            set(gcf,'units','normalized','outerposition',[0 0 .5 .5])
            figure(7)
            set(gcf,'units','normalized','outerposition',[0.5 0 .5 .5])
            figure(4)
            set(gcf,'units','normalized','outerposition',[0.5 0.5 .5 .5])
%             figure(1)
%             set(gcf, 'Position', get(0, 'Screensize'));
%             load(['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer1_1.mat'])
%             x1 = movementData.targetPosition(:,1);
%             y1 = movementData.targetPosition(:,2);
%             load(['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer1_2.mat'])
%             x2 = movementData.targetPosition(:,1);
%             y2 = movementData.targetPosition(:,2);
%             load(['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer1_3.mat'])
%             x3 = movementData.targetPosition(:,1);
%             y3 = movementData.targetPosition(:,2);
%             subplot(5,1,1)
%             plot(x1)
%             hold on
%             plot(x2)
%             plot(x3)
%             plot(mean([x1 x2 x3],2),'k')
%             hold off
%             subplot(5,1,2)
%             plot(y1)
%             hold on
%             plot(y2)
%             plot(y3)
%             plot(mean([y1 y2 y3],2),'k')
%             hold off
%             
%             load(['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer2_1.mat'])
%             x1 = movementData.targetPosition(:,1);
%             y1 = movementData.targetPosition(:,2);
%             load(['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer2_2.mat'])
%             x2 = movementData.targetPosition(:,1);
%             y2 = movementData.targetPosition(:,2);
%             load(['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer2_3.mat'])
%             x3 = movementData.targetPosition(:,1);
%             y3 = movementData.targetPosition(:,2);
%             subplot(5,1,3)
%             plot(x1)
%             hold on
%             plot(x2)
%             plot(x3)
%             plot(mean([x1 x2 x3],2),'k')
%             hold off
%             subplot(5,1,4)
%             plot(y1)
%             hold on
%             plot(y2)
%             plot(y3)
%             plot(mean([y1 y2 y3],2),'k')
%             hold off
%             
%             subplot(5,1,5)
%             plot(movementData.ballData(:,2))
%             xlim([0 size(movementData.ballData(:,2),1)])
            
            goodRun = 2;
            while goodRun >= 2
                goodRun = input('');
                if goodRun == 0
                    rerunCell{end+1} = [dateCell{n} '_' runNumberStr];
                elseif goodRun == 2 % orthogonal vector
                    close all
                    plotMotionTracking2Layer(['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer2_combined.mat'],['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer1_combined.mat'],true,false,false)
                    figure(7)
                elseif goodRun == 3 % switch vector direction
                    close all
                    plotMotionTracking2Layer(['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer2_combined.mat'],['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer1_combined.mat'],false,true,false)
                    figure(7)
                elseif goodRun == 4 % switch layers
                    close all
                    plotMotionTracking2Layer(['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer2_combined.mat'],['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer1_combined.mat'],false,false,true)
                    figure(7)
                elseif goodRun == 5 % orthogonal vector and switch layers
                    close all
                    plotMotionTracking2Layer(['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer2_combined.mat'],['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer1_combined.mat'],true,false,true)
                    figure(7)
                elseif goodRun == 6 % switch vector direction and switch layers
                    close all
                    plotMotionTracking2Layer(['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer2_combined.mat'],['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer1_combined.mat'],false,true,true)
                    figure(7)
                elseif goodRun == 7 % orthogonal vector and switch vector direction
                    close all
                    plotMotionTracking2Layer(['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer2_combined.mat'],['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer1_combined.mat'],true,true,false)
                    figure(7)
                elseif goodRun == 8 % switch vector direction and switch layers and switch layers
                    close all
                    plotMotionTracking2Layer(['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer2_combined.mat'],['H:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer1_combined.mat'],true,true,true)
                    figure(7)
                end
            end
            close all
        end
    end
end