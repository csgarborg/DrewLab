% dateCell = {'211022','211101','211102','211105','211109','211112','211116','211117','211119','211203','211216','220203','220209','220210','220211','220214','220221','220223','220303','220308','220309','220314','220318','220404','220406','220407','220429','220509','220511','220712','220714','220718','220719','220808','220809','220813','220815','220816','220822','220823'};
% dateCell = {'230404','230405','230410','230414','230421','230524','230527','230620','230622','230623','230624','230627','230630','230802'};
dateCell = {'211216'};
% '211021',
close all
for n = 1:size(dateCell,2)
    folderName = ['D:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/'];
    fileList = dir(folderName);
    fileNamesCell = struct2cell(fileList);
    fileNames = fileNamesCell(1,:);
    maxRun = 0;
    for i = 1:size(fileNames,2)
        if contains(fileNames{i},dateCell{n}) && str2double(fileNames{i}(8:10)) > maxRun
             maxRun = str2double(fileNames{i}(8:10));
        end
    end
    for i = 2%i = 1:maxRun
        if i > 9
            runNumberStr = num2str(i);
        else
            runNumberStr = ['0' num2str(i)];
        end
        if exist(['D:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer1_3.mat'],'file') && exist(['D:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer2_3.mat'],'file')
            for j = 1:5
                load(['D:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer1_' num2str(j) '.mat']);
                dataFieldNames = fieldnames(movementData);
                if any(strcmp('targetPositionNoBaseline',dataFieldNames))
                    movementData.targetPosition = movementData.targetPositionNoBaseline;
                end
                figure(1);
                set(gcf, 'Position', get(0, 'Screensize'));
                subplot(2,1,1)
                plot(movementData.targetPosition(:,1))
                title(['x baseline shift Layer 1 D:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer1_' num2str(j) '.mat'])
                subplot(2,1,2)
                plot(movementData.ballData(:,2))
                offsetx = ginput(1);
                if offsetx(1) < 0
                    offsetx(2) = 0;
%                     disp('No change')
                end
                close
                figure(1)
                set(gcf, 'Position', get(0, 'Screensize'));
                subplot(2,1,1)
                plot(movementData.targetPosition(:,2))
                title(['y baseline shift Layer 1 D:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer1_' num2str(j) '.mat'])
                subplot(2,1,2)
                plot(movementData.ballData(:,2))
                offsety = ginput(1);
                if offsety(1) < 0
                    offsety(2) = 0;
%                     disp('No change')
                end
                close
                movementData.targetPositionNoBaseline = movementData.targetPosition;
                movementData.targetPosition = [movementData.targetPosition(:,1) - offsetx(2), movementData.targetPosition(:,2) - offsety(2)];
                save(['D:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer1_' num2str(j) '.mat'],'movementData');
            end
            for j = 1:5
                load(['D:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer2_' num2str(j) '.mat']);
                figure(1);
                set(gcf, 'Position', get(0, 'Screensize'));
                subplot(2,1,1)
                plot(movementData.targetPosition(:,1))
                title(['x baseline shift Layer 2 D:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer2_' num2str(j) '.mat'])
                subplot(2,1,2)
                plot(movementData.ballData(:,2))
                offsetx = ginput(1);
                if offsetx(1) < 0
                    offsetx(2) = 0;
%                     disp('No change')
                end
                close
                figure(1);
                set(gcf, 'Position', get(0, 'Screensize'));
                subplot(2,1,1)
                plot(movementData.targetPosition(:,2))
                title(['y baseline shift Layer 2 D:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer2_' num2str(j) '.mat'])
                subplot(2,1,2)
                plot(movementData.ballData(:,2))
                offsety = ginput(1);
                if offsety(1) < 0
                    offsety(2) = 0;
%                     disp('No change')
                end
                close
                movementData.targetPositionNoBaseline = movementData.targetPosition;
                movementData.targetPosition = [movementData.targetPosition(:,1) - offsetx(2), movementData.targetPosition(:,2) - offsety(2)];
                save(['D:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer2_' num2str(j) '.mat'],'movementData');
            end
        else
            disp([dateCell{n} '_' num2str(i)])
        end
    end
end
