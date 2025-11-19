close all
% dateCell = {'211021','211022','211101','211102','211105','211109','211112','211116','211117','211119','211203','211216','220203','220209','220210','220211','220214','220221','220223','220303','220308','220309','220314','220318','220404','220406','220407','220429','220509','220511','220712','220714','220718','220719','220808','220809','220813','220815','220816','220822','220823','221205','221207','221208','221213'};
figure(1)
load('movementDataLog.mat')
moveDataMatSkull = moveDataMat;
for n = 1:size(moveDataMat,1)
    %     folderName = ['I:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/'];
    %     disp(folderName)
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
    %         fileName = ['I:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processed_Layer2_combined.mat'];
    dateCell = moveDataMat{n,1};
    if length(dateCell) > 8
        runNumberStr = num2str(dateCell(8:9));
    else
        runNumberStr = ['0' num2str(dateCell(end))];
    end
    fileName = ['D:/' dateCell(1:2) '-' dateCell(3:4) '-' dateCell(5:6) '_MouseExp/' dateCell(1:6) '_0' runNumberStr '_processed_Layer2_combined.mat'];
    if exist(fileName,'file')
        load(fileName)
        motionVec = pcaMotionAnalysis(movementData.targetPosition);
        plot(movementData.xLocMicrons,-movementData.yLocMicrons,'ko','MarkerSize',5)
        hold on
        quiver(movementData.xLocMicrons,-movementData.yLocMicrons,motionVec(1),motionVec(2),150,'Color','k','LineWidth',2.5,'MaxHeadSize',40)
        moveDataMatSkull{n,7} = motionVec;
    else
        disp('not found')
    end
    %     end
end
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
xlim([-8000 8000])
ylim([-8000 8000])
hold off