dateCell = {'211021','211022','211101','211102','211105','211109','211112','211116','211117','211119','211203','211216','220203','220209','220210','220211','220214','220221','220223','220303','220308','220309','220314','220318','220404','220406','220407','220429','220509','220511','220712','220714','220718','220719','220808','220809','220813','220815','220816','220822','220823','221205','221207','221208','221213'};
close all
dataCell = {};
for n = 1:size(dateCell,2)
    folderName = ['I:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/'];
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
        if exist(['I:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processe_2layerBrainInSkullDataFinal.mat'],'file')
            load(['I:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processe_2layerBrainInSkullDataFinal.mat'])
            figure('units','normalized','outerposition',[0 0 1 1])
            subplot(3,1,1)
            plot(movementData.emgData(:,1),movementData.emgData(:,2))
            subplot(3,1,2)
            plot((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame,movementData.targetPosition(:,1))
            subplot(3,1,3)
            plot((1:size(movementData.targetPosition,1))*movementData.secondsPerFrame,movementData.targetPosition(:,2))
            title(['I:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processe_2layerBrainInSkullDataFinal.mat'])
            goodRun = input('');
            if goodRun > 0
                dataCell(end+1,1:5) = {['I:/' dateCell{n}(1:2) '-' dateCell{n}(3:4) '-' dateCell{n}(5:6) '_MouseExp/' dateCell{n} '_0' runNumberStr '_processe_2layerBrainInSkullDataFinal.mat'],movementData.emgData,movementData.targetPosition,movementData.secondsPerFrame,goodRun};
            end
            close all
        end
    end
end




% mouseInfo = {'18','M','12/1/20','60.5';'21','F','12/1/20','36';'24','F','6/14/21','31.2';'25','F','6/14/21','42.2';'27','M','6/14/21','33.5';'28','M','6/14/21','40.8';'34','M','9/20/21','53.7';'35','M','9/20/21','51.1';'39','F','9/20/21','31';'40','F','9/20/21','40.1';'44','M','12/10/21','43.1';'45','M','12/10/21','40.2';'46','M','12/10/21','40.5';'48','F','12/10/21','42.2';'50','F','12/10/21','28.5';'53','M','3/4/22','53.1';'54','M','3/4/22','48';'59','F','3/4/22','34.1';'60','F','3/4/22','40.2';'61','F','3/4/22','47.7';'62','F','3/4/22','33.7';'63','M','6/3/22','49';'65','M','6/3/22','39.2';'68','F','6/3/22','47.1'};
% 
% moveDataMatNew = cell(316,18);
% moveDataMatNew(1:316,1:15) = moveDataMat;
% for n = 1:316
%     cmpre = strcmp(moveDataMat{n,3},mouseInfo(:,1));
%     moveDataMatNew(n,16:18) = mouseInfo(cmpre,2:4);
% end