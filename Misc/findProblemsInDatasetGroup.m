load('movementDataLog.mat')
close all
for n = 1:size(moveDataMat,1)
    if length(moveDataMat{n,1}) == 8
        fileName = ['H:\' moveDataMat{n,1}(1:2) '-' moveDataMat{n,1}(3:4) '-' moveDataMat{n,1}(5:6) '_MouseExp\' moveDataMat{n,1}(1:end-1) '00' moveDataMat{n,1}(end) '_processe_2layerBrainInSkullDataFinal.mat'];
    else
        fileName = ['H:\' moveDataMat{n,1}(1:2) '-' moveDataMat{n,1}(3:4) '-' moveDataMat{n,1}(5:6) '_MouseExp\' moveDataMat{n,1}(1:end-2) '0' moveDataMat{n,1}(end-1:end) '_processe_2layerBrainInSkullDataFinal.mat'];
    end
    load(fileName)
    if any(movementData.ballData(1:91,2)>.01)
        disp(moveDataMat{n,1})
        figure(1)
        plot(movementData.ballData(:,1),movementData.ballData(:,2))
        close all
    end
    clear movementData
end
