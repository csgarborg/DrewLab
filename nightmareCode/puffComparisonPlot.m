function puffComparisonPlot(tdmsDataStruct)

close all

if ~exist('tdmsDataStruct','var')
    tdmsDataStruct = openTDMS;
end

if all(tdmsDataStruct.Digital_Data.Puff == 0)
    disp('No puffs')
    return
end

tdmsDataStruct.Analog_Data.EMG = filterNeckEMGData(tdmsDataStruct.Analog_Data.EMG,str2double(tdmsDataStruct.AnalogSamplingRate_Hz_));
tdmsDataStruct.Analog_Data.Force_Sensor = filterForceSensorData(tdmsDataStruct.Analog_Data.Force_Sensor,tdmsDataStruct.AnalogSamplingRate_Hz_);

% if str2double(tdmsDataStruct.PuffDuration_s_) == 2
%     disp('skip 2s')
%     return
% end
secBeforePuff = 8;
secAfterPuff = str2double(tdmsDataStruct.PuffDuration_s_) + 12;

digitalPuffIdx = find(diff(tdmsDataStruct.Digital_Data.Puff > 0) == 1) + 1;
analogPuffIdx = round(digitalPuffIdx * (str2double(tdmsDataStruct.AnalogSamplingRate_Hz_) / str2double(tdmsDataStruct.CameraFrameratePerSecond)));

digitalDataIdxBeforePuff = secBeforePuff*str2double(tdmsDataStruct.CameraFrameratePerSecond);
digitalDataIdxAfterPuff = secAfterPuff*str2double(tdmsDataStruct.CameraFrameratePerSecond);
digitalIdxLength = digitalDataIdxBeforePuff + digitalDataIdxAfterPuff + 1;

analogDataIdxBeforePuff = secBeforePuff*str2double(tdmsDataStruct.AnalogSamplingRate_Hz_);
analoglDataIdxAfterPuff = secAfterPuff*str2double(tdmsDataStruct.AnalogSamplingRate_Hz_);
analogIdxLength = analogDataIdxBeforePuff + analoglDataIdxAfterPuff + 1;

fields = fieldnames(tdmsDataStruct.Digital_Data);
fields = setdiff(fields, {'Puff','Respiration_Sum'}, 'stable');

for i = 1:length(fields)
    figure(i)
    avgCalcMat = [];
    for n = 1:length(digitalPuffIdx)
        plot(linspace(-secBeforePuff,secAfterPuff,digitalIdxLength),tdmsDataStruct.Digital_Data.(fields{i})(digitalPuffIdx(n)-digitalDataIdxBeforePuff:digitalPuffIdx(n)+digitalDataIdxAfterPuff))
        avgCalcMat(n,:) = tdmsDataStruct.Digital_Data.(fields{i})(digitalPuffIdx(n)-digitalDataIdxBeforePuff:digitalPuffIdx(n)+digitalDataIdxAfterPuff);
        hold on
    end
    ylimAxisVals = ylim;
    rectangle('Position', [0, 0, str2double(tdmsDataStruct.PuffDuration_s_), ylimAxisVals(2)], ...
          'EdgeColor', 'r', ...
          'LineWidth', 1, ...
          'FaceColor', [1 0 0], ...
          'FaceAlpha', 0.1);
    plot(linspace(-secBeforePuff,secAfterPuff,digitalIdxLength),mean(avgCalcMat),'k')
    xlabel('Time (s)')
    ylabel(strrep(fields{i}, '_', ' '))
    puffDataStruct.digitalData.(['mean_' fields{i}]) = mean(avgCalcMat);
end

fields = fieldnames(tdmsDataStruct.Analog_Data);
fields = setdiff(fields, {'ECoG'}, 'stable');

for i = 1:length(fields)
    figure(numel(findall(0, 'Type', 'figure')) + 1)
    avgCalcMat = [];
    for n = 1:length(analogPuffIdx)
        plot(linspace(-secBeforePuff,secAfterPuff,analogIdxLength),tdmsDataStruct.Analog_Data.(fields{i})(analogPuffIdx(n)-analogDataIdxBeforePuff:analogPuffIdx(n)+analoglDataIdxAfterPuff))
        avgCalcMat(n,:) = tdmsDataStruct.Analog_Data.(fields{i})(analogPuffIdx(n)-analogDataIdxBeforePuff:analogPuffIdx(n)+analoglDataIdxAfterPuff);
        hold on
    end
    ylimAxisVals = ylim;
    rectangle('Position', [0, 0, str2double(tdmsDataStruct.PuffDuration_s_), ylimAxisVals(2)], ...
          'EdgeColor', 'r', ...
          'LineWidth', 1, ...
          'FaceColor', [1 0 0], ...
          'FaceAlpha', 0.1);
    plot(linspace(-secBeforePuff,secAfterPuff,analogIdxLength),mean(avgCalcMat),'k')
    xlabel('Time (s)')
    ylabel(strrep(fields{i}, '_', ' '))
    puffDataStruct.analogData.(['mean_' fields{i}]) = mean(avgCalcMat);
end

puffDataStruct.numPuffs = length(digitalPuffIdx);
puffDataStruct.AnimalID = tdmsDataStruct.AnimalID;
puffDataStruct.digitalTimescaleS = linspace(-secBeforePuff,secAfterPuff,digitalIdxLength);
puffDataStruct.analogTimescaleS = linspace(-secBeforePuff,secAfterPuff,analogIdxLength);
puffDataStruct.puffDurationS = tdmsDataStruct.PuffDuration_s_;
puffDataStruct.secBeforePuff = secBeforePuff;
puffDataStruct.secAfterPuff = secAfterPuff;
puffDataStruct.name = strrep(tdmsDataStruct.name, '_converted', '');

saveFileString = [tdmsDataStruct.filePath strrep(tdmsDataStruct.name, '_converted', '')];

save([saveFileString '_puffDataStruct.mat'],'puffDataStruct')

figs = findall(0, 'Type', 'figure');
savefig(figs,[saveFileString '_puffFigs.fig'])