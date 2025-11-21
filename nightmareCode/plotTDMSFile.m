function plotTDMSFile(tdmsPath)

close all

if ~exist('tdmsPath','var')
    [tdmsFile,tdmsFolder] = uigetfile('*.tdms');
    tdmsPath = fullfile(tdmsFolder,tdmsFile);
end
tdmsDataStruct = openTDMS(tdmsPath);

figure(1)
title('Digital Data')
fieldNames = fieldnames(tdmsDataStruct.Digital_Data);
timeVector = (1:1:length(tdmsDataStruct.Digital_Data.(fieldNames{1})))./str2double(tdmsDataStruct.CameraFrameratePerSecond);
for n = 1:length(fieldNames)
    subplot(length(fieldNames),1,n)
    plot(timeVector,tdmsDataStruct.Digital_Data.(fieldNames{n}))
    xlabel('Time (s)')
    xlim([0 3600])
    ylabel(fieldNames{n})
end

figure(2)
title('Analog and Serial Data')
fieldNames = fieldnames(tdmsDataStruct.Analog_Data);
timeVector = (1:1:length(tdmsDataStruct.Analog_Data.(fieldNames{1})))./str2double(tdmsDataStruct.AnalogSamplingRate_Hz_);
for n = 1:length(fieldNames)
    subplot(length(fieldNames)+1,1,n)
    plot(timeVector,tdmsDataStruct.Analog_Data.(fieldNames{n}))
    xlabel('Time (s)')
    xlim([0 3600])
    ylabel(fieldNames{n})
end
timeVector = (1:1:length(tdmsDataStruct.Serial_Data.Temperature))./str2double(tdmsDataStruct.TemperatureSamplingRate_Hz_);
subplot(length(fieldNames)+1,1,length(fieldNames)+1)
plot(timeVector,tdmsDataStruct.Serial_Data.Temperature)
xlabel('Time (s)')
xlim([0 3600])
ylabel('Temperature')
end