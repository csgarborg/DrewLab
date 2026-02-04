function procData = filterForceSensorData(forceSensorData,samplingRate)

cutoffFrequency = .5;
order = 4;

[b, a] = butter(order, cutoffFrequency/(str2double(samplingRate)/2), 'high');
forceSensorData = forceSensorData - mean(forceSensorData);
procData = abs(filtfilt(b, a, forceSensorData));

% plot(procData-2)
% hold on
% plot(forceSensorData)