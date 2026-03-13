function plotDigitalEventCoincidence(tdmsPath)

if exist('tdmsPath','var')
    tdmsDataStruct = openTDMS(tdmsPath);
else
    tdmsDataStruct = openTDMS;
end

Fs = str2double(tdmsDataStruct.CameraFrameratePerSecond);

signalsStruct = tdmsDataStruct.Digital_Data;
signalNames = fieldnames(signalsStruct);

removeFields = {'Puff','Respiration_Sum','Respiration'};
signalNames = setdiff(signalNames, removeFields, 'stable');

nSignals = length(signalNames);

durationSec = 15*60;
samplesToUse = round(durationSec * Fs);

signals = cell(nSignals,1);

for i = 1:nSignals
    data = signalsStruct.(signalNames{i});
    data = data(:);
    signals{i} = data(1:min(samplesToUse,length(data)));
end

%% detect events

events = cell(nSignals,1);

for i = 1:nSignals
    
    s = signals{i};
    
    thresh = mean(s) + 2*std(s);
    
    events{i} = find(s > thresh);
    
end

%% coincidence matrix

coincMatrix = zeros(nSignals);

window = round(Fs*0.1); % 100 ms window

for i = 1:nSignals
    for j = 1:nSignals
        
        e1 = events{i};
        e2 = events{j};
        
        matches = 0;
        
        for k = 1:length(e1)
            
            if any(abs(e2 - e1(k)) < window)
                matches = matches + 1;
            end
            
        end
        
        coincMatrix(i,j) = matches/length(e1);
        
    end
end

%% plot

figure
imagesc(coincMatrix)
for i = 1:nSignals
    for j = 1:nSignals
        text(j,i,sprintf('%.2f',coincMatrix(i,j)), ...
            'HorizontalAlignment','center', ...
            'Color','k');
    end
end
colorbar
axis square
clim([0 1])

ax = gca;
ax.XTick = 1:nSignals;
ax.YTick = 1:nSignals;
ax.XTickLabel = signalNames;
ax.YTickLabel = signalNames;
ax.TickLabelInterpreter = 'none';

title('Event Coincidence Matrix')

end