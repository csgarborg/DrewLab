function plotDigitalSignalCoherence(method,tdmsPath)

if ~exist('method','var')
    method = 'welch'; % default
end

if exist('tdmsPath','var')
    tdmsDataStruct = openTDMS(tdmsPath);
else
    tdmsDataStruct = openTDMS;
end

% Get sampling rate
Fs = str2double(tdmsDataStruct.CameraFrameratePerSecond);

% Extract signals
signalsStruct = tdmsDataStruct.Digital_Data;
signalNames = fieldnames(signalsStruct);

% Remove unwanted signals
removeFields = {'Puff','Respiration_Sum','Respiration'};
signalNames = setdiff(signalNames, removeFields, 'stable');

nSignals = length(signalNames);

% Duration to analyze (15 minutes)
durationSec = 15*60;
samplesToUse = round(durationSec * Fs);

signals = cell(nSignals,1);

for i = 1:nSignals
    data = signalsStruct.(signalNames{i});
    signals{i} = data(1:min(samplesToUse,length(data)));
end

% Welch parameters
window = hamming(1024);
noverlap = 512;
nfft = 1024;

% Chronux parameters
params.Fs = Fs;
params.tapers = [3 5]; % time-bandwidth 3, 5 tapers
params.pad = 0;
params.fpass = [0 5];

figure
tiledlayout(nSignals,nSignals,'TileSpacing','compact')

for i = 1:nSignals
    for j = 1:nSignals
        
        nexttile
        
        x = signals{i};
        y = signals{j};
        
        L = min(length(x),length(y));
        x = x(1:L);
        y = y(1:L);
        
        if strcmp(method,'welch')
            
            [Cxy,f] = mscohere(x,y,window,noverlap,nfft,Fs);
            idx = f <= 5;
            f = f(idx);
            Cxy = Cxy(idx);
            
        elseif strcmp(method,'chronux')
            
            [Cxy,~,~,~,f] = coherencyc(x,y,params);
            
        end
        
        plot(f,Cxy,'w','LineWidth',1)
        set(gca,'XScale','log')
        xlim([0.1 5])
        ylim([0 1])
        
        if i == nSignals
            xlabel('Frequency (Hz)')
        end
        
        if j == 1
            ylabel(signalNames{i},'Interpreter','none')
        end
        
        if i == 1
            title(signalNames{j},'Interpreter','none')
        end
        
    end
end

sgtitle(['Pairwise Signal Coherence (',method,')'])

end