function plotDigitalSignalPowerSpectra(tdmsPath)

if exist('tdmsPath','var')
    tdmsDataStruct = openTDMS(tdmsPath);
else
    tdmsDataStruct = openTDMS;
end

Fs = str2double(tdmsDataStruct.CameraFrameratePerSecond);

signalsStruct = tdmsDataStruct.Digital_Data;
signalNames = fieldnames(signalsStruct);

removeFields = {'Puff','Respiration_Sum'};
signalNames = setdiff(signalNames, removeFields, 'stable');

nSignals = length(signalNames);

durationSec = 15*60;
samplesToUse = round(durationSec*Fs);

signals = cell(nSignals,1);

for i = 1:nSignals
    
    data = signalsStruct.(signalNames{i});
    data = data(:);
    
    signals{i} = data(1:min(samplesToUse,length(data)));
    
end

%% Welch parameters

window = hamming(4096);
noverlap = 2048;
nfft = 4096;

%% Compute spectra

figure
hold on

spectra = cell(nSignals,1);

for i = 1:nSignals
    
    x = signals{i};
    
    [Pxx,f] = pwelch(x,window,noverlap,nfft,Fs);
    
    spectra{i} = Pxx;
    
    plot(f,Pxx,'LineWidth',1.5)
    
end

set(gca,'XScale','log')
set(gca,'YScale','log')

xlim([0.1 10])

xlabel('Frequency (Hz)')
ylabel('Power')

plotNames = strrep(signalNames,'_',' ');
legend(plotNames,'Location','best')

title('Signal Power Spectrum Comparison')

%% Stacked spectra

figure
tiledlayout(nSignals,1)

for i=1:nSignals
    
    nexttile
    
    plot(f,spectra{i},'LineWidth',2)
    
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    
    xlim([0.1 10])
    
    ylabel(plotNames{i})
    
    if i==nSignals
        xlabel('Frequency (Hz)')
    end
    
end

sgtitle('Signal Power Spectra')

end