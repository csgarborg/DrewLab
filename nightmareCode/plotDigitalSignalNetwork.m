function plotDigitalSignalNetwork(tdmsPath)

if exist('tdmsPath','var')
    tdmsDataStruct = openTDMS(tdmsPath);
else
    tdmsDataStruct = openTDMS;
end

% Sampling rate
Fs = str2double(tdmsDataStruct.CameraFrameratePerSecond);

signalsStruct = tdmsDataStruct.Digital_Data;
signalNames = fieldnames(signalsStruct);

removeFields = {'Puff','Respiration_Sum'};
signalNames = setdiff(signalNames, removeFields, 'stable');

nSignals = length(signalNames);

durationSec = 15*60;
samplesToUse = round(durationSec * Fs);

X = zeros(samplesToUse,nSignals);

for i = 1:nSignals
    
    data = signalsStruct.(signalNames{i});
    data = data(:);
    
    X(:,i) = data(1:min(samplesToUse,length(data)));
    
end

%% Correlation matrix
C = corrcoef(X);

%% Threshold
threshold = 0.6;

% Build adjacency using correlation values
A = C;
A(A < threshold) = 0;
A = A - diag(diag(A)); % remove self edges

%% Graph object
G = graph(A,signalNames,'upper');

%% Plot graph
figure

p = plot(G,'Layout','force');

% Edge weights
weights = G.Edges.Weight;

% Scale thickness
p.LineWidth = 1 + 6*weights;

% Edge labels
p.EdgeLabel = round(weights,2);

% Node formatting
p.NodeFontSize = 12;
p.MarkerSize = 8;

title(sprintf('Signal Similarity Network (r > %.2f)',threshold))

end