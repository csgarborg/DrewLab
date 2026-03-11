function analyzeDigitalSignals(tdmsPath)

close all

fprintf('\nLoading TDMS file...\n')

if ~exist('tdmsPath','var')
    [file, loc] = uigetfile('H:\SGNMData\*.tdms', 'Select a File');
    tdmsPath = fullfile(loc,file);
end

%% Coherence plots
fprintf('Generating coherence spectra...\n')
plotDigitalSignalCoherence('welch',tdmsPath)
% plotDigitalSignalCoherence('chronux',tdmsPath)

fprintf('Generating coherence heatmap...\n')
plotDigitalSignalCoherenceHeatmap('welch',tdmsPath)
% plotDigitalSignalCoherenceHeatmap('chronux',tdmsPath)

%% Correlation + cross correlation
fprintf('Generating correlation and cross-correlation plots...\n')
plotDigitalSignalCorrelation(tdmsPath)

%% Event coincidence
fprintf('Generating event coincidence matrix...\n')
plotDigitalEventCoincidence(tdmsPath)

%% PCA analysis
fprintf('Generating PCA plots...\n')
plotDigitalSignalPCA(tdmsPath)

%% Clustering dendrogram
fprintf('Generating signal similarity dendrogram...\n')
plotDigitalSignalDendrogram(tdmsPath)

%% Clustered correlation heatmap
fprintf('Generating clustered correlation heatmap...\n')
plotDigitalSignalClusteredCorrelation(tdmsPath)

%% Signal similarity network
fprintf('Generating signal similarity network...\n')
plotDigitalSignalNetwork(tdmsPath)

%% Signal overlay
fprintf('Generating z-scored signal overlay...\n')
plotDigitalSignalOverlay(tdmsPath)

%% Lag matrix
fprintf('Generating lag matrix...\n')
plotDigitalSignalLagMatrix(tdmsPath)

%% Motion triggered average
fprintf('Generating motion triggered averages...\n')
plotDigitalMotionTriggeredAverage(tdmsPath)

%% Power spectrum comparison
fprintf('Generating power spectrum comparison...\n')
plotDigitalSignalPowerSpectra(tdmsPath)

fprintf('\nAll analyses complete.\n')
end