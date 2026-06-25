function plotROIMotionAndEMG(videoFile,roiPath,procDataFile,tStart,tStop,selectedROIs)

%% ==============================
%% User settings
%% ==============================
spacing = 8;
traceScale = 1;
% artifactThresh = 5;   % SD threshold for single-frame artifacts

%% ==============================
%% Load ROI masks
%% ==============================
load(roiPath,'roiStruct')

roiMasks  = roiStruct.roiMasks;
roiLabels = roiStruct.roiLabels;

% Keep only requested ROIs
if nargin >= 6 && ~isempty(selectedROIs)

    keepIdx = find(ismember(roiLabels,selectedROIs));

    if isempty(keepIdx)
        error('None of the selected ROIs were found.')
    end

    roiMasks  = roiMasks(keepIdx);
    roiLabels = roiLabels(keepIdx);

end

nROIs = numel(roiMasks);

%% ==============================
%% Open video
%% ==============================
v = VideoReader(videoFile);
v.CurrentTime = tStart;

prevFrame = [];

motionMat = [];
tROI = [];

fprintf('Processing ROI motion...\n')

h = waitbar(0,'Computing ROI motion...');

frameCount = 0;

while hasFrame(v) && v.CurrentTime <= tStop

    frame = readFrame(v);

    gray = double(rgb2gray(frame));

    % Skip black frames
    if mean(gray(:)) < 5
        prevFrame = [];
        continue
    end

    if isempty(prevFrame)
        prevFrame = gray;
        continue
    end

    diffFrame = abs(gray - prevFrame);

    roiVals = zeros(nROIs,1,'single');

    for r = 1:nROIs
        roiVals(r) = sum(diffFrame(roiMasks{r}),'all');
    end

    motionMat(:,end+1) = roiVals;
    tROI(end+1) = v.CurrentTime;

    prevFrame = gray;

    frameCount = frameCount + 1;

    if mod(frameCount,100) == 0

        frac = (v.CurrentTime - tStart) ./ (tStop - tStart);

        waitbar(min(frac,1),h,...
            sprintf('ROI motion %.0f %%',100*frac));

    end

end

close(h)

%% ==============================
%% Remove single-frame artifacts
%% ==============================
for r = 1:nROIs

    x = double(motionMat(r,:));

    x = medfilt1(x,3);

    % xMed = medfilt1(x,3);
    % 
    % artifactIdx = abs(x - xMed) > ...
    %     artifactThresh*std(x);
    % 
    % x(artifactIdx) = xMed(artifactIdx);

    motionMat(r,:) = x;

end

%% ==============================
%% Load EMG power
%% ==============================
S = load(procDataFile,'ProcData');
ProcData = S.ProcData;

Fs = ProcData.notes.dsFs;

if isfield(ProcData,'EMG') && ...
        isfield(ProcData.EMG,'emgPower_norm')

    emgPower = ProcData.EMG.emgPower_norm(:);

elseif isfield(ProcData,'EMG') && ...
        isfield(ProcData.EMG,'emgPower')

    emgPower = ProcData.EMG.emgPower(:);

else

    error('No EMG power found.')

end

tEMG = (0:length(emgPower)-1)./Fs;

emgMask = tEMG >= tStart & tEMG <= tStop;

tEMG = tEMG(emgMask);
emgPower = emgPower(emgMask);

%% ==============================
%% Normalize ROI traces
%% ==============================
for r = 1:nROIs

    x = motionMat(r,:);

    x = log10(x + 1);
    
    x = x - median(x);

    sd = std(x);

    if sd > 0
        x = x ./ sd;
    end

    motionMat(r,:) = traceScale*x;

end

%% ==============================
%% Normalize EMG
%% ==============================
emgPower = emgPower - median(emgPower);

sd = std(emgPower);

if sd > 0
    emgPower = emgPower ./ sd;
end

emgPower = traceScale*emgPower;

%% ==============================
%% Plot
%% ==============================
emgOffset = 0;
roiOffsets = (1:nROIs)*spacing;

figure
hold on

% EMG on bottom
plot(tEMG,...
     emgPower + emgOffset,...
     'm',...
     'LineWidth',1.5)

% ROI traces
for r = 1:nROIs

    plot(tROI,...
         motionMat(r,:) + roiOffsets(r),...
         'LineWidth',1)

end

%% ==============================
%% Labels
%% ==============================
yticks([emgOffset roiOffsets])

yticklabels([{'EMG Power'}; ...
    strrep(roiLabels(:),'_',' ')])

xlabel('Time (s)')
ylabel('Signal')

title(sprintf('ROI Motion + EMG Power (%.1f - %.1f s)',...
    tStart,tStop))

xlim([tStart tStop])

grid on
box off

end