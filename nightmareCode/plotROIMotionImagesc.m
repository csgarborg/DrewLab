function plotROIMotionImagesc(matFile)
S = load(matFile);
roiMotion = S.roiMotion;

fields = fieldnames(roiMotion);
fields(strcmp(fields,'time')) = [];

try
    t = roiMotion.time;
catch
    t = (1:length(roiMotion.(fields{1})))./60;
end
nROIs = length(fields);

%% ==============================
%% Build matrix
%% ==============================
dataMat = zeros(nROIs, length(t));

for i = 1:nROIs
    dataMat(i,:) = roiMotion.(fields{i});
end

%% ==============================
%% Scaling
%% ==============================
gain = 3;  % adjust as needed
offsetStep = prctile(dataMat(:), 99) * gain * 5;

offsets = (0:nROIs-1) * offsetStep;

%% ==============================
%% Plot
%% ==============================
figure; hold on

for i = 1:nROIs
    plot(t, gain * dataMat(i,:) + offsets(i), 'LineWidth', 1.2)
end

%% ==============================
%% Axis labeling (KEY PART)
%% ==============================
yticks(offsets)
yticklabels(strrep(fields,'_',' '))

xlabel('Time (s)')
ylabel('ROI')
title('ROI Motion Traces (Raw)')

box off

end