function plotTDMSFileWithVideo(tdmsPath,timeStart,timeStop)

% Wide Digital + Analog Viewer
% Large synchronized video display

%% -------------------- Input Handling --------------------

if nargin < 3
    timeStart = 0;
    timeStop  = inf;
end

if nargin < 1 || isempty(tdmsPath)
    [file,folder] = uigetfile('*.tdms','Select TDMS File');
    if isequal(file,0)
        return
    end
    tdmsPath = fullfile(folder,file);
end

%% -------------------- Locate Matching MP4 --------------------

[tdmsFolder,baseName,~] = fileparts(tdmsPath);
videoPath = fullfile(tdmsFolder, baseName + ".mp4");

if ~isfile(videoPath)
    error("Matching MP4 file not found: %s", videoPath)
end

%% -------------------- Load Data --------------------

tdmsDataStruct = openTDMS(tdmsPath);
vidObj = VideoReader(videoPath);

%% -------------------- Create UI --------------------

fig = uifigure('Name',"TDMS + Video Viewer — " + baseName,...
               'Position',[50 50 1800 1000]);

% Vertical layout: Signals (top) + Video (bottom)
mainLayout = uigridlayout(fig,[2 1]);
mainLayout.RowHeight = {'2x','3x'};   % Video larger than signals

%% ==================== TOP: SIGNALS ====================

plotPanel = uipanel(mainLayout,'Title','Digital + Analog Signals');
plotLayout = tiledlayout(plotPanel,'flow','TileSpacing','compact');

axesList = [];

%% ---------- Digital ----------
digitalFields = fieldnames(tdmsDataStruct.Digital_Data);
digitalFields = digitalFields(2:end-1);

camRate = str2double(tdmsDataStruct.CameraFrameratePerSecond);
digitalLength = length(tdmsDataStruct.Digital_Data.(digitalFields{1}));
tDigital = (0:digitalLength-1)/camRate;

for k = 1:length(digitalFields)
    ax = nexttile(plotLayout);
    plot(ax,tDigital,...
         tdmsDataStruct.Digital_Data.(digitalFields{k}),...
         'LineWidth',1)
    ylabel(ax,strrep(digitalFields{k},'_',' '))
    xlim(ax,[timeStart min(timeStop,max(tDigital))])
    axesList = [axesList ax];
end

%% ---------- Analog ----------
analogFields = fieldnames(tdmsDataStruct.Analog_Data);
analogRate = str2double(tdmsDataStruct.AnalogSamplingRate_Hz_);
analogLength = length(tdmsDataStruct.Analog_Data.(analogFields{1}));
tAnalog = (0:analogLength-1)/analogRate;

for k = 1:length(analogFields)
    ax = nexttile(plotLayout);
    plot(ax,tAnalog,...
         tdmsDataStruct.Analog_Data.(analogFields{k}),...
         'LineWidth',1)
    ylabel(ax,strrep(analogFields{k},'_',' '))
    xlim(ax,[timeStart min(timeStop,max(tAnalog))])
    axesList = [axesList ax];
end

xlabel(axesList(end),'Time (s)')
linkaxes(axesList,'x')

%% ==================== BOTTOM: VIDEO ====================

videoPanel = uipanel(mainLayout,'Title','Video Frame');

videoLayout = uigridlayout(videoPanel,[1 1]);
axVideo = uiaxes(videoLayout);
axis(axVideo,'off')

% Preserve aspect ratio properly
axVideo.DataAspectRatioMode = 'auto';
axVideo.PlotBoxAspectRatioMode = 'auto';

vidObj.CurrentTime = 0;
frame = readFrame(vidObj);
imshow(frame,'Parent',axVideo)

%% -------------------- Click-To-Show Video Frame --------------------

masterAx = axesList(1);

% Create vertical indicator line (initially invisible)
hCursor = xline(masterAx,timeStart,'r','LineWidth',2);
hCursor.Visible = 'off';

% Copy cursor to other axes
cursorCopies = gobjects(length(axesList)-1,1);
for k = 2:length(axesList)
    cursorCopies(k-1) = xline(axesList(k),timeStart,'r','LineWidth',2);
    cursorCopies(k-1).Visible = 'off';
end

% Set click callback on the figure
fig.WindowButtonDownFcn = @showFrameAtClick;

    function showFrameAtClick(~,~)

        % Get click location
        cp = masterAx.CurrentPoint;
        clickTime = cp(1,1);

        % Make sure click is inside axes
        if clickTime < masterAx.XLim(1) || clickTime > masterAx.XLim(2)
            return
        end

        % Clamp to video duration
        clickTime = max(0,min(clickTime,vidObj.Duration));

        % Show vertical indicator
        hCursor.Value = clickTime;
        hCursor.Visible = 'on';

        for i = 1:length(cursorCopies)
            cursorCopies(i).Value = clickTime;
            cursorCopies(i).Visible = 'on';
        end

        % Jump to frame (more stable than readFrame)
        frameNumber = max(1, ...
            min(floor(clickTime * vidObj.FrameRate), ...
                floor(vidObj.Duration * vidObj.FrameRate)));

        frame = read(vidObj,frameNumber);

        imshow(frame,'Parent',axVideo)
    end

%% -------------------- Callback --------------------

    function updateVideo()
        currentTime = cursor.Position(1,1);
        currentTime = max(0,min(currentTime,vidObj.Duration));

        vidObj.CurrentTime = currentTime;

        try
            frame = readFrame(vidObj);
            imshow(frame,'Parent',axVideo)
        catch
        end
    end

end
