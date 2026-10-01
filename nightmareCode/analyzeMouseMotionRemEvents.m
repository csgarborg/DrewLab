function analyzeMouseMotionRemEvents(excelPath, varargin)
%ANALYZEREMEVENTS  Run motion analysis (std-dev / motion-energy, overall
%   and windowed) on only the REM-scored segments of each recording listed
%   in an Excel scoring file.
%
%   analyzeRemEvents(excelPath) reads Sheet 2 of excelPath (same format
%   used elsewhere in your pipeline: col 1 = TDMS path, remaining columns
%   = flattened start/stop second pairs for each REM bout), converts each
%   TDMS path to the matching .mp4 path, and for every REM segment in
%   every recording computes:
%       - an overall std-dev projection across the whole segment
%       - an overall motion-energy map (std of frame differences) across
%         the whole segment
%       - (optionally) windowed versions of both, for sub-segment
%         time-resolved detail
%   All outputs are written as 32-bit float TIFFs, directly openable in
%   ImageJ/Fiji, using a memory-safe single-pass streaming approach (no
%   full frame stack is ever held in memory).
%
%   FOLDER LAYOUT
%   All output goes in a folder next to the Excel file, named after the
%   Excel file itself:
%       <excelFolder>/<excelBaseName>/
%           stdDev_<mp4name>_<startSec>_<stopSec>.tif       (overall, per event)
%           diffStdDev_<mp4name>_<startSec>_<stopSec>.tif   (overall, per event)
%           <mp4name>_<startSec>_<stopSec>/                 (only if WindowSize set)
%               stdDev_####_frames_X-Y.tif
%               diffStdDev_####_diffs_X-Y.tif
%
%   Name-Value options:
%       'WindowSize'    - if set (e.g. 300), also computes windowed
%                         std-dev / motion-energy maps within each REM
%                         segment, saved into a per-event subfolder as
%                         described above. Default: [] (off)
%       'WakeTF'        - mirrors the wakeTF flag from your existing
%                         parsing code (adds a stop+30 column). NOTE: only
%                         columns 1-2 of each segment row (start/stop) are
%                         used to define the REM window actually analyzed
%                         here; the wakeTF column is parsed for
%                         compatibility with your table format but is not
%                         separately processed as its own segment. Flag
%                         if you need that handled differently.
%                         Default: false
%       'ProgressEvery' - print progress every N frames (default: 2000)
%       'Precision'     - 'single' or 'double' for accumulators
%                         (default: 'single')
%
%   Example:
%       analyzeRemEvents('C:\data\REM_scoring.xlsx', 'WindowSize', 300);

    p = inputParser;
    addParameter(p, 'WindowSize', []);
    addParameter(p, 'WakeTF', false);
    addParameter(p, 'ProgressEvery', 2000);
    addParameter(p, 'Precision', 'single');
    parse(p, varargin{:});
    opts = p.Results;
    wakeTF = opts.WakeTF;

    if ~isfile(excelPath)
        error('analyzeRemEvents:fileNotFound', 'Cannot find Excel file: %s', excelPath);
    end

    % ---- Set up main output folder (named after the Excel file) ----
    [excelFolder, excelBaseName, ~] = fileparts(excelPath);
    outputRoot = fullfile(excelFolder, excelBaseName);
    if ~isfolder(outputRoot)
        mkdir(outputRoot);
    end
    fprintf('Output folder:\n  %s\n\n', outputRoot);

    % ---- Parse the scoring table (user-provided logic) ----
    T = readtable(excelPath, 'Sheet', 2, 'ReadVariableNames', false);
    if strcmp(T{1,1}{1}, 'file name')
        T = T(2:end, :);
    end
    nRec = height(T);

    for r = 1:nRec
        tdmsPath = T{r,1}{1};
        segFlat = table2array(T(r, 2:end));
        segFlat = segFlat(~isnan(segFlat));
        segments = reshape(segFlat, 2, [])';

        if wakeTF
            segments = [segments segments(:,2) + 30]; %#ok<AGROW> % preserved from original code
        end

        % ---- Convert TDMS path -> video path ----
        [folder, name, ~] = fileparts(tdmsPath);
        videoPath = fullfile(folder, [name '.mp4']); % <-- adjust extension if needed

        if ~isfile(videoPath)
            warning('analyzeRemEvents:videoNotFound', 'Skipping (video not found): %s', videoPath);
            continue;
        end

        nSegments = size(segments, 1);
        fprintf('Recording %d/%d: %s (%d REM segments)\n', r, nRec, name, nSegments);

        for s = 1:nSegments
            startSec = segments(s, 1);
            stopSec = segments(s, 2);
            if stopSec <= startSec
                warning('analyzeRemEvents:badSegment', ...
                    '  Skipping segment %d (stop <= start: %.2f, %.2f)', s, startSec, stopSec);
                continue;
            end

            fprintf('  Segment %d/%d: %.2f-%.2f s\n', s, nSegments, startSec, stopSec);
            processRemSegment(videoPath, name, startSec, stopSec, outputRoot, opts);
        end
    end

    fprintf('\nAll recordings processed.\n');
end

function processRemSegment(videoPath, mp4Name, startSec, stopSec, outputRoot, opts)
%PROCESSREMSEGMENT  Stream through one REM segment of one video, computing
%   overall std-dev / motion-energy maps (and optional windowed maps),
%   writing everything as 32-bit float TIFFs.

    startTag = formatSecTag(startSec);
    stopTag = formatSecTag(stopSec);

    vr = VideoReader(videoPath);
    fps = vr.FrameRate;
    castFn = str2func(opts.Precision);

    startFrame = max(1, round(startSec * fps) + 1);
    nFramesToRead = round((stopSec - startSec) * fps);
    if nFramesToRead < 2
        warning('analyzeRemEvents:tooShort', ...
            '    Segment too short after frame rounding, skipping: %s %.2f-%.2f', ...
            mp4Name, startSec, stopSec);
        return;
    end

    % Seek to the segment start
    vr.CurrentTime = startSec;

    % ---- Set up windowed output subfolder if requested ----
    useWindows = ~isempty(opts.WindowSize);
    if useWindows
        winSubfolder = fullfile(outputRoot, sprintf('%s_%s_%s', mp4Name, startTag, stopTag),'stdDev');
        if ~isfolder(winSubfolder)
            mkdir(winSubfolder);
        end
        winSubfolder = fullfile(outputRoot, sprintf('%s_%s_%s', mp4Name, startTag, stopTag),'diffStdDev');
        if ~isfolder(winSubfolder)
            mkdir(winSubfolder);
        end
        winSubfolder = fullfile(outputRoot, sprintf('%s_%s_%s', mp4Name, startTag, stopTag));
    end

    % ---- Accumulator init ----
    count = 0;   meanImg = [];   M2 = [];
    countD = 0;  meanD = [];     M2D = [];
    wCount = 0;  wMean = [];     wM2 = [];
    wCountD = 0; wMeanD = [];    wM2D = [];
    windowIdx = 0;
    meWindowIdx = 0;
    prevFrame = [];
    localFrameIdx = 0;
    localDiffIdx = 0;

    for i = 1:nFramesToRead
        if ~hasFrame(vr)
            break;
        end
        frame = readFrame(vr);
        if size(frame, 3) == 3
            frameGray = rgb2gray(frame);
        else
            frameGray = frame;
        end
        frameGray = castFn(im2double(frameGray));

        if isempty(meanImg)
            [h, w] = size(frameGray);
            meanImg = zeros(h, w, opts.Precision);
            M2      = zeros(h, w, opts.Precision);
            meanD   = zeros(h, w, opts.Precision);
            M2D     = zeros(h, w, opts.Precision);
            if useWindows
                wMean  = zeros(h, w, opts.Precision);
                wM2    = zeros(h, w, opts.Precision);
                wMeanD = zeros(h, w, opts.Precision);
                wM2D   = zeros(h, w, opts.Precision);
            end
        end

        localFrameIdx = localFrameIdx + 1;
        absFrameIdx = startFrame + localFrameIdx - 1;

        % ---- Welford update: overall std-dev accumulator ----
        count = count + 1;
        delta = frameGray - meanImg;
        meanImg = meanImg + delta / count;
        M2 = M2 + delta .* (frameGray - meanImg);

        % ---- Windowed std-dev accumulator ----
        if useWindows
            wCount = wCount + 1;
            wDelta = frameGray - wMean;
            wMean = wMean + wDelta / wCount;
            wM2 = wM2 + wDelta .* (frameGray - wMean);

            if wCount == opts.WindowSize
                windowIdx = windowIdx + 1;
                windowStd = sqrt(wM2 / (wCount - 1));
                outPath = fullfile(winSubfolder,'stdDev', ...
                    sprintf('stdDev_%04d_frames_%d-%d.tif', windowIdx, ...
                    absFrameIdx - opts.WindowSize + 1, absFrameIdx));
                writeFloatTiff(outPath, windowStd);
                wCount = 0;
                wMean(:) = 0;
                wM2(:) = 0;
            end
        end

        % ---- Frame-differencing / motion-energy accumulator ----
        if ~isempty(prevFrame)
            d = abs(frameGray - prevFrame);
            localDiffIdx = localDiffIdx + 1;

            countD = countD + 1;
            deltaD = d - meanD;
            meanD = meanD + deltaD / countD;
            M2D = M2D + deltaD .* (d - meanD);

            if useWindows
                wCountD = wCountD + 1;
                wDeltaD = d - wMeanD;
                wMeanD = wMeanD + wDeltaD / wCountD;
                wM2D = wM2D + wDeltaD .* (d - wMeanD);

                if wCountD == opts.WindowSize
                    meWindowIdx = meWindowIdx + 1;
                    windowDiffStd = sqrt(wM2D / (wCountD - 1));
                    outPath = fullfile(winSubfolder,'diffStdDev', ...
                        sprintf('diffStdDev_%04d_diffs_%d-%d.tif', meWindowIdx, ...
                        absFrameIdx - opts.WindowSize + 1, absFrameIdx));
                    writeFloatTiff(outPath, windowDiffStd);
                    wCountD = 0;
                    wMeanD(:) = 0;
                    wM2D(:) = 0;
                end
            end
        end
        prevFrame = frameGray;

        if mod(localFrameIdx, opts.ProgressEvery) == 0
            fprintf('    processed %d/%d frames\n', localFrameIdx, nFramesToRead);
        end
    end

    if count < 2
        warning('analyzeRemEvents:tooFewFrames', ...
            '    Fewer than 2 frames actually read for %s %.2f-%.2f, skipping overall maps.', ...
            mp4Name, startSec, stopSec);
        return;
    end

    stdMap = sqrt(double(M2) / double(count - 1));
    motionEnergyMap = sqrt(double(M2D) / double(countD - 1));

    if ~isfolder(fullfile(outputRoot,'stdDev'))
        mkdir(fullfile(outputRoot,'stdDev'));
    end
    if ~isfolder(fullfile(outputRoot,'diffStdDev'))
        mkdir(fullfile(outputRoot,'diffStdDev'));
    end
    stdOutPath = fullfile(outputRoot,'stdDev', sprintf('stdDev_%s_%s_%s.tif', mp4Name, startTag, stopTag));
    diffOutPath = fullfile(outputRoot,'diffStdDev', sprintf('diffStdDev_%s_%s_%s.tif', mp4Name, startTag, stopTag));
    writeFloatTiff(stdOutPath, stdMap);
    writeFloatTiff(diffOutPath, motionEnergyMap);

    fprintf('    wrote overall stdDev + diffStdDev (%d frames)\n', count);
    if useWindows
        fprintf('    wrote %d windowed stdDev + %d windowed diffStdDev TIFFs to:\n      %s\n', ...
            windowIdx, meWindowIdx, winSubfolder);
    end
end

function tag = formatSecTag(secVal)
%FORMATSECTAG  Turn a seconds value into a filename-safe string, e.g.
%   123.5 -> '123p5', 90 -> '90'.
    tag = strrep(sprintf('%g', secVal), '.', 'p');
    tag = strrep(tag, '-', 'neg');
end

function writeFloatTiff(outPath, img)
%WRITEFLOATTIFF  Write a single HxW image as a 32-bit float TIFF, openable
%   directly in ImageJ/Fiji with real (non-rescaled) pixel values.
    img = single(img);
    [h, w] = size(img);

    t = Tiff(outPath, 'w');
    tagstruct = struct();
    tagstruct.ImageLength = h;
    tagstruct.ImageWidth = w;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 32;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';
    t.setTag(tagstruct);
    t.write(img);
    t.close();
end