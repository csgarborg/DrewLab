function results = analyzeMouseMotion(videoPath, varargin)
%ANALYZEMOUSEMOTION  Quantify and visualize motion in a head-fixed mouse video.
%   Memory-safe, single-pass streaming version for large videos (does NOT
%   load the whole frame stack into memory -- uses Welford's online
%   algorithm to accumulate per-pixel mean/variance frame by frame).
%
%   results = analyzeMouseMotion(videoPath) reads the video at videoPath,
%   converts every frame to grayscale, and computes:
%       1. A per-pixel standard-deviation projection across the whole
%          video (bright = high temporal variance = lots of motion there)
%       2. A frame-differencing motion trace (mean |frame(t)-frame(t-1)|
%          per frame) showing WHEN motion happens over time
%       3. A per-pixel motion-energy map (std of frame differences)
%
%   Name-Value options:
%       'WindowSize'      - if set (e.g. 300), also computes:
%                              - a std-dev projection per non-overlapping
%                                window of raw frames
%                              - a motion-energy map (std of frame diffs)
%                                per non-overlapping window of diffs
%                            Both are WRITTEN TO DISK as each window
%                            finishes -- as individual 32-bit float TIFF
%                            files, directly openable in ImageJ/Fiji, not
%                            held in RAM.
%                            Default: [] (off)
%       'WindowOutputDir' - folder to save windowed maps into. Default:
%                            a folder named '<video>_windows' in pwd.
%                            Std-dev windows are saved as
%                              stdDev_####_frames_X-Y.tif
%                            Motion-energy (frame-diff std) windows are
%                            saved as
%                              diffStdDev_####_diffs_X-Y.tif
%       'Show'            - true/false, whether to plot results (default: true)
%       'ProgressEvery'   - print progress every N frames (default: 5000)
%       'Precision'       - 'single' or 'double' for accumulators
%                            (default: 'single')
%
%   OUTPUT (struct):
%       results.stdMap           - HxW std-dev projection (full video)
%       results.motionTrace      - 1x(N-1) frame differencing trace
%       results.motionEnergyMap  - HxW std of frame differences (full video)
%       results.fps              - detected frame rate
%       results.nFrames          - number of frames processed
%       results.windowOutputDir  - folder containing per-window .tif files
%                                  (only set if WindowSize was used)
%       results.nWindows             - number of stdDev window TIFFs written
%       results.nMotionEnergyWindows - number of diffStdDev window TIFFs written
%
%   Example:
%       r = analyzeMouseMotion('mouse_face.mp4', 'WindowSize', 300);

    p = inputParser;
    addParameter(p, 'WindowSize', []);
    addParameter(p, 'WindowOutputDir', '');
    addParameter(p, 'Show', true);
    addParameter(p, 'ProgressEvery', 5000);
    addParameter(p, 'Precision', 'single');
    parse(p, varargin{:});
    opts = p.Results;

    if ~isfile(videoPath)
        error('analyzeMouseMotion:fileNotFound', 'Cannot find file: %s', videoPath);
    end

    castFn = str2func(opts.Precision);

    % ---- Set up video reader ----
    fprintf('Opening video: %s\n', videoPath);
    vr = VideoReader(videoPath);
    fps = vr.FrameRate;
    nFramesEst = floor(vr.Duration * vr.FrameRate);
    fprintf('  ~%d frames at %.2f fps (%.1f s)\n', nFramesEst, fps, vr.Duration);

    % ---- Set up windowed-map output folder if requested ----
    useWindows = ~isempty(opts.WindowSize);
    if useWindows
        if isempty(opts.WindowOutputDir)
            [~, vname] = fileparts(videoPath);
            opts.WindowOutputDir = fullfile(pwd, [vname '_windows']);
        end
        if ~isfolder(opts.WindowOutputDir)
            mkdir(opts.WindowOutputDir);
        end
        fprintf('  Windowed stdDev and diffStdDev TIFFs (window=%d) will be saved to:\n    %s\n', ...
            opts.WindowSize, opts.WindowOutputDir);
    end

    % ---- Preallocate motion trace (grow if underestimated) ----
    traceCap = max(nFramesEst, 100);
    motionTrace = zeros(1, traceCap, 'double');

    % ---- Accumulator init (deferred until first frame gives us H,W) ----
    count = 0;   meanImg = [];   M2 = [];
    countD = 0;  meanD = [];     M2D = [];
    wCount = 0;  wMean = [];     wM2 = [];
    wCountD = 0; wMeanD = [];    wM2D = [];
    windowIdx = 0;
    meWindowIdx = 0;
    prevFrame = [];
    frameIdx = 0;
    diffIdx = 0;

    tic;
    while hasFrame(vr)
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

        % ---- Welford update: global std-dev accumulator ----
        frameIdx = frameIdx + 1;
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
                outPath = fullfile(opts.WindowOutputDir,'stdDev', ...
                    sprintf('stdDev_%04d_frames_%d-%d.tif', windowIdx, ...
                    frameIdx - opts.WindowSize + 1, frameIdx));
                writeFloatTiff(outPath, windowStd);
                wCount = 0;
                wMean(:) = 0;
                wM2(:) = 0;
            end
        end

        % ---- Frame-differencing motion trace + motion energy accumulator ----
        if ~isempty(prevFrame)
            d = abs(frameGray - prevFrame);
            diffIdx = diffIdx + 1;
            if diffIdx > numel(motionTrace)
                motionTrace(end + traceCap) = 0;
            end
            motionTrace(diffIdx) = mean(d(:));

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
                    outPath = fullfile(opts.WindowOutputDir,'diffStdDev', ...
                        sprintf('diffStdDev_%04d_diffs_%d-%d.tif', meWindowIdx, ...
                        diffIdx - opts.WindowSize + 1, diffIdx));
                    writeFloatTiff(outPath, windowDiffStd);
                    wCountD = 0;
                    wMeanD(:) = 0;
                    wM2D(:) = 0;
                end
            end
        end
        prevFrame = frameGray;

        if mod(frameIdx, opts.ProgressEvery) == 0
            fprintf('  processed %d frames (%.1f s elapsed)\n', frameIdx, toc);
        end
    end

    motionTrace = motionTrace(1:diffIdx);

    if count < 2
        error('analyzeMouseMotion:tooFewFrames', 'Need at least 2 frames to compute motion.');
    end

    stdMap = sqrt(double(M2) / double(count - 1));
    motionEnergyMap = sqrt(double(M2D) / double(countD - 1));

    fprintf('Done. Processed %d frames in %.1f s.\n', frameIdx, toc);

    results = struct();
    results.stdMap = stdMap;
    results.motionTrace = motionTrace;
    results.motionEnergyMap = motionEnergyMap;
    results.fps = fps;
    results.nFrames = frameIdx;
    if useWindows
        results.windowOutputDir = opts.WindowOutputDir;
        results.nWindows = windowIdx;
        results.nMotionEnergyWindows = meWindowIdx;
    end

    if opts.Show
        plotResults(results);
    end

    save(fullfile(opts.WindowOutputDir,[vname '.mat']),'results')
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

function plotResults(results)
    figure('Name', 'Mouse Motion Analysis', 'Color', 'w', 'Position', [100 100 1100 800]);

    subplot(2,2,1);
    imagesc(results.stdMap);
    axis image off; colormap(gca, hot); colorbar;
    title('Std-Dev Projection (whole video)');

    subplot(2,2,2);
    imagesc(results.motionEnergyMap);
    axis image off; colormap(gca, hot); colorbar;
    title('Motion Energy Map (std of frame diffs)');

    subplot(2,2,[3 4]);
    t = (1:length(results.motionTrace)) / results.fps;
    plot(t, results.motionTrace, 'k', 'LineWidth', 1);
    xlabel('Time (s)'); ylabel('Mean |\Deltapixel|');
    title('Frame-Differencing Motion Trace');
    xlim([0 t(end)]);
    grid on;

    if isfield(results, 'windowOutputDir')
        fprintf('Windowed TIFFs saved individually in:\n  %s\n', results.windowOutputDir);
        fprintf('  Std-dev windows:       stdDev_####_frames_X-Y.tif\n');
        fprintf('  Motion-energy windows: diffStdDev_####_diffs_X-Y.tif\n');

        dStd = dir(fullfile(results.windowOutputDir, 'stdDev_*.tif'));
        if ~isempty(dStd)
            nShow = min(16, numel(dStd));
            showIdx = round(linspace(1, numel(dStd), nShow));
            figure('Name', 'Windowed Std-Dev Preview', 'Color', 'w');
            for i = 1:nShow
                img = imread(fullfile(dStd(showIdx(i)).folder, dStd(showIdx(i)).name));
                subplot(ceil(sqrt(nShow)), ceil(sqrt(nShow)), i);
                imagesc(img); axis image off; colormap(hot);
                title(sprintf('win %d', showIdx(i)), 'FontSize', 8);
            end
            sgtitle('Sample of windowed std-dev projections over time');
        end

        dME = dir(fullfile(results.windowOutputDir, 'diffStdDev_*.tif'));
        if ~isempty(dME)
            nShow = min(16, numel(dME));
            showIdx = round(linspace(1, numel(dME), nShow));
            figure('Name', 'Windowed Motion-Energy Preview', 'Color', 'w');
            for i = 1:nShow
                img = imread(fullfile(dME(showIdx(i)).folder, dME(showIdx(i)).name));
                subplot(ceil(sqrt(nShow)), ceil(sqrt(nShow)), i);
                imagesc(img); axis image off; colormap(hot);
                title(sprintf('win %d', showIdx(i)), 'FontSize', 8);
            end
            sgtitle('Sample of windowed motion-energy (diff std) maps over time');
        end
    end
end