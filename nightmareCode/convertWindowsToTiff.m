function convertWindowsToTiff(windowOutputDir, outputTiffPath, varargin)
%CONVERTWINDOWSTOTIFF  Combine windowed std-dev .mat files into one
%   multi-page 32-bit float TIFF stack that ImageJ/Fiji can open directly.
%
%   convertWindowsToTiff(windowOutputDir, outputTiffPath)
%   convertWindowsToTiff(..., 'Transform', 'log', 'Epsilon', 1e-4)
%
%   windowOutputDir - folder containing 'window_####_frames_X-Y.mat' files
%   outputTiffPath  - path to write, e.g. 'windowed_std_stack.tif'
%
%   Name-Value options:
%       'Transform' - 'none' (default), 'log', 'sqrt', or 'gamma'
%                     Since std-dev maps are heavily right-skewed (mostly
%                     near-zero pixels, a few large ones), a linear scale
%                     crushes small motions near black. These transforms
%                     boost small values relative to large ones so subtle
%                     motion is visible instead of washed out:
%                       'log'   - log(img + Epsilon). Strongest contrast
%                                 boost for small values; needs Epsilon
%                                 since log(0) = -Inf.
%                       'sqrt'  - sqrt(img). Gentler than log, no epsilon
%                                 needed (sqrt(0) = 0).
%                       'gamma' - img.^Gamma (Gamma < 1, e.g. 0.4-0.5).
%                                 Tunable strength between none and log,
%                                 stays bounded at 0 (no epsilon needed).
%       'Epsilon'   - offset added before log, default: 1e-4 (only used
%                     if Transform='log'). Increase if the very bottom
%                     of the range looks like uniform noise; decrease if
%                     you're not seeing enough contrast boost.
%       'Gamma'     - exponent for Transform='gamma', default: 0.5
%
%   Pixel values are written as real 32-bit floats (not rescaled to
%   8-bit), so you still get quantitative values in ImageJ -- just on
%   whichever scale you picked. Use ImageJ's Brightness/Contrast > Auto
%   after opening to get a good default display range.
%
%   Example:
%       convertWindowsToTiff('mouse_video_windows', 'windowed_std_stack_log.tif', ...
%                             'Transform', 'log', 'Epsilon', 1e-4);

    p = inputParser;
    addParameter(p, 'Transform', 'none');
    addParameter(p, 'Epsilon', 1e-4);
    addParameter(p, 'Gamma', 0.5);
    parse(p, varargin{:});
    opts = p.Results;

    if ~isfolder(windowOutputDir)
        error('convertWindowsToTiff:noDir', 'Folder not found: %s', windowOutputDir);
    end

    d = dir(fullfile(windowOutputDir, 'window_*.mat'));
    if isempty(d)
        error('convertWindowsToTiff:noFiles', 'No window_*.mat files found in %s', windowOutputDir);
    end

    names = {d.name};
    idxNums = zeros(1, numel(names));
    for i = 1:numel(names)
        tok = regexp(names{i}, 'window_(\d+)_', 'tokens');
        idxNums(i) = str2double(tok{1}{1});
    end
    [~, order] = sort(idxNums);
    d = d(order);

    fprintf('Found %d window files. Transform=%s. Writing TIFF stack to:\n  %s\n', ...
        numel(d), opts.Transform, outputTiffPath);

    if isfile(outputTiffPath)
        delete(outputTiffPath);
    end

    t = Tiff(outputTiffPath, 'w');

    for i = 1:numel(d)
        s = load(fullfile(d(i).folder, d(i).name), 'windowStd');
        img = single(s.windowStd);

        switch lower(opts.Transform)
            case 'none'
                % no-op
            case 'log'
                img = log(img + single(opts.Epsilon));
            case 'sqrt'
                img = sqrt(img);
            case 'gamma'
                img = img .^ single(opts.Gamma);
            otherwise
                error('convertWindowsToTiff:badTransform', ...
                    'Unknown Transform: %s (use none, log, sqrt, or gamma)', opts.Transform);
        end

        [h, w] = size(img);
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

        if i < numel(d)
            t.writeDirectory();
        end

        if mod(i, 50) == 0 || i == numel(d)
            fprintf('  wrote %d / %d\n', i, numel(d));
        end
    end

    t.close();
    fprintf('Done. Open %s in ImageJ/Fiji (File > Open...), then Image > Adjust > Brightness/Contrast > Auto.\n', outputTiffPath);
end