function batchAnnotateAndCompressMouseAVI(hbPreset)
% Batch annotate AVI files using matching CSVs, compress with HandBrake.
% Features:
% - Per-file waitbars
% - Cancel button
% - Throttled UI updates
% - Annotation ETA
% - Compression ETA
% - Total batch ETA

    %% Hard-coded HandBrakeCLI path
    hbExe = '"C:\Users\csg178\handbrakeCLI\HandBrakeCLI.exe"';

    %% Select folder
    folderPath = uigetdir('C:\Data\Spencer\LabVIEWTestData', 'Select folder containing AVI and CSV files');
    if folderPath == 0
        disp('No folder selected. Exiting.');
        return;
    end

    aviFiles = dir(fullfile(folderPath, '*.avi'));
    if isempty(aviFiles)
        error('No AVI files found.');
    end

    numFiles = numel(aviFiles);
    batchStart = tic;

    completedFiles = 0;
    completedAnnotTimes = [];
    completedCompTimes = [];

    for k = 1:numFiles
        aviName = aviFiles(k).name;
        [~, baseName, ~] = fileparts(aviName);

        aviPath = fullfile(folderPath, aviName);
        csvPath = fullfile(folderPath, [baseName '.csv']);
        tempAvi = fullfile(folderPath, [baseName '_temp.avi']);
        mp4Out  = fullfile(folderPath, [baseName '.mp4']);

        if ~isfile(csvPath)
            warning('CSV not found for %s — skipping.', aviName);
            continue;
        end

        %% Read CSV
        opts = detectImportOptions(csvPath, ...
            'Delimiter', ',', ...
            'VariableNamingRule', 'preserve');
        opts = setvartype(opts, 'char');   % force everything to text
        tbl = readtable(csvPath, opts);
        csvText = string(tbl{:,1});        % full, unmodified lines

        %% Video setup
        vr = VideoReader(aviPath);
        vw = VideoWriter(tempAvi, 'Motion JPEG AVI');
        vw.FrameRate = vr.FrameRate;
        open(vw);

        totalFrames = floor(vr.Duration * vr.FrameRate);
        updateEvery = max(1, round(totalFrames / 100)); % ~100 updates/file
        frameIdx = 1;

        %% Waitbar (per file)
        hWait = waitbar(0, ...
            sprintf('Annotating %s', aviName), ...
            'Name', sprintf('File %d of %d', k, numFiles), ...
            'CreateCancelBtn', 'setappdata(gcbf,''cancel'',true)');
        setappdata(hWait, 'cancel', false);
        set(findall(hWait,'Type','text'), 'Interpreter','none');
        fig = ancestor(hWait, 'figure');

        % Make the window taller so text isn't clipped
        fig.Position(4) = fig.Position(4) + 40;   % add height (pixels)

        tAnnot = tic;
        cancelFlag = false;

        %% ---------- Annotation loop ----------
        while hasFrame(vr)
            if getappdata(hWait, 'cancel')
                cancelFlag = true;
                break;
            end

            frame = readFrame(vr);

            if frameIdx <= numel(csvText)
                txt = csvText{frameIdx};
            else
                txt = csvText{end};
            end

            frame = insertText(frame, [1 1], txt, ...
                'FontSize', 12, ...
                'TextColor', 'white', ...
                'BoxColor', 'black', ...
                'BoxOpacity', 0.8, ...
                'AnchorPoint', 'LeftTop');

            writeVideo(vw, frame);

            % Throttled UI update
            if mod(frameIdx, updateEvery) == 0 || frameIdx == totalFrames
                progress = frameIdx / totalFrames;
                elapsedAnnot = toc(tAnnot);
                etaAnnot = elapsedAnnot / progress - elapsedAnnot;

                % Batch ETA estimate
                if ~isempty(completedAnnotTimes)
                    avgFileTime = mean(completedAnnotTimes + completedCompTimes);
                    remainingFiles = numFiles - completedFiles - 1;
                    etaBatch = remainingFiles * avgFileTime + etaAnnot;
                else
                    etaBatch = NaN;
                end

                waitbar(progress, hWait, sprintf( ...
                    '%s\nAnnotating: %d / %d frames\nAnnotation ETA: %s\nBatch ETA: %s', ...
                    aviName, frameIdx, totalFrames, ...
                    formatETA(etaAnnot), formatETA(etaBatch)));
            end

            frameIdx = frameIdx + 1;
        end

        close(vw);

        if cancelFlag
            delete(hWait);
            if isfile(tempAvi), delete(tempAvi); end
            disp('Processing cancelled by user.');
            return;
        end

        annotTime = toc(tAnnot);

        %% ---------- Compression ----------
        tComp = tic;

        % Search for preset
        presetDir = 'C:\Data\Spencer\LabVIEWTestData\handbrakeEncodeFiles'; % Folder containing HandBrake preset JSON files

        % Look for JSON with same base name
        jsonPath = fullfile(presetDir, [hbPreset '.json']);
        if isfile(jsonPath)
            fprintf('Found preset JSON:\n%s\n', jsonPath);
        else
            warning('No matching preset JSON found for %s', baseName);
            jsonPath = '';
        end

        if nargin < 1 || isempty(hbPreset)
            disp('Preparing to encode default preset')
            hbCmd = sprintf([ ...
                '%s -i "%s" -o "%s" ' ...
                '-e x264 -q 20 --encoder-preset slow ' ...
                '--encoder-profile high --encoder-level 4.1 ' ...
                '--optimize'], ...
                hbExe, tempAvi, mp4Out);
        elseif ~isempty(jsonPath)
            disp('Preparing to encode with custom preset')
            hbCmd = sprintf( ...
                '%s --preset-import-file "%s" --preset "%s" -i "%s" -o "%s"', ...
                hbExe, ...
                jsonPath, ...
                hbPreset, ...
                tempAvi, mp4Out);
        else
            disp('Preparing to encode selected standard preset')
            hbCmd = sprintf( ...
                '%s -i "%s" -o "%s" --preset "%s"', ...
                hbExe, tempAvi, mp4Out, hbPreset);
        end

        % Compression ETA estimate
        if ~isempty(completedCompTimes)
            etaComp = mean(completedCompTimes);
        else
            etaComp = NaN;
        end

        waitbar(1, hWait, sprintf( ...
            '%s\nCompressing...\nEst. Compression ETA: %s', ...
            aviName, formatETA(etaComp)));

        status = system(hbCmd);
        if status ~= 0
            warning('HandBrake failed for %s', aviName);
            delete(hWait);
            continue;
        end

        compTime = toc(tComp);

        %% Cleanup + stats
        if isfile(tempAvi), delete(tempAvi); end
        delete(hWait);

        completedFiles = completedFiles + 1;
        completedAnnotTimes(end+1) = annotTime;
        completedCompTimes(end+1)  = compTime;

        fprintf('Finished %s (%d/%d)\n', aviName, completedFiles, numFiles);
    end

    fprintf('\nAll files processed.\n');
end

%% ---------- Helper ----------
function str = formatETA(seconds)
    if isnan(seconds) || seconds < 0 || isinf(seconds)
        str = 'calculating...';
        return;
    end
    hrs = floor(seconds / 3600);
    mins = floor(mod(seconds, 3600) / 60);
    secs = floor(mod(seconds, 60));
    if hrs > 0
        str = sprintf('%dh %dm %ds', hrs, mins, secs);
    elseif mins > 0
        str = sprintf('%dm %ds', mins, secs);
    else
        str = sprintf('%ds', secs);
    end
end
