function batchAnnotateAndCompressMouseAVI(hbPreset)
% batchAnnotateAndCompressAVI
% Selects a folder, annotates all AVI files using matching CSV files,
% compresses with HandBrake, and deletes temp AVIs.
%
% CSV must have same base filename as AVI.
% Example: video01.avi <-> video01.csv

    %% Select folder
    folderPath = uigetdir(pwd, 'Select folder containing AVI and CSV files');
    if folderPath == 0
        disp('No folder selected. Exiting.');
        return;
    end

    %% Find AVI files
    aviFiles = dir(fullfile(folderPath, '*.avi'));
    if isempty(aviFiles)
        error('No AVI files found in selected folder.');
    end

    for k = 1:numel(aviFiles)
        aviName = aviFiles(k).name;
        [~, baseName, ~] = fileparts(aviName);

        aviPath  = fullfile(folderPath, aviName);
        csvPath  = fullfile(folderPath, [baseName '.csv']);
        tempAvi  = fullfile(folderPath, [baseName '_temp.avi']);
        mp4Out   = fullfile(folderPath, [baseName '.mp4']);

        fprintf('\nProcessing: %s\n', aviName);

        if ~isfile(csvPath)
            warning('CSV not found for %s — skipping.', aviName);
            continue;
        end

        %% Read CSV
        csvData = readcell(csvPath);
        csvText = cellfun(@string, csvData(:,1), 'UniformOutput', false);

        %% Video setup
        vr = VideoReader(aviPath);
        vw = VideoWriter(tempAvi, 'Motion JPEG AVI');
        vw.FrameRate = vr.FrameRate;
        open(vw);

        frameIdx = 1;

        %% Annotate frames
        while hasFrame(vr)
            frame = readFrame(vr);

            if frameIdx <= numel(csvText)
                txt = csvText{frameIdx};
            else
                txt = csvText{end};
            end

            frame = insertText(frame, [10 10], txt, ...
                'FontSize', 24, ...
                'TextColor', 'white', ...
                'BoxColor', 'black', ...
                'BoxOpacity', 0.6);

            writeVideo(vw, frame);
            frameIdx = frameIdx + 1;
        end

        close(vw);

        %% HandBrake compression
        if nargin < 1 || isempty(hbPreset)
            hbCmd = sprintf([ ...
                'HandBrakeCLI -i "%s" -o "%s" ' ...
                '-e x264 -q 20 --encoder-preset slow ' ...
                '--encoder-profile high --encoder-level 4.1 ' ...
                '--optimize'], ...
                tempAvi, mp4Out);
        else
            hbCmd = sprintf( ...
                'HandBrakeCLI -i "%s" -o "%s" --preset "%s"', ...
                tempAvi, mp4Out, hbPreset);
        end

        status = system(hbCmd);
        if status ~= 0
            warning('HandBrake failed for %s', aviName);
            continue;
        end

        %% Delete temporary AVI
        if isfile(tempAvi)
            delete(tempAvi);
        end

        fprintf('Finished: %s\n', mp4Out);
    end

    fprintf('\nAll done.\n');
end