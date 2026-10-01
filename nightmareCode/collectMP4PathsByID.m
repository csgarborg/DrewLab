function collectMP4PathsByID(inputFolder, idString)
% collectMP4PathsByID
%
% Searches each folder directly inside inputFolder for a text file whose
% filename contains idString. For every matching folder, all .mp4 files
% in that folder are collected and written to an Excel file.
%
% INPUTS:
%   inputFolder - Folder containing the folders to search
%   idString    - String to search for in text filenames
%                 Example: 'SGNM017'
%
% EXAMPLE:
%   collectMP4PathsByID('D:\MouseData', 'SGNM017');

    %% Validate inputs
    if ~isfolder(inputFolder)
        error('Input folder does not exist: %s', inputFolder);
    end

    if ~(ischar(idString) || isstring(idString))
        error('idString must be a character vector or string.');
    end

    idString = char(idString);

    %% Get folders inside input folder
    folderInfo = dir(inputFolder);
    folderInfo = folderInfo([folderInfo.isdir]);

    % Remove "." and ".."
    folderInfo = folderInfo(~ismember({folderInfo.name}, {'.', '..'}));

    %% Initialize output
    mp4Paths = {};

    %% Search each folder
    for i = 1:numel(folderInfo)

        currentFolder = fullfile(inputFolder, folderInfo(i).name);

        % Look for text files containing the ID string
        txtFiles = dir(fullfile(currentFolder, '*.txt'));

        matchFound = false;

        for j = 1:numel(txtFiles)
            if contains(txtFiles(j).name, idString, 'IgnoreCase', true)
                matchFound = true;
                break;
            end
        end

        % If this folder does not contain a matching text file, skip it
        if ~matchFound
            continue;
        end

        fprintf('Match found: %s\n', currentFolder);

        %% Get all MP4 files in this folder
        mp4Files = dir(fullfile(currentFolder, '*.mp4'));

        for k = 1:numel(mp4Files)
            mp4Path = fullfile(currentFolder, mp4Files(k).name);

            % Add path as a new row in the cell vector
            mp4Paths{end+1,1} = mp4Path;
        end
    end

    %% Check whether anything was found
    if isempty(mp4Paths)
        warning('No matching folders or MP4 files were found for ID "%s".', idString);
        return;
    end

    %% Ask user where to save the Excel file
    defaultName = sprintf('%s_MP4_paths.xlsx', idString);

    [fileName, saveFolder] = uiputfile( ...
        {'*.xlsx', 'Excel Workbook (*.xlsx)'}, ...
        'Save MP4 path list', ...
        defaultName);

    % User cancelled
    if isequal(fileName, 0) || isequal(saveFolder, 0)
        fprintf('Excel file save cancelled.\n');
        return;
    end

    outputFile = fullfile(saveFolder, fileName);

    %% Write paths to Excel
    writecell(mp4Paths, outputFile);

    fprintf('\nFound %d MP4 files.\n', numel(mp4Paths));
    fprintf('Saved to:\n%s\n', outputFile);

end