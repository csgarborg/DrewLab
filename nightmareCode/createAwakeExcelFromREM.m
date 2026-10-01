function createAwakeExcelFromREM(excelPath)
% createAwakeExcelFromREM
%
% Reads an Excel file containing .mp4 paths in column 1 of Sheet 1.
%
% For each unique date folder (yymmdd), the function:
%   1. Finds the first .mp4 file for that date.
%   2. Uses that file as the Awake file for that day.
%
% The output Excel file:
%   - Uses the same filename as the input Excel file.
%   - Replaces 'REM' with 'Awake' in the Excel filename.
%   - Sheet 1 contains the selected Awake .mp4 paths in column 1.
%   - Sheet 2 contains:
%       file name, start1_s, stop1_s, ..., start5_s, stop5_s
%     followed by one row for each Awake file.
%
% INPUT:
%   excelPath - Full path to the input Excel file.
%
% EXAMPLE:
%   createAwakeExcelFromREM('D:\Data\SGNM017_REM.xlsx');

    %% Validate input
    if ~isfile(excelPath)
        error('Excel file does not exist:\n%s', excelPath);
    end

    %% Read Sheet 1
    data = readcell(excelPath, 'Sheet', 1);

    if isempty(data)
        error('Sheet 1 is empty.');
    end

    % First column contains MP4 paths
    filePaths = data(:,1);

    % Remove empty cells
    validRows = ~cellfun(@(x) isempty(x) || ...
        (isstring(x) && strlength(x) == 0) || ...
        (ischar(x) && isempty(strtrim(x))), filePaths);

    filePaths = filePaths(validRows);

    % Convert everything to character vectors
    filePaths = cellfun(@char, filePaths, 'UniformOutput', false);

    %% Extract date folders
    %
    % Expected structure:
    %
    %   ...\yymmdd\yymmdd_hh_mm_ssdd.mp4
    %
    % Example:
    %
    %   D:\MouseData\260929\260929_10_35_1234.mp4

    dateStrings = cell(size(filePaths));

    for i = 1:numel(filePaths)

        [parentFolder, ~, ~] = fileparts(filePaths{i});
        [~, dateFolder, ~] = fileparts(parentFolder);

        % Make sure the folder looks like yymmdd
        if isempty(regexp(dateFolder, '^\d{6}$', 'once'))
            error(['Could not identify a yymmdd date folder for:\n' ...
                   '%s\n\nExpected the parent folder to be something like 260929.'], ...
                   filePaths{i});
        end

        dateStrings{i} = dateFolder;
    end

    %% Find unique dates
    uniqueDates = unique(dateStrings, 'stable');

    awakeFiles = cell(numel(uniqueDates),1);

    %% Find the first MP4 for each date
    for i = 1:numel(uniqueDates)

        thisDate = uniqueDates{i};

        % Find all entries belonging to this date
        dateIdx = strcmp(dateStrings, thisDate);

        dateFiles = filePaths(dateIdx);

        % Sort alphabetically by filename.
        % Because filenames begin with yymmdd_hh_mm_ssdd,
        % this also sorts chronologically.
        [~, sortIdx] = sort(lower(dateFiles));
        dateFiles = dateFiles(sortIdx);

        % Select the first file for this date
        awakeFiles{i} = dateFiles{1};
    end

    %% Create output Excel filename
    %
    % Example:
    %   SGNM017_REM.xlsx
    %
    % becomes:
    %   SGNM017_Awake.xlsx

    [inputFolder, inputName, inputExt] = fileparts(excelPath);

    outputName = strrep(inputName, 'REM', 'Awake');

    % If REM wasn't present, append Awake instead
    if strcmp(outputName, inputName)
        outputName = [inputName '_Awake'];
    end

    outputExcelPath = fullfile(inputFolder, ...
        [outputName inputExt]);

    %% Sheet 1
    %
    % First column = first MP4 file for each date

    writecell(awakeFiles, outputExcelPath, 'Sheet', 1, 'Range', 'A1');

    %% Sheet 2 header
    headers = { ...
        'file name', ...
        'start1_s', 'stop1_s', ...
        'start2_s', 'stop2_s', ...
        'start3_s', 'stop3_s', ...
        'start4_s', 'stop4_s', ...
        'start5_s', 'stop5_s'};

    writecell(headers, outputExcelPath, 'Sheet', 2, 'Range', 'A1');

    %% Create Sheet 2 data
    %
    % Column 1:
    %   MP4 path -> TDMS path
    %
    % Column 2:
    %   0
    %
    % Column 3:
    %   900
    %
    % Remaining start/stop columns remain empty.

    sheet2Data = cell(numel(awakeFiles), 11);

    for i = 1:numel(awakeFiles)

        % Replace .mp4 with .tdms
        [folder, fileName, ~] = fileparts(awakeFiles{i});
        tdmsPath = fullfile(folder, [fileName '.tdms']);

        sheet2Data{i,1} = tdmsPath;
        sheet2Data{i,2} = 0;
        sheet2Data{i,3} = 900;

    end

    writecell(sheet2Data, outputExcelPath, ...
        'Sheet', 2, 'Range', 'A2');

    %% Display result
    fprintf('\nCreated Awake Excel file:\n%s\n', outputExcelPath);
    fprintf('Unique dates found: %d\n', numel(uniqueDates));
    fprintf('Awake files written: %d\n', numel(awakeFiles));

end