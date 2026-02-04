function drawROIBoxes(txtFile,imageFile)

if ~exist('imageFile','var')
    txtFile   = 'H:\SGNMData\260114\260114_SGNM005_ROI.txt';     % your text file
    imageFile = 'C:\Users\csg178\Pictures\vlcsnap-2026-01-21-21h39m45s621.png';     % image to draw on
end


% Read image
img = imread(imageFile);

figure;
imshow(img);
hold on;

% Open and read file
fid = fopen(txtFile, 'r');

% Read and discard first line (header)
fgetl(fid);

% Loop through remaining lines
while ~feof(fid)
    line = fgetl(fid);
    if isempty(line)
        continue
    end

    % Split by commas
    parts = strsplit(line, ',');

    % Parse data
    label = strtrim(parts{1});
    x1 = str2double(parts{2});
    y1 = str2double(parts{3});
    x2 = str2double(parts{4});
    y2 = str2double(parts{5});

    % Convert to rectangle format
    width  = x2 - x1;
    height = y2 - y1;

    % Draw rectangle
    rectangle('Position', [x1, y1, width, height], ...
              'EdgeColor', 'r', ...
              'LineWidth', 2);

    % Add label
    text(x1, y1 - 20, label, ...
        'Color', 'r', ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Interpreter', 'none');
end

fclose(fid);
hold off;