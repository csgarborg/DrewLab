function tdmsFilePaths = findTDMSFiles(dateList)

if ~exist('dateList','var')
    % dateList = {'251201','251202','251203','251204'};
    dateList = {'251201'};
end
dataDir = 'H:\SGNMData\';
tdmsFilePaths = {};
for n = 1:length(dateList)
    files = dir([dataDir dateList{n} '\*.tdms']);
    fileNames = {files(~[files.isdir]).name};
    filteredFiles = fileNames(~contains(lower(fileNames), 'converted'));
    for i = 1:length(filteredFiles)
        filteredFiles{i} = [dataDir dateList{n} '\' filteredFiles{i}];
    end
    tdmsFilePaths = [tdmsFilePaths,filteredFiles];
end
end