function [tdmsDataStruct] = openTDMS(tdmsPath)

try
    tdmsData = tdmsread(tdmsPath);
    tdmsInfo = tdmsinfo(tdmsPath);
    tdmsProp = tdmsreadprop(tdmsPath);
catch
    newtdmsPath = [tdmsPath(1:end-5),'_converted',tdmsPath(end-4:end)];
    if ~exist(newtdmsPath,'file')
        tdmsconvert(tdmsPath,newtdmsPath);
    end
    tdmsData = tdmsread(newtdmsPath);
    tdmsInfo = tdmsinfo(newtdmsPath);
    tdmsProp = tdmsreadprop(newtdmsPath);
end
tdmsDataStruct = table2struct(tdmsProp);
tdmsDataStruct.filePath = tdmsInfo.Path;
channelList = tdmsInfo.ChannelList{:,4};
for i = 1:size(tdmsData,2)
    headerList = tdmsData{1,i}.Properties.VariableNames;
    for n = 1:length(headerList)
        header = headerList{n};
        channelIndex = strcmp(channelList,header);
        tdmsDataStruct.(replace(replace(tdmsInfo.ChannelList{channelIndex,2}, ' ', '_'), '-', '')).(replace(replace(tdmsInfo.ChannelList{channelIndex,4},' ', '_'), '-', '')) = table2array(tdmsData{1,i}(:,n));
    end
end
