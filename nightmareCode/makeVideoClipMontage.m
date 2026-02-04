function makeVideoClipMontage(dateList,puffLength,videoOutPath)

tdmsFilePaths = findTDMSFiles(dateList);
secBeforePuff = 2;
secAfterPuff = 14;

% Create VideoWriter object
vw = VideoWriter(videoOutPath, 'MPEG-4');
tdmsDataStruct = openTDMS(tdmsFilePaths{1});
vw.FrameRate = str2double(tdmsDataStruct.CameraFrameratePerSecond);   % Match original frame rate
open(vw);

for n = 1:length(tdmsFilePaths)
    if ~exist(strrep(tdmsFilePaths{n},'.tdms','_puffDataStruct.mat'))
        continue
    end
    tdmsDataStruct = openTDMS(tdmsFilePaths{n});
    if puffLength ~= str2double(tdmsDataStruct.PuffDuration_s_)
        continue
    end
    digitalPuffIdx = find(diff(tdmsDataStruct.Digital_Data.Puff > 0) == 1) + 1;
    for i = 1:length(digitalPuffIdx)
        digitalFramePuffTimeS = digitalPuffIdx(i) / str2double(tdmsDataStruct.CameraFrameratePerSecond);
        startTime = digitalFramePuffTimeS - secBeforePuff;
        stopTime = digitalFramePuffTimeS + secAfterPuff;
        videoFilePath = strrep(tdmsFilePaths{n},'.tdms','.mp4');
        videoImages = readMP4Video(false,startTime,stopTime,videoFilePath);
        for k = 1:length(videoImages)
            writeVideo(vw, videoImages{k});
        end
        for k = 1:5
            writeVideo(vw, videoImages{end});
        end
    end
end

close(vw)
end