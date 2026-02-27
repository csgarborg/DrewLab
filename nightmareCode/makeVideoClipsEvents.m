function makeVideoClipsEvents(dateList,puffLength,videoOutPath)

tdmsFilePaths = findTDMSFiles(dateList);

videoIdx = 1;
for n = 1:length(tdmsFilePaths)
    if ~exist(strrep(tdmsFilePaths{n},'.tdms','_puffDataStruct.mat'))
        continue
    end
    tdmsDataStruct = openTDMS(tdmsFilePaths{n});
    secBeforePuff = 2;
    secAfterPuff = str2double(tdmsDataStruct.PuffDuration_s_) + 8;
    if puffLength ~= str2double(tdmsDataStruct.PuffDuration_s_)
        continue
    end
    digitalPuffIdx = find(diff(tdmsDataStruct.Digital_Data.Puff > 0) == 1) + 1;
    if isempty(digitalPuffIdx)
        continue
    end

    for i = 1:length(digitalPuffIdx)
        % Create VideoWriter object
        if ~exist(videoOutPath,'dir')
            mkdir(videoOutPath)
        end
        vw = VideoWriter([videoOutPath '\' num2str(videoIdx) '.mp4'], 'MPEG-4');
        tdmsDataStruct = openTDMS(tdmsFilePaths{1});
        vw.FrameRate = str2double(tdmsDataStruct.CameraFrameratePerSecond);   % Match original frame rate
        open(vw);
        videoIdx = videoIdx + 1;

        digitalFramePuffTimeS = digitalPuffIdx(i) / str2double(tdmsDataStruct.CameraFrameratePerSecond);
        startTime = digitalFramePuffTimeS - secBeforePuff;
        stopTime = digitalFramePuffTimeS + secAfterPuff;
        videoFilePath = strrep(tdmsFilePaths{n},'.tdms','.mp4');
        videoImages = readMP4Video(false,startTime,stopTime,videoFilePath);
        for k = 1:length(videoImages)
            writeVideo(vw, videoImages{k});
        end
        close(vw)
    end
end
end