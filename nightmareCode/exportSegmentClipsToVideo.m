function exportSegmentClipsToVideo(excelPath)

%% Load Excel
T = readtable(excelPath,'Sheet',2,'ReadVariableNames',false);

nRec = height(T);

% Output video path
[excelFolder,excelName,~] = fileparts(excelPath);
outputPath = fullfile(excelFolder,[excelName '.mp4']);

vWriter = VideoWriter(outputPath,'MPEG-4'); % H264
vWriter.FrameRate = 60; % will adjust dynamically if needed
open(vWriter)

fprintf('Writing output video: %s\n',outputPath)

for r = 1:nRec
    
    tdmsPath = T{r,1}{1};
    [folder,name,~] = fileparts(tdmsPath);
    
    videoPath = fullfile(folder,[name '.mp4']);
    
    if ~exist(videoPath,'file')
        warning('Video not found: %s',videoPath)
        continue
    end
    
    fprintf('Processing: %s\n',name)
    
    vReader = VideoReader(videoPath);
    
    % Extract segments
    segments = table2array(T(r,2:end));
    segments = segments(~isnan(segments));
    
    if mod(length(segments),2) ~= 0
        warning('Skipping row %d (uneven segments)',r)
        continue
    end
    
    segments = reshape(segments,2,[])';
    
    % Update writer framerate to match input (once)
    % vWriter.FrameRate = vReader.FrameRate;
    
    for s = 1:size(segments,1)
        
        startTime = segments(s,1);
        stopTime  = segments(s,2);
        
        vReader.CurrentTime = startTime;
        
        while vReader.CurrentTime < stopTime && hasFrame(vReader)
            
            frame = readFrame(vReader);
            
            % Overlay text
            label = sprintf('%s | %.2f - %.2f s',name,startTime,stopTime);
            
            frame = insertText(frame,[10 10],label, ...
                'FontSize',18,'BoxColor','black','TextColor','white','BoxOpacity',0.6);
            
            writeVideo(vWriter,frame);
            
        end
        
        % Add 5 black frames between clips
        blackFrame = zeros(size(frame),'uint8');
        
        for k = 1:5
            writeVideo(vWriter,blackFrame);
        end
        
    end
    
end

close(vWriter)

fprintf('Done!\n')

end