function computePixelDiffVideo(inputPath)
% COMPUTEPIXELDIFFVIDEO Computes absolute pixel difference between consecutive
% frames of a video and saves the result as a new MP4.
%
% Usage:
%   computePixelDiffVideo('C:\path\to\video.mp4')

    % --- Build output path ---
    [folder, name, ~] = fileparts(inputPath);
    outputPath = fullfile(folder, [name '_PixelDiff.mp4']);

    % --- Open input video ---
    vReader = VideoReader(inputPath);
    fps     = vReader.FrameRate;
    
    % --- Open output video ---
    vWriter = VideoWriter(outputPath, 'MPEG-4');
    vWriter.FrameRate = fps;
    open(vWriter);

    % --- Read first frame ---
    prevFrame = im2double(rgb2gray(readFrame(vReader)));

    % --- Process frames ---
    frameCount = 0;
    while hasFrame(vReader)
        currFrame = im2double(rgb2gray(readFrame(vReader)));

        % Absolute difference, scaled to [0, 1] for visibility
        diffFrame = abs(currFrame - prevFrame);

        % Optional: enhance contrast so subtle changes are visible
        % diffFrame = imadjust(diffFrame, [0, 0.1], [0, 1]);

        writeVideo(vWriter, im2uint8(diffFrame));

        prevFrame  = currFrame;
        frameCount = frameCount + 1;

        if mod(frameCount, 100) == 0
            fprintf('Processed %d frames...\n', frameCount);
        end
    end

    close(vWriter);
    fprintf('Done. %d difference frames written to:\n  %s\n', frameCount, outputPath);
end