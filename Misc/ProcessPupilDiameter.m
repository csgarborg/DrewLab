function ProcessPupilDiameter(pupilCamFileList, procDataFileList)
close all
%     for a = 1:size(pupilCamFileList)
        %% Initialize data variables
%         fileID = pupilCamFileList(a,:);
fileID = 'C:\Users\csg178\Documents\MATLAB\190625_13_50_50_PupilCam.bin';
%         procDataFileID = procDataFileList(a,:);
%         load(procDataFileID)   % this loads the variable/struct ProcData which is already created in previous analysis. The pupilDiameter should be added to it and then saved to cd.

        % Do the stuff in however many subfunctions you'd like to ultimately obtain pupilDiameter
        
        %% Load pupil vid
        imageStack = openPupilCamVid(fileID, 400, 400);
        
        %% Create template
        hTM = vision.TemplateMatcher('ROIInputPort', true, 'BestMatchNeighborhoodOutputPort', true);
        
        %% Create motion tracking output display
        % Create a System object to display the original video and the stabilized
        % video.
        hVideoOut = vision.VideoPlayer('Name', 'Video Stabilization');
        hVideoOut.Position(1) = round(0.4*hVideoOut.Position(1));
        hVideoOut.Position(2) = round(.5*(hVideoOut.Position(2)));
        hVideoOut.Position(3:4) = [1050 550];
        
        %% Select tracking area
        % Get target window
        imshow(imageStack(:,:,1));
        title('Select upper left, then lower right target corners and press enter');
        [inputCoordTargetX,inputCoordTargetY] = getpts(gcf);
        close(gcf);
        inputCoordTargetX = round(inputCoordTargetX);
        inputCoordTargetY = round(inputCoordTargetY);
        pos.template_orig = [inputCoordTargetX(1) inputCoordTargetY(1)]; % [x y] upper left corner
        pos.template_size = [inputCoordTargetX(2) inputCoordTargetY(2)] - [inputCoordTargetX(1) inputCoordTargetY(1)];   % [width height]
        
        % Get search window
        imshow(imageStack(:,:,1));
        title('Select upper left search area corner (target box pictured) and press enter');
        rectangle('Position',[pos.template_orig(1) pos.template_orig(2) pos.template_size(1) pos.template_size(2)],'EdgeColor','w');
        [inputCoordSearchX,inputCoordSearchY] = getpts(gcf);
        close(gcf);
        inputCoordSearchX = round(inputCoordSearchX);
        inputCoordSearchY = round(inputCoordSearchY);
        pos.search_border = [abs(inputCoordSearchX(1) - inputCoordTargetX(1)),abs(inputCoordSearchY(1) - inputCoordTargetY(1))];   % max horizontal and vertical displacement
        
        %% Initialize variables for processing loop
        pos.template_center = floor((pos.template_size-1)/2);
        pos.template_center_pos = (pos.template_orig + pos.template_center - 1);
        W = 400; % Width of video in pixels
        H = 400; % Height of video in pixels
        sz = [400 400];
        BorderCols = [1:pos.search_border(1)+4 W-pos.search_border(1)+4:W];
        BorderRows = [1:pos.search_border(2)+4 H-pos.search_border(2)+4:H];
        TargetRowIndices = ...
            pos.template_orig(2)-1:pos.template_orig(2)+pos.template_size(2)-2;
        TargetColIndices = ...
            pos.template_orig(1)-1:pos.template_orig(1)+pos.template_size(1)-2;
        SearchRowIndices = ...
            pos.template_orig(2)-1-pos.search_border(2):pos.template_orig(2)+pos.template_size(2)-2+pos.search_border(2);
        SearchColIndices = ...
            pos.template_orig(1)-1-pos.search_border(1):pos.template_orig(1)+pos.template_size(1)-2+pos.search_border(1);
        SearchRegion = pos.template_orig - pos.search_border - 1;
        Offset = [0 0];
        Target = zeros(length(TargetRowIndices),length(TargetColIndices));
        targetStack = [];
        firstTime = true;
        n = 1;
        targetNum = 0;
        targetAvgNum = 18;
        targetSum = zeros(length(TargetRowIndices),length(TargetColIndices));
        moveDist = [];
        velocity = [];
%         tifLength = tifFrameBounds(2)-tifFrameBounds(1)+1;
%         imageStack = cell(1,tifLength);
        midlineX = round(W/2);
        midlineY = round(H/2);
%         secondsPerFrame = 1/framesPerSecond;
        if ~exist('commentString','var') || isempty(commentString)
            commentString = 'No comments';
        end
        updateSearchTF = false;
        targetFWHM = [];
        targetFWHM2 = [];
        
        %% Stream Processing Loop
        % This is the main processing loop which uses the objects we instantiated
        % above to stabilize the input video.
        for j = 1:size(imageStack,3)
            input = imageStack(:,:,j);
            
            % Find location of Target in the input video frame
            if firstTime
                Idx = int32(pos.template_center_pos);
                motionVector = [0 0];
            else
                IdxPrev = Idx;
                % IdxPrev = int32(pos.template_center_pos);
                
                ROI = [SearchRegion, pos.template_size+2*pos.search_border];
                Idx = hTM(input,Target,ROI);
                
                motionVector = double(Idx-IdxPrev);
            end
            
            if updateSearchTF
                [Offset, SearchRegion] = updatesearch(sz, motionVector, ...
                    SearchRegion, Offset, pos);
            else
                [Offset] = updatesearch(sz, motionVector, ...
                    SearchRegion, Offset, pos);
            end
            
            % Translate video frame to offset the camera motion
            Stabilized = imtranslate(input, Offset, 'linear');
            
            targetMean = mean(Target);
            targetMean(targetMean > .2) = .2;
            targetMean2 = mean(Target,2);
            targetMean2(targetMean2 > .2) = .2;
            plot(targetMean,'b')
            hold on
            plot(targetMean2,'r')
            hold off
            try
                targetFWHM(end+1) = fwhm(1:size(Target,2),targetMean);
                targetFWHM2(end+1) = fwhm(1:size(Target,2),targetMean);
                % Target = Stabilized(round(TargetRowIndices), round(TargetColIndices));
                targetStack(:,:,end+1) = Stabilized(round(TargetRowIndices), round(TargetColIndices));
                if size(targetStack,3) > targetAvgNum
                    targetStack (:,:,1) = [];
                end
                targetSum = sum(targetStack,3);
                Target = targetSum ./ size(targetStack,3);
            catch
                targetFWHM(end+1) = NaN;
            end
            
%             [centersDark, radiiDark] = imfindcircles(input(round(SearchRowIndices), round(SearchColIndices)),[40 70],'ObjectPolarity','dark');
            
            
%             try
%                 currTargetImageProcessed = imfill(bwareaopen(imbinarize(input(round(SearchRowIndices), round(SearchColIndices)),.2),30),'holes');
%                 [~,currTargetImageProcessed] = bwboundaries(currTargetImageProcessed,'noholes');
%                 areaStruct = regionprops(currTargetImageProcessed,'Area');
%                 targetArea(end+1,1) = areaStruct(1).Area;
%             catch
%                 targetArea(end+1,1) = -1;
%             end

            % Add black border for display
            Stabilized(:, BorderCols) = 0;
            Stabilized(BorderRows, :) = 0;
            
            TargetRect = [pos.template_orig-Offset, pos.template_size];
            SearchRegionRect = [SearchRegion, pos.template_size + 2*pos.search_border];
            
            % Draw rectangles on input to show target and search region
%             viscircles(centersDark, radiiDark,'Color','b');
            input = insertShape(input, 'Rectangle', [TargetRect; SearchRegionRect],...
                'Color', 'white');
            % Display the offset (displacement) values on the input image
            txt = sprintf('(%+05.1f pixels,%+05.1f pixels)', Offset);
            input = insertText(input(:,:,1),[1 1],txt,'FontSize',16, ...
                'TextColor', 'white', 'BoxOpacity', 0);
            % Display video
            hVideoOut([input(:,:,1) Stabilized]);
            
%             % Save output video to variable
%             if saveOutputVideoTF
%                 imageStack{1,n} = [input(:,:,1) Stabilized];
%                 n = n + 1;
%             end
            
            % Add pixel motion to data
            if firstTime
                targetPositionPixel = [double(Idx(1)), H-double(Idx(2))+1];
                targetPosition = [0 0];
                firstTime = false;
            else
                motionVector(2) = -motionVector(2);
                moveDist = [moveDist;motionVector];
                targetPositionPixel = [targetPositionPixel;targetPositionPixel(end,:)+motionVector];
            end
        end
        %% Save data
        h(1) = figure('Color','White');
        plot(targetFWHM);
        % Generate plots
        targetPosition(:,1) = medfilt1(targetPosition(:,1));
        targetPosition(:,2) = medfilt1(targetPosition(:,2));
        subtitle = 'Filler subtitle, can fill with info later';
        h(2) = figure('Color','White');
        subplot(2,1,1)
        plot(1:size(moveDist,1),moveDist(:,1),'r')
        title(['\fontsize{20pt}\bf{Object Movement Between Frames}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
        xlabel('Frame')
        ylabel('X Movement (Pixels)')
        grid on
        axis([1 size(moveDist,1) floor(min([moveDist(:,1);moveDist(:,2)])/10)*10 ceil(max([moveDist(:,1);moveDist(:,2)])/10)*10])
        subplot(2,1,2)
        plot(1:size(moveDist,1),moveDist(:,2),'b')
        xlabel('Frame')
        ylabel('Y Movement (Pixels)')
        grid on
        axis([1 size(moveDist,1) floor(min([moveDist(:,1);moveDist(:,2)])/10)*10 ceil(max([moveDist(:,1);moveDist(:,2)])/10)*10])
        %         ProcData.data.pupilDiameter = pupilDiameter;
        %         save(procDataFileID, 'ProcData')   % saves the updated structure to the current folder where all the data is
        %     end
end