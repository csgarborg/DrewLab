
% 1000 (lines per inch) Mesh Copper Grid (2145C from SPI supplies): bar size = 6 microns, hole size = 19 microns
% 2000 (lines per inch) Mesh Copper Grid (2155C from SPI supplies): bar size = 5 microns, hole size = 7.5 microns
function motionTrackingCalibrationGrid2P(filename,barSize,holeSize,objMag,digMag,turnabout,commentString)
close all
% read in and filter image
a = imread(filename);
filtImg = medfilt2(a);
b = filtImg;
c = b;
widthData = [];
heightData = [];
% W = 1;
% run loop while new selections are made, finish calibration if none
newSelections = true;
while newSelections
    try
        % remove NaNs to prevent error
        if ~isempty(widthData)
            widthData(mean(isnan(widthData),2)>0,:) = [];
        end
        if ~isempty(heightData)
            heightData(mean(isnan(heightData),2)>0,:) = [];
        end
        % draw new height and width spots on the image
        if isempty(widthData)
            bShow = insertText(b(:,:,1),[1 1],'Horizontal','FontSize',16, 'TextColor', 'white', 'BoxOpacity', 0);
            %             %         imshow(b);
        else
            bShow = insertText(b(:,:,1),[1 1],'Horizontal','FontSize',16, 'TextColor', 'white', 'BoxOpacity', 0);
            bShow = insertShape(bShow, 'Line', [[widthData(:,4),widthData(:,2),widthData(:,5),widthData(:,2)];[widthData(:,4),widthData(:,2)-5,widthData(:,4),widthData(:,2)+5];[widthData(:,5),widthData(:,2)-5,widthData(:,5),widthData(:,2)+5]], 'Color', 'yellow');
            bShow = insertShape(bShow, 'FilledCircle', [widthData(:,1),widthData(:,2),2*ones(size(widthData,1),1)], 'Color', 'yellow');
            %         imshow(b)
        end
        if isempty(heightData)
            cShow = insertText(c(:,:,1),[1 1],'Vertical','FontSize',16, 'TextColor', 'white', 'BoxOpacity', 0);
            %         imshow(c);
        else
            cShow = insertText(c(:,:,1),[1 1],'Vertical','FontSize',16, 'TextColor', 'white', 'BoxOpacity', 0);
            cShow = insertShape(cShow, 'Line', [[heightData(:,1),heightData(:,4),heightData(:,1),heightData(:,5)];[heightData(:,1)-5,heightData(:,4),heightData(:,1)+5,heightData(:,4)];[heightData(:,1)-5,heightData(:,5),heightData(:,1)+5,heightData(:,5)]], 'Color', 'yellow');
            cShow = insertShape(cShow, 'FilledCircle', [heightData(:,1),heightData(:,2),2*ones(size(heightData,1),1)], 'Color', 'yellow');
            %         imshow(c)
        end
        imshow([bShow cShow]);
        title('Select upper left, then lower right target corners on left picture and press enter. Add a third point to remove previous selections.');
        [inputCoordTargetX,inputCoordTargetY] = getpts(gcf);
        close(gcf);
        % if no new points, finish calibration
        if isempty(inputCoordTargetX)
            newSelections = false;
            break
        end
        inputCoordTargetX = round(inputCoordTargetX);
        inputCoordTargetY = round(inputCoordTargetY);
        % horizontal
        if inputCoordTargetX(2)-inputCoordTargetX(1) > inputCoordTargetY(2)-inputCoordTargetY(1)
            if length(inputCoordTargetX) > 2 && ~isempty(widthData) 
                widthData = widthData(1:end-newWidth,:);
            end
            pixelBlock = filtImg(inputCoordTargetY(1):inputCoordTargetY(2),inputCoordTargetX(1):inputCoordTargetX(2));
            pixelRangeX = inputCoordTargetX(1):inputCoordTargetX(2);
            pixelRangeInt = mean(pixelBlock);
            pixelMidpointY = inputCoordTargetY(1)+((inputCoordTargetY(2)-inputCoordTargetY(1))/2);
            peakSelection = [];
            cont = true;
            [pks,locs] = findpeaks(pixelRangeInt,'MinPeakDistance',20);
            plot(pixelRangeX(locs),pks,'o')
            hold on
            plot(pixelRangeX,pixelRangeInt)
            title('If automation looks correct, enter 1, if manual is needed, enter any other number (automatic peak selection can be fine tuned by chnanging MinPeakDistance to match data')
            autoTF = input('Enter 1 if auto is correct: ');
            close
            if autoTF == 1
                peakSelection(1) = pixelRangeX(locs(1));
                for i = 2:length(locs)-1
                    peakSelection(end+1:end+2) = pixelRangeX(locs(i));
                end
                peakSelection(end+1) = pixelRangeX(locs(end));
            else
                h(1) = figure('Color','White','Name','Select frame starts and stops, select behind line to erase, right and left arrow key to move plot, press enter when finished','NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
                while cont
                    plot(pixelRangeX,pixelRangeInt,'b');
                    hold on
                    %     if length(frameSelection) == 1
                    %         line([frameSelection(1) frameSelection(1)],[min(targetPosition(:,1)) max(targetPosition(:,1))],'Color','g','LineStyle','--');
                    %     end
                    for n = 1:length(peakSelection)
                        line([peakSelection(n) peakSelection(n)],[0 max(pixelRangeInt)],'Color','k','LineStyle','--');
                    end
                    hold off
                    %     title(['\fontsize{20pt}\bf{Object Position per Frame}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
                    xlabel('X Position (Pixels)')
                    ylabel('Mean Pixel Intensity')
                    axis([0 512 0 ceil(max(pixelRangeInt)/100)*100])
                    grid on
                    
                    [selectedCoordX,~,arrowKey] = ginput(1);
                    
                    if isempty(arrowKey)
                        cont = false;
                        close(h(1))
                    elseif arrowKey == 1
                        if ~isempty(peakSelection) && round(selectedCoordX) <= peakSelection(end)
                            peakSelection = peakSelection(1:end-1);
                        else
                            peakSelection(end+1) = round(selectedCoordX);
                        end
                    end
                end
            end
            newWidth = 0;
            for n = 1:2:length(peakSelection)-1
                pixelChunkInt = pixelRangeInt(peakSelection(n) <= pixelRangeX & pixelRangeX <= peakSelection(n+1));
                [width,leftBound,rightBound,polarity] = fwhm(peakSelection(n):peakSelection(n+1),pixelChunkInt);
                widthData(end+1,:) = [leftBound+((rightBound-leftBound)/2), pixelMidpointY, width, leftBound, rightBound, polarity];
                newWidth = newWidth + 1;
            end
            % vertical
        else
            if length(inputCoordTargetX) > 2 && ~isempty(heightData) 
                heightData = heightData(1:end-newHeight,:);
            end
            pixelBlock = filtImg(inputCoordTargetY(1):inputCoordTargetY(2),inputCoordTargetX(1):inputCoordTargetX(2));
            pixelRangeY = inputCoordTargetY(1):inputCoordTargetY(2);
            pixelRangeInt = mean(pixelBlock,2);
            pixelMidpointX = inputCoordTargetX(1)+((inputCoordTargetX(2)-inputCoordTargetX(1))/2);
            peakSelection = [];
            cont = true;
            [pks,locs] = findpeaks(pixelRangeInt,'MinPeakDistance',20);
            plot(pixelRangeY(locs),pks,'o')
            hold on
            plot(pixelRangeY,pixelRangeInt)
            title('If automation looks correct, enter 1, if manual is needed, enter any other number (automatic peak selection can be fine tuned by chnanging MinPeakDistance to match data')
            autoTF = input('Enter 1 if auto is correct: ');
            close
            if autoTF == 1
                peakSelection(1) = pixelRangeY(locs(1));
                for i = 2:length(locs)-1
                    peakSelection(end+1:end+2) = pixelRangeY(locs(i));
                end
                peakSelection(end+1) = pixelRangeY(locs(end));
            else
                h(1) = figure('Color','White','Name','Select frame starts and stops, select behind line to erase, right and left arrow key to move plot, press enter when finished','NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
                while cont
                    plot(pixelRangeY,pixelRangeInt,'b');
                    hold on
                    %     if length(frameSelection) == 1
                    %         line([frameSelection(1) frameSelection(1)],[min(targetPosition(:,1)) max(targetPosition(:,1))],'Color','g','LineStyle','--');
                    %     end
                    for n = 1:length(peakSelection)
                        line([peakSelection(n) peakSelection(n)],[0 max(pixelRangeInt)],'Color','k','LineStyle','--');
                    end
                    hold off
                    %     title(['\fontsize{20pt}\bf{Object Position per Frame}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
                    xlabel('Y Position (Pixels)')
                    ylabel('Mean Pixel Intensity')
                    axis([0 512 0 ceil(max(pixelRangeInt)/100)*100])
                    grid on
                    
                    [selectedCoordX,~,arrowKey] = ginput(1);
                    
                    if isempty(arrowKey)
                        cont = false;
                        close(h(1))
                    elseif arrowKey == 1
                        if ~isempty(peakSelection) && round(selectedCoordX) <= peakSelection(end)
                            peakSelection = peakSelection(1:end-1);
                        else
                            peakSelection(end+1) = round(selectedCoordX);
                        end
                    end
                end
            end
            newHeight = 0;
            for n = 1:2:length(peakSelection)-1
                pixelChunkInt = pixelRangeInt(peakSelection(n) <= pixelRangeY & pixelRangeY <= peakSelection(n+1));
                [height,topBound,bottomBound,polarity] = fwhm(peakSelection(n):peakSelection(n+1),pixelChunkInt);
                heightData(end+1,:) = [pixelMidpointX, topBound+((bottomBound-topBound)/2), height, topBound, bottomBound, polarity];
                newHeight = newHeight + 1;
            end
        end
    catch
        disp('Error occured, please retry')
    end
end
close all

subtitle = [num2str(barSize) ' ' char(181) 'm Bar Size, ' num2str(holeSize) ' ' char(181) 'm Hole Size, ' num2str(objMag*digMag) 'x Magnification (' num2str(objMag) 'x Objective, ' num2str(digMag) 'x Digital), Turnabout = ' num2str(turnabout)];
if ~exist('commentString','var') || isempty(commentString)
    commentString = 'No comments';
end
W = size(a,2);
H = size(a,1);
h(1) = figure('Color','White');
micronLengthValWidth = ones(size(widthData,1),1);
micronLengthValWidth(widthData(:,6) == 1,1) = barSize;
micronLengthValWidth(widthData(:,6) == -1,1) = holeSize;
surfaceCalibFitX = fit(widthData(:,1:2),micronLengthValWidth./widthData(:,3),'poly55');
plot(surfaceCalibFitX,widthData(:,1:2),micronLengthValWidth./widthData(:,3));
axis([1 W 1 H 0 ceil(max(micronLengthValWidth./widthData(:,3)))])
title(['\fontsize{20pt}\bf{Surface Calibration Fit (X)}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixels (X)')
ylabel('Pixels (Y)')
zlabel('\mum/Pixel')

h(2) = figure('Color','White');
micronLengthValHeight = ones(size(heightData,1),1);
micronLengthValHeight(heightData(:,6) == 1,1) = barSize;
micronLengthValHeight(heightData(:,6) == -1,1) = holeSize;
surfaceCalibFitY = fit(heightData(:,1:2),micronLengthValHeight./heightData(:,3),'poly55');
plot(surfaceCalibFitY,heightData(:,1:2),micronLengthValHeight./heightData(:,3));
axis([1 W 1 H 0 ceil(max(micronLengthValHeight./heightData(:,3)))])
title(['\fontsize{20pt}\bf{Surface Calibration Fit (Y)}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixels (X)')
ylabel('Pixels (Y)')
zlabel('\mum/Pixel')

h(3) = figure('Color','White');
scatter(widthData(:,1),widthData(:,2),50,micronLengthValWidth./widthData(:,3),'filled');
colorbar;
k = colorbar;
ylabel(k,'\mum/pixel','FontSize',15)
hold off
axis equal
axis([1 W 1 H])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title(['\fontsize{20pt}\bf{Position of X Calibration Measurements}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixels')
ylabel('Pixels')

h(4) = figure('Color','White');
scatter(heightData(:,1),heightData(:,2),50,micronLengthValHeight./heightData(:,3),'filled');
colorbar;
k = colorbar;
ylabel(k,'\mum/pixel','FontSize',15)
hold off
axis equal
axis([1 W 1 H])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title(['\fontsize{20pt}\bf{Position of Y Calibration Measurements}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' commentString '}'])
xlabel('Pixels')
ylabel('Pixels')

h(5) = figure;
imshow([bShow cShow]);

% Save figures to single .fig file
savefig(h,[filename(1:end-4) '_calibrationPlots.fig']);

% Save calibration values
if exist('C:\Workspace\Code\DrewLab\calibrationValues.mat')
    load('C:\Workspace\Code\DrewLab\calibrationValues.mat');
end

fileName = filename(end-13:end-4);
if size(widthData,2) >= 4
    calibrationValues.(['file_' fileName]).pixelMidpointXWidth = widthData(:,1);
    calibrationValues.(['file_' fileName]).pixelMidpointYWidth = widthData(:,2);
    calibrationValues.(['file_' fileName]).width = widthData(:,3);
    calibrationValues.(['file_' fileName]).leftBoundWidth = widthData(:,4);
    calibrationValues.(['file_' fileName]).rightBoundWidth = widthData(:,5);
    calibrationValues.(['file_' fileName]).polarityWidth = widthData(:,6);
    calibrationValues.(['file_' fileName]).pixelMidpointXHeight = heightData(:,1);
    calibrationValues.(['file_' fileName]).pixelMidpointYHeight = heightData(:,2);
    calibrationValues.(['file_' fileName]).height = heightData(:,3);
    calibrationValues.(['file_' fileName]).topBoundHeight = heightData(:,4);
    calibrationValues.(['file_' fileName]).bottomBoundHeight = heightData(:,5);
    calibrationValues.(['file_' fileName]).polarityHeight = heightData(:,6);
    calibrationValues.(['file_' fileName]).positivePolarityLengthInMicrons = barSize;
    calibrationValues.(['file_' fileName]).negativePolarityLengthInMicrons = holeSize;
    calibrationValues.(['file_' fileName]).surfaceCalibFitX = surfaceCalibFitX;
    calibrationValues.(['file_' fileName]).surfaceCalibFitY = surfaceCalibFitY;
    calibrationValues.(['file_' fileName]).imageWidthPixels = W;
    calibrationValues.(['file_' fileName]).imageHeightPixels = H;
    
    save('C:\Workspace\Code\DrewLab\calibrationValues.mat','calibrationValues');
    
    if ~exist('C:\Workspace\Code\DrewLab\calibrationLog.csv')
        dlmwrite('C:\Workspace\Code\DrewLab\calibrationLog.csv',[]);
        fid = fopen('C:\Workspace\Code\DrewLab\calibrationLog.csv','wt');
        fprintf(fid,'%s, %s, %s, %s, %s, %s\n','Date','Filename','Objective Zoom','Digital Zoom','Frame Height (Pixels)','Turnabout');
        fidClose = fclose(fid);
    end
    fid = fopen('C:\Workspace\Code\DrewLab\calibrationLog.csv','a');
    fprintf(fid,'%s, %s, %.0f, %.0f, %.0f, %.0f\n',datestr(now,'yyyy/mm/dd HH:MM'),fileName,objMag,digMag,H,turnabout);
    fidClose = fclose(fid);
end
