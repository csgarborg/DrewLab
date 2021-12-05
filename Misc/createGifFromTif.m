filename = 'C:\Users\csgar\OneDrive - The Pennsylvania State University\Records\DrewLabMeeting\21-10-05\duralPialVesselShift.gif';
dataFileStr = 'C:\Users\csgar\OneDrive - The Pennsylvania State University\Records\DrewLabMeeting\21-10-05\210824_001_Layer2Filt';
for n = 1:130
    % imshow(im2double(imread([dataFileStr '.tif'],n)));
    % F = getframe(gcf);
    image = uint8(imread([dataFileStr '.tif'],n));
    % [image, ~] = frame2im(F);
%     [imind,cm] = rgb2ind(image,256);
%     if n == 1
%         imwrite(imind,cm,filename,'gif','DelayTime',.1, 'Loopcount',inf);
%     else
%         imwrite(imind,cm,filename,'gif','DelayTime',.1, 'WriteMode','append');
%     end
    if n == 1
        imwrite(image,filename,'gif','DelayTime',.001, 'Loopcount',inf);
    else
        imwrite(image,filename,'gif','DelayTime',.001, 'WriteMode','append');
    end
    close all
end


