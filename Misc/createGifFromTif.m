filename = 'H:\21-02-17_MouseExp\movementExample.gif';
dataFileStr = 'H:\21-02-17_MouseExp\210217_006_ProcessedMovementReducedFlattened';
for n = 1:3999
    % imshow(im2double(imread([dataFileStr '.tif'],n)));
    % F = getframe(gcf);
    image = uint8(imread([dataFileStr '.tif'],n));
%     image = im2double(imread([dataFileStr '.tif'],n));
    % [image, ~] = frame2im(F);
%     [imind,cm] = rgb2ind(image,256);
%     if n == 1
%         imwrite(imind,cm,filename,'gif','DelayTime',.1, 'Loopcount',inf);
%     else
%         imwrite(imind,cm,filename,'gif','DelayTime',.1, 'WriteMode','append');
%     end
    if n == 1
        imwrite(image(:,:,2),filename,'gif','DelayTime',.001, 'Loopcount',inf);
    else
        imwrite(image(:,:,2),filename,'gif','DelayTime',.001, 'WriteMode','append');
    end
    close all
end


