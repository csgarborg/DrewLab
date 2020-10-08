filename = 'Y:\Figures\PSFBeadMovement3D.gif';
dataFileStr = 'Y:\Figures\PSFBeadMovement3D';
for n = 1:77
    % imshow(im2double(imread([dataFileStr '.tif'],n)));
    % F = getframe(gcf);
    image = im2double(imread([dataFileStr '.tif'],n));
    % [image, ~] = frame2im(F);
    [imind,cm] = rgb2ind(image,256);
    if n == 1
        imwrite(imind,cm,filename,'gif','DelayTime',.1, 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','DelayTime',.1, 'WriteMode','append');
    end
    close all
end


