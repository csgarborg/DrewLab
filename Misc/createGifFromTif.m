filename = 'E:\20-01-15_MouseExp\200115_003_2LMovement.gif';
dataFileStr = 'E:\20-01-15_MouseExp\200115_003_2LMovement';
for n = 1:1945
    % imshow(im2double(imread([dataFileStr '.tif'],n)));
    % F = getframe(gcf);
    image = im2double(imread([dataFileStr '.tif'],n));
    % [image, ~] = frame2im(F);
    [imind,cm] = rgb2ind(image,256);
    if n == 1
        imwrite(imind,cm,filename,'gif','DelayTime',.02, 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','DelayTime',.02, 'WriteMode','append');
    end
    close all
end


