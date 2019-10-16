filename = 'F:\DrewLabMeeting\19-10-07\gif10.gif';
dataFileStr = 'E:\19-08-14_MouseExp\190814_002';
for n = 1:300
imshow(im2double(imread([dataFileStr '.tif'],n+800)));
F = getframe(gcf);
[image, ~] = frame2im(F);
[imind,cm] = rgb2ind(image,256);
if n == 1
imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
else
imwrite(imind,cm,filename,'gif','WriteMode','append');
end
close all
end


