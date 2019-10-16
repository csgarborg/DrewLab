close all

layer1 = zeros(512,512);
vec = 7:6:200;
for n = 1:length(vec)
layer1 = layer1 + im2double(medfilt2(imread('F:\19-09-24_MouseExp\190924_002.tif', vec(n))));
end
layer1avg = layer1./length(vec);
% imshow(layer1avg)
layer2 = zeros(512,512);
vec = 8:6:200;
for n = 1:length(vec)
layer2 = layer2 + im2double(medfilt2(imread('F:\19-09-24_MouseExp\190924_002.tif', vec(n))));
end
layer2avg = layer2./length(vec);
% imshow(layer2avg)
layer3 = zeros(512,512);
vec = 9:6:200;
for n = 1:length(vec)
layer3 = layer3 + im2double(medfilt2(imread('F:\19-09-24_MouseExp\190924_002.tif', vec(n))));
end
layer3avg = layer3./length(vec);
% imshow(layer3avg)
layer4 = zeros(512,512);
vec = 10:6:200;
for n = 1:length(vec)
layer4 = layer4 + im2double(medfilt2(imread('F:\19-09-24_MouseExp\190924_002.tif', vec(n))));
end
layer4avg = layer4./length(vec);
% imshow(layer4avg)
layer5 = zeros(512,512);
vec = 11:6:200;
for n = 1:length(vec)
layer5 = layer5 + im2double(medfilt2(imread('F:\19-09-24_MouseExp\190924_002.tif', vec(n))));
end
layer5avg = layer5./length(vec);
% imshow(layer5avg)

volMat = interp3(cat(3,layer1avg,layer2avg,layer3avg,layer4avg,layer5avg),1.5);
figure('units','normalized','outerposition',[0 0 1 1])
h = vol3d('CData',volMat);
colormap hsv

view(70,20)

alphamap('rampup')
alphamap('decrease',.05)

grid off
axis off

daspect([1 1 0.2])
pbaspect([1 1 0.1])

filename = 'F:\DrewLabMeeting\19-10-07\gif7.gif';
dataFileStr = 'F:\19-09-24_MouseExp\190924_002';
camzoom(2.5)
for ii = 1:90
    camorbit(-2,0,'data',[0 0 1])
    F = getframe(gcf);
    [image, ~] = frame2im(F);
    [imind,cm] = rgb2ind(image,256);
    if ii == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end