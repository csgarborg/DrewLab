close all

layer1 = zeros(512,512);
vec = 2:5:81;
for n = 1:length(vec)
layer1 = layer1 + im2double(medfilt2(imread('F:\19-09-24_MouseExp\190819_032.tif', vec(n))));
end
layer1avg = layer1./length(vec);
% imshow(layer1avg)
layer2 = zeros(512,512);
vec = 3:5:81;
for n = 1:length(vec)
layer2 = layer2 + im2double(medfilt2(imread('F:\19-09-24_MouseExp\190819_032.tif', vec(n))));
end
layer2avg = layer2./length(vec);
% imshow(layer2avg)
layer3 = zeros(512,512);
vec = 4:5:81;
for n = 1:length(vec)
layer3 = layer3 + im2double(medfilt2(imread('F:\19-09-24_MouseExp\190819_032.tif', vec(n))));
end
layer3avg = layer3./length(vec);
% imshow(layer3avg)
layer4 = zeros(512,512);
vec = 5:5:81;
for n = 1:length(vec)
layer4 = layer4 + im2double(medfilt2(imread('F:\19-09-24_MouseExp\190819_032.tif', vec(n))));
end
layer4avg = layer4./length(vec);
% imshow(layer4avg)
layer5 = zeros(512,512);
vec = 6:5:81;
for n = 1:length(vec)
layer5 = layer5 + im2double(medfilt2(imread('F:\19-09-24_MouseExp\190819_032.tif', vec(n))));
end
layer5avg = layer5./length(vec);
% imshow(layer5avg)

volMat = interp3(cat(3,layer1avg,layer2avg,layer3avg,layer4avg,layer5avg),1.5);
figure('units','normalized','outerposition',[0 0 1 1])
h = vol3d('CData',volMat);
colormap hsv

alphamap('rampup')
alphamap('decrease',.005)

grid off
axis off

daspect([1 1 0.2])
pbaspect([1 1 0.5])

view(70,20)


filename = 'F:\DrewLabMeeting\19-10-07\gif6.gif';
dataFileStr = 'F:\19-09-24_MouseExp\190819_032';
camzoom(.8)
for ii = 1:60
    camorbit(-6,0,'data',[0 0 1])
    F = getframe(gcf);
    [image, ~] = frame2im(F);
    [imind,cm] = rgb2ind(image,256);
    if ii == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end