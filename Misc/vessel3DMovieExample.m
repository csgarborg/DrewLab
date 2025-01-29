close all

dataFileStr = 'F:\19-09-24_MouseExp\190924_002';
vec = 7:6:200-6;

for i = 1:length(vec)
    
    layer1 = zeros(512,512);
    layer1 = im2double(medfilt2(imread('F:\19-09-24_MouseExp\190924_002.tif', vec(i))));
    
    layer2 = zeros(512,512);
    layer2 = im2double(medfilt2(imread('F:\19-09-24_MouseExp\190924_002.tif', vec(i)+1)));
    
    layer3 = zeros(512,512);
    layer3 = im2double(medfilt2(imread('F:\19-09-24_MouseExp\190924_002.tif', vec(i)+2)));
    
    layer4 = zeros(512,512);
    layer4 = im2double(medfilt2(imread('F:\19-09-24_MouseExp\190924_002.tif', vec(i)+3)));
    
    layer5 = zeros(512,512);
    layer5 = im2double(medfilt2(imread('F:\19-09-24_MouseExp\190924_002.tif', vec(i)+4)));
    
    
    volMat = interp3(cat(3,layer1,layer2,layer3,layer4,layer5),1.5);
    close all
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
    
    F = getframe(gcf);
    [image, ~] = frame2im(F);
    if n == 1
        imwrite(image,[dataFileStr '_3dMovie.tif']);
    else
        imwrite(image,[dataFileStr '_3dMovie.tif'],'WriteMode','append');
    end
end
close all