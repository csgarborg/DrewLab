cd('C:\Users\wdh130\Documents\MATLAB\Spencer')
clear
clc

fname = '181219_014.tif';
info = imfinfo(fname);
num_images = numel(info);
num_images = 300;
for k = 1:num_images
    A{1,k} = imread(fname, k);
    % ... Do something with image A ...
    %imwrite(A{k},'test.tif','WriteMode','append')
end

B = cell2mat(A);
B = double(reshape(B,512,512,num_images));

fs = 30;
[y1 d1] = lowpass(squeeze(B(250,250,:)),100,fs);

figure
pspectrum(y1,fs)


y1 = squeeze(B(350,350,:));

lpFilt = designfilt('lowpassiir','FilterOrder',8, ...
         'PassbandFrequency',1,'PassbandRipple',0.2, ...
         'SampleRate',30);
fvtool(lpFilt)

output = filter(lpFilt,y1);
figure(1), subplot(2,2,[1 2]), hold on
plot(y1)
plot(output)
subplot(2,2,3)
pspectrum(y1,fs)
subplot(2,2,4)
pspectrum(output,fs)

output = zeros(512,512,300);
for xx = 1:512
    for yy = 1:512
        y1 = squeeze(B(xx,yy,:));
        output(xx,yy,:) = filter(lpFilt,y1);
    end
end

v = VideoWriter('Fourier_filtered_and_original.avi');
open(v);
%Generate initial data and set axes and figure properties.
fig = figure;
for k = 1:300
    subplot(1,2,2)
    imagesc(output(:,:,k))
    subplot(1,2,1)
    imagesc(B(:,:,k))
    frame = getframe(fig);
    writeVideo(v,frame);
end
% close(v);
% 
% v = VideoWriter('raw.avi');
% open(v);
% %Generate initial data and set axes and figure properties.
% figure
% for k = 1:num_images
%    imagesc(B(:,:,k))
%    frame = getframe;
%    writeVideo(v,frame);
% end
% close(v);
% 
% 
% figure
% subplot(1,2,1)
% imagesc(output(:,:,50))
% subplot(1,2,2)
% imagesc(output(:,:,500))
% 
% for kk = 1:5;
% blur_value = 10;
% 
% original = A{kk};
% blurred = imgaussfilt(original,blur_value);
% 
% sub_im = original - blurred;
% 
% test = original - sub_im;
% 
% figure(kk)
% subplot(2,2,1)
% imshow(original), title('original')
% subplot(2,2,2)
% imshow(blurred), title('blurred')
% subplot(2,2,3)
% imshow(sub_im), title('subtracted')
% subplot(2,2,4)
% imshow(test), title('test')
% end
