
addpath('C:\Users\wdh130\Downloads\Original_resolution_2\Original_resolution_2')
clear t1 t2
tic
f = waitbar(0,'rendering...')
for ii = 1:50
    t1(ii) = toc;
    A_hold = importdata(['Z' num2str(ii+800) '_ch02.tif']);
    A(1:floor(size(A_hold,1)/2),1:floor(size(A_hold,2)/2),ii) = A_hold(1:floor(size(A_hold,1)/2),1:floor(size(A_hold,2)/2));
    t2(ii) = toc;
    t_end = (51-ii)*mean(t2-t1)/60;
    waitbar(ii/51,f,['rendering... eta ' num2str(round(t_end)) ' min'])
end
close(f)

%%
figure(1)
imagesc(sum(A,3))
%vol3d('CData',A(4000:4100,4000:4100,:))

axis equal
hBox = imrect;
roiPosition = wait(hBox);  % Wait for user to double-click
rectangle('Position',roiPosition)

xCoords = [round(roiPosition(1)), round(roiPosition(1)+roiPosition(3))];
yCoords = [round(roiPosition(2)), round(roiPosition(2)+roiPosition(4))];
croppingRectangle = roiPosition;
close

% imshow(Im_crop)
clear B
for ii = 1:46
B(:,:,ii) = A(yCoords(1):yCoords(2),xCoords(1):xCoords(2),ii);
end
% clear A

figure
imagesc(sum(B,3))
axis equal
%caxis([0 max(max(max(B)))])
[BW,xi,yi] = roipoly;
close

clear B2 Z
for ii = 1:size(B,3)
B2(:,:,ii) = B(:,:,ii).*uint16(BW);
Y = conv2(B2(:,:,ii),ones(5),'valid');
Z(:,:,ii) = Y(1:5:end,1:5:end)/25;
end
clear Y A

%%

figure
imagesc(squeeze(sum(B2,2)))
axis equal
%caxis([0 max(max(max(B)))])
[BW,xi,yi] = roipoly;
close

clear B3 Z
for ii = 1:size(B,3)
B3(:,ii,:) = squeeze(B2(:,ii,:)).*uint16(BW);
% Y = conv2(B2(:,:,ii),ones(3),'valid');
% Z(:,:,ii) = Y(1:3:end,1:3:end)/4;
end
%%
%% orient image
B = (Z);
%B = B2;
%B = B3;

figure
dr = 1;
dr_2 = 1;

[x,y,z] = meshgrid(1:size(B,2),1:size(B,1),1:dr_2:30);
v = double(B(1:dr:end,1:dr:end,1:dr_2:end));

%h = slice(x,y,z,v,[],[round(xCoords(3)) round(xCoords(5)+(xCoords(3)-xCoords(5))/2) round(xCoords(5))],[1 25 51]);

h = slice(x,y,z,v,[1:size(B,1)],[1:size(B,2)],[1:dr_2:30]);
set(h,'EdgeColor','none',...
'FaceColor','interp',...
'FaceAlpha','interp')
alpha('color')

alphamap('rampup')
alphamap('decrease',.075)
% alphamap('rampdown')
% alphamap('increase',.025)
% colormap hsv
colormap jet

a = 0.0;
b = -90;

daspect([1 1 0.2])
pbaspect([1 1 0.2])
%set(gca,'BoxStyle','full','Box','on')


grid off
axis off

%[a b] = view;
i = 1;
view(a,b)
%%
B = Z;
B = B2;
B_interp = interp3(single(B),1.2);

figure
h = vol3d('CData',B_interp)
colormap hsv

alphamap('rampup')
alphamap('decrease',.005)

grid off
axis off

daspect([1 1 0.2])
pbaspect([1 1 0.2])

view(151.4043,88.916)



figure
imagesc(sum(B_interp,3))
axis equal
colormap bone

daspect([1 1 0.2])
pbaspect([1 1 0.2])
a = 0.0;
b = -90;
view(a,b)

view(130.4043,81.9375)

camorbit(-42,0,'data',[0 0 1])
camorbit(5,0,'data',[0 0 1])
camorbit(0,90,'camera')
camorbit(0,5,'camera')
camorbit(0,5,'camera')
camorbit(0,5,'camera')
camorbit(0,5,'camera')
camorbit(0,3,'camera')
camorbit(0,-90,'camera')

alphamap('rampup')
alphamap('decrease',.005)

grid off
axis off

%%
for ii = 1:60:360
     camorbit(10,0,'camera')
     pause(0.5)
end

%%
camorbit(-42,0,'data',[0 0 1])
camorbit(-5,0,'data',[0 0 1])
camorbit(0,102,'camera')
camorbit(0,-20,'camera')
camorbit(0,10,'camera')
camorbit(0,1,'camera')
camorbit(0,-1,'camera')
camorbit(0,-1,'camera')
camorbit(0,90,'camera')
camzoom(1.5)
% %%
% camorbit(-1,0,'data',[0 0 1])
% camdolly(-0.30,-0.25,0)
% camzoom(0.9)
% 
% %%
% camorbit(-24,0,'data',[0 0 1])
% camorbit(0,90,'camera')
% camorbit(0,13,'camera')
% camorbit(0,-90,'camera')

%% video

clear t1 t2
f = waitbar(0,'rendering...')
i = 1;
tic
for ii = 1:3:360
    t1(ii) = toc;
    camorbit(3,0,'camera')
    F2(i) = getframe(gcf);
    i = i+1;
    pause(1)
    t2(ii) = toc;
    t_end = (360-ii)*mean(t2-t1)/60;
    waitbar(ii/360,f,['rendering... eta ' num2str(round(t_end)) ' min'])
end
close(f)

video = VideoWriter('vascular_orbit_S1BF_un_18.avi');
video.FrameRate = 30;
open(video);
writeVideo(video, F2);
close(video);

%%

figure
imagesc(squeeze(sum(B2,3)))
axis equal
colormap bone