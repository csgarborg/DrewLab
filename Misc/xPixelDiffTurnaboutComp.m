close all

load('calibrationValues.mat')

% Turnabout: 200
i1 = 150 <= calibrationValues.file_190311_001.pixelDiffPosX(:,1) & calibrationValues.file_190311_001.pixelDiffPosX(:,1) <= 350;
p1 = polyfit(calibrationValues.file_190311_001.pixelDiffPosX(i1,1),calibrationValues.file_190311_001.pixelDiffX(i1),2);
x1 = 150:.5:350;
y1 = polyval(p1,x1);
figure(1)
plot(calibrationValues.file_190311_001.pixelDiffPosX(:,1),calibrationValues.file_190311_001.pixelDiffX,'*',x1,y1,'k')
title('200')
axis([1 512 0 35])

% Turnabout: 204
i9 = 150 <= calibrationValues.file_190319_002.pixelDiffPosX(:,1) & calibrationValues.file_190319_002.pixelDiffPosX(:,1) <= 350;
p9 = polyfit(calibrationValues.file_190319_002.pixelDiffPosX(i9,1),calibrationValues.file_190319_002.pixelDiffX(i9),2);
x9 = 150:.5:350;
y9 = polyval(p9,x9);
figure(2)
plot(calibrationValues.file_190319_002.pixelDiffPosX(:,1),calibrationValues.file_190319_002.pixelDiffX,'*',x9,y9,'k')
title('204')
axis([1 512 0 35])
% Turnabout: 208
i10 = 150 <= calibrationValues.file_190319_001.pixelDiffPosX(:,1) & calibrationValues.file_190319_001.pixelDiffPosX(:,1) <= 350;
p10 = polyfit(calibrationValues.file_190319_001.pixelDiffPosX(i10,1),calibrationValues.file_190319_001.pixelDiffX(i10),2);
x10 = 150:.5:350;
y10 = polyval(p10,x10);
figure(3)
plot(calibrationValues.file_190319_001.pixelDiffPosX(:,1),calibrationValues.file_190319_001.pixelDiffX,'*',x10,y10,'k')
title('208')
axis([1 512 0 35])

% Turnabout: 212
i8 = 150 <= calibrationValues.file_190315_008.pixelDiffPosX(:,1) & calibrationValues.file_190315_008.pixelDiffPosX(:,1) <= 350;
p8 = polyfit(calibrationValues.file_190315_008.pixelDiffPosX(i8,1),calibrationValues.file_190315_008.pixelDiffX(i8),2);
x8 = 150:.5:350;
y8 = polyval(p8,x8);
figure(4)
plot(calibrationValues.file_190315_008.pixelDiffPosX(:,1),calibrationValues.file_190315_008.pixelDiffX,'*',x8,y8,'k')
title('212')
axis([1 512 0 35])

% Turnabout: 225
i7 = 150 <= calibrationValues.file_190315_007.pixelDiffPosX(:,1) & calibrationValues.file_190315_007.pixelDiffPosX(:,1) <= 350;
p7 = polyfit(calibrationValues.file_190315_007.pixelDiffPosX(i7,1),calibrationValues.file_190315_007.pixelDiffX(i7),2);
x7 = 150:.5:350;
y7 = polyval(p7,x7);
figure(5)
plot(calibrationValues.file_190315_007.pixelDiffPosX(:,1),calibrationValues.file_190315_007.pixelDiffX,'*',x7,y7,'k')
title('225')
axis([1 512 0 35])

% Turnabout: 250
i3 = 150 <= calibrationValues.file_190315_002.pixelDiffPosX(:,1) & calibrationValues.file_190315_002.pixelDiffPosX(:,1) <= 350;
p3 = polyfit(calibrationValues.file_190315_002.pixelDiffPosX(i3,1),calibrationValues.file_190315_002.pixelDiffX(i3),2);
x3 = 150:.5:350;
y3 = polyval(p3,x3);
figure(6)
plot(calibrationValues.file_190315_002.pixelDiffPosX(:,1),calibrationValues.file_190315_002.pixelDiffX,'*',x3,y3,'k')
title('250')
axis([1 512 0 35])

% Turnabout: 265
i2 = 150 <= calibrationValues.file_190315_001.pixelDiffPosX(:,1) & calibrationValues.file_190315_001.pixelDiffPosX(:,1) <= 350;
p2 = polyfit(calibrationValues.file_190315_001.pixelDiffPosX(i2,1),calibrationValues.file_190315_001.pixelDiffX(i2),2);
x2 = 150:.5:350;
y2 = polyval(p2,x2);
figure(7)
plot(calibrationValues.file_190315_001.pixelDiffPosX(:,1),calibrationValues.file_190315_001.pixelDiffX,'*',x2,y2,'k')
title('265')
axis([1 512 0 35])

% Turnabout: 278
i4 = 150 <= calibrationValues.file_190315_003.pixelDiffPosX(:,1) & calibrationValues.file_190315_003.pixelDiffPosX(:,1) <= 350;
p4 = polyfit(calibrationValues.file_190315_003.pixelDiffPosX(i4,1),calibrationValues.file_190315_003.pixelDiffX(i4),2);
x4 = 150:.5:350;
y4 = polyval(p4,x4);
figure(8)
plot(calibrationValues.file_190315_003.pixelDiffPosX(:,1),calibrationValues.file_190315_003.pixelDiffX,'*',x4,y4,'k')
title('278')
axis([1 512 0 35])

% Turnabout: 300
i5 = 150 <= calibrationValues.file_190315_005.pixelDiffPosX(:,1) & calibrationValues.file_190315_005.pixelDiffPosX(:,1) <= 350;
p5 = polyfit(calibrationValues.file_190315_005.pixelDiffPosX(i5,1),calibrationValues.file_190315_005.pixelDiffX(i5),2);
x5 = 150:.5:350;
y5 = polyval(p5,x5);
figure(9)
plot(calibrationValues.file_190315_005.pixelDiffPosX(:,1),calibrationValues.file_190315_005.pixelDiffX,'*',x5,y5,'k')
title('300')
axis([1 512 0 35])

% Turnabout: 350
i6 = 150 <= calibrationValues.file_190315_006.pixelDiffPosX(:,1) & calibrationValues.file_190315_006.pixelDiffPosX(:,1) <= 350;
p6 = polyfit(calibrationValues.file_190315_006.pixelDiffPosX(i6,1),calibrationValues.file_190315_006.pixelDiffX(i6),2);
x6 = 150:.5:350;
y6 = polyval(p6,x6);
figure(10)
plot(calibrationValues.file_190315_006.pixelDiffPosX(:,1),calibrationValues.file_190315_006.pixelDiffX,'*',x6,y6,'k')
title('350')
axis([1 512 0 35])

figure(11)
plot(x1,y1,'r',x9,y9,'k-.',x10,y10,'b-.',x8,y8,'k--',x7,y7,'g--',x3,y3,'g',x2,y2,'b',x4,y4,'k',x5,y5,'r--',x6,y6,'b--')
legend('200','204','208','212','225','250','265','278','300','350')