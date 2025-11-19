close all
fileName = 'C:\Workspace\Code\DrewLab\movementDataLog.mat';
load(fileName)

moveMag = [];
weigth = [];
age = [];
sex = [];
for n = 1:size(moveDataMat,1)
    moveMag(n) = sqrt(moveDataMat{n,7}(1)^2 + moveDataMat{n,7}(2)^2);
    weight(n) = str2double(moveDataMat{n,18});
    age(n) = (datenum(moveDataMat{n,1}(1:6),'yymmdd')-datenum(moveDataMat{n,17},'mm/dd/yy'))/7;
    if strcmp(moveDataMat{n,16},'M')
        sex(n) = 1;
    else
        sex(n) = 0;
    end
end

figure(1)
subplot(2,3,1)
scatter(weight,moveMag)
title('Mouse Weight vs Brain Shift')
xlabel('Mouse Weight (g)')
ylabel('Brain Shift (\mum)')

subplot(2,3,2)
scatter(age,moveMag)
title('Mouse Age vs Brain Shift')
xlabel('Mouse Age (Weeks)')
ylabel('Brain Shift (\mum)')
xlim([0 70])

subplot(2,3,3)
plotSpread({moveMag(~logical(sex)),moveMag(logical(sex))},'categoryIdx',[zeros(length(moveMag(~logical(sex))),1);ones(length(moveMag(logical(sex))),1)],'categoryColors',{'r','b'},'categoryLabels',{'F','M'});
title('Mouse Sex vs Brain Shift')
ylabel('Brain Shift (\mum)')

subplot(2,3,4)
removeDup = unique([age',weight'],'rows');
scatter(removeDup(:,1),removeDup(:,2))
title('Mouse Age vs Weight')
xlabel('Mouse Age (Weeks)')
ylabel('Mouse Weight (g)')
xlim([0 70])

subplot(2,3,5)
removeDup = unique([sex',weight'],'rows');
sexRem = removeDup(:,1);
weightRem = removeDup(:,2);
plotSpread({weightRem(~logical(sexRem)),weightRem(logical(sexRem))},'categoryIdx',[zeros(length(weightRem(~logical(sexRem))),1);ones(length(weightRem(logical(sexRem))),1)],'categoryColors',{'r','b'},'categoryLabels',{'F','M'});
title('Mouse Sex vs Weight')
ylabel('Mouse Weight (g)')