%http://mres.uni-potsdam.de/index.php/2017/09/14/principal-component-analysis-in-6-steps/

clear, clc, close all

rng(0)
data(:,1) = randn(30,1);
data(:,2) = 3.4 + 1.2 * data(:,1);
data(:,2) = data(:,2) + 0.2*randn(size(data(:,1)));
data = sortrows(data,1);

figure
axes('LineWidth',0.6,...
    'FontName','Helvetica',...
    'FontSize',8,...
    'XAxisLocation','Origin',...
    'YAxisLocation','Origin');
line(data(:,1),data(:,2),...
    'LineStyle','None',...
    'Marker','o');
axis equal

data(:,1) = data(:,1)-mean(data(:,1));
data(:,2) = data(:,2)-mean(data(:,2));

C = cov(data);

[V,D] = eig(C);

figure
axes('LineWidth',0.6,...
    'FontName','Helvetica',...
    'FontSize',8,...
    'XAxisLocation','Origin',...
    'YAxisLocation','Origin');
line(data(:,1),data(:,2),...
    'LineStyle','None',...
    'Marker','o');
line([0 V(1,1)],[0 V(2,1)],...
    'Color',[0.8 0.5 0.3],...
    'LineWidth',0.75);
% Line I want showing direction of most variance
line([0 V(1,2)],[0 V(2,2)],...
    'Color',[0.8 0.5 0.3],...
    'LineWidth',0.75);
axis equal

norm(V(:,1))
norm(V(:,2))
dot(V(:,1),V(:,2))

newdata = V * data';
newdata = newdata';
newdata = fliplr(newdata);

var(newdata)
var(newdata)/sum(var(newdata))
variance = D / sum(D(:))

figure
axes('LineWidth',0.6,...
    'FontName','Helvetica',...
    'FontSize',8,...
    'XAxisLocation','Origin',...
    'YAxisLocation','Origin')
line(newdata(:,1),newdata(:,2),...
    'LineStyle','None',...
    'Marker','o');
axis equal

[coeff,newdata,latend,tsd,variance] = pca(data);
newdata
variance

corrcoef(newdata)