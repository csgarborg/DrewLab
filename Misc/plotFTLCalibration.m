%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:    plotFTLCalibration
%
% FUNCTION:         plotFTLCalibration
%
% DESCRIPTION:      Plots all saved z calibrations of FTL
%
% INPUT:
%
% VARIABLES:
%
% OUTPUT:
%
% FUNCTIONS USED:
%
% LIBARIES USED:
%
% NOTES:
%
% WRITTEN BY:       Spencer Garborg 11/08/19
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotFTLCalibration(inputZVec)

close all

load('C:\Workspace\Code\DrewLab\calibrationValuesFTL.mat');

fn = fieldnames(calibrationValues);

allDataX = [];
allDataY = [];
for n = 1:numel(fn)
    allDataX = [allDataX calibrationValues.(fn{n}).diopterVals];
    allDataY = [allDataY calibrationValues.(fn{n}).zMatchMicrons];
end
coefficients = polyfit(allDataX, allDataY, 3);
xFit = linspace(min(allDataX), max(allDataX), 1000);
yFit = polyval(coefficients , xFit);
yZeroDiff = yFit(1);
yFit = yFit - yZeroDiff;
uniqueX = unique(allDataX);
for n = 1:length(uniqueX)
    i = allDataX == uniqueX(n);
    matchingY = allDataY(i);
    stdX(n) = uniqueX(n);
    stdYPlus(n) =  polyval(coefficients , uniqueX(n)) - yZeroDiff + std(matchingY);
    stdYMinus(n) = polyval(coefficients , uniqueX(n)) - yZeroDiff - std(matchingY);
end

figure(1);
subplot(2,1,1)
for n = 1:numel(fn)
%     plot(calibrationValues.(fn{n}).diopterVals,calibrationValues.(fn{n}).zMatchMicrons - yZeroDiff,'-*');
    hold on
end

plot(xFit, yFit, '--', 'LineWidth', 2);
f = fill([stdX, fliplr(stdX)], [stdYPlus, fliplr(stdYMinus)], 'r','Linestyle','none');
set(f,'facea',[.2]);
hold off
xlabel('Diopters (meters^{-1})')
ylabel('\Delta Z (microns)')
title('Diopter Values vs. Focal Plane Position in Z')
grid on
% axis([-2.1 1.5 0 40])


subplot(2,1,2)
plot(xFit, yFit, '--', 'LineWidth', 2);
hold on
if exist('inputZVec','var')
    inputZVec = inputZVec + xFit(1);
    yFitInput = polyval(coefficients , inputZVec);
    for n = 1:length(yFitInput)
        plot(inputZVec(n),yFitInput(n),'k*','MarkerSize',15);
    end
    hold off
    text(.25,ceil(max(yFitInput))-.2,['Diopter Values Output: ' num2str(round(yFitInput,2))])
end
ylabel('Diopters (meters^{-1})')
xlabel('Z (microns)')
title('Diopter Values vs. Focal Plane Position in Z')
grid on
if exist('inputZVec','var')
    axis([0 ceil(max(inputZVec)) -2.1 ceil(max(yFitInput))])
else
    axis([0 40 -2.1 1.5])
end
end
