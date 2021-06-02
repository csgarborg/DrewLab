function [yMean cIntFillPts] = getCIntMeanAndFillPts(yMat,CIntPercentage)

N = size(yMat,1);                                      % Number of ‘Experiments’ In Data Set
yMean = sgolayfilt(mean(yMat),3,13);                                    % Mean Of All Experiments At Each Value Of ‘x’
ySEM = std(yMat)/sqrt(N);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
cIntVal = (1 - (CIntPercentage/100)) * 0.5;
CI95 = tinv([cIntVal 1-cIntVal], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
cIntFillPts = sgolayfilt([(yMean + yCI95(1,:)) (flip(yMean) + flip(yCI95(2,:)))],3,13);
end