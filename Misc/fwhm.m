function [width,ttrail,tlead,Pol] = fwhm(x,y)

% function width = fwhm(x,y)
%
% Full-Width at Half-Maximum (FWHM) of the waveform y(x)
% and its polarity.
% The FWHM result in 'width' will be in units of 'x'
%
%
% Rev 1.2, April 2006 (Patrick Egan)


yOrig = y;
y = y / max(y);
N = length(y);
lev50 = 0.5;
if (y(1) < lev50) && (median(y) > lev50)                 % find index of center (max or min) of pulse
    [~,centerindex]=max(y);
    Pol = +1;
    y1 = y;
    y2 = y;
%     disp('Pulse Polarity = Positive')
else
    [~,centerindex]=min(y);
    Pol = -1;
    y1 = yOrig/max(yOrig(1:centerindex));
    y2 = yOrig/max(yOrig(centerindex:end));
%     disp('Pulse Polarity = Negative')
end
i = 2;
while sign(y1(i)-lev50) == sign(y1(i-1)-lev50)
    i = i+1;
end                                   %first crossing is between v(i-1) & v(i)
interp = (lev50-y1(i-1)) / (y1(i)-y1(i-1));
tlead = x(i-1) + interp*(x(i)-x(i-1));
i = centerindex+1;                    %start search for next crossing at center
while ((sign(y2(i)-lev50) == sign(y2(i-1)-lev50)) && (i <= N-1))
    i = i+1;
end
if i ~= N || (sign(y2(i)-lev50) ~= sign(y2(i-1)-lev50))
    Ptype = 1;  
%     disp('Pulse is Impulse or Rectangular with 2 edges')
    interp = (lev50-y2(i-1)) / (y2(i)-y2(i-1));
    ttrail = x(i-1) + interp*(x(i)-x(i-1));
    width = ttrail - tlead;
else
    Ptype = 2; 
%     disp('Step-Like Pulse, no second edge')
    ttrail = NaN;
    width = NaN;
end
