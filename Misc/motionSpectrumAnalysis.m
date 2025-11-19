%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:    motionSpectrumAnalysis
%
% FUNCTION:         motionSpectrumAnalysis(positionData)
%
% DESCRIPTION:      Uses mtspectrumc to look at frequency domain of
%                   movement data
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
% WRITTEN BY:       Spencer Garborg 4/17/20
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PowerX,HzX,ErrorX,PowerY,HzY,ErrorY] = motionSpectrumAnalysis(positionData, plotTF)
dt = 1/30;
tapers = [2 3];
params.Fs = 1/dt;
params.tapers = tapers;
params.fpass = [0.025 7.5];
params.err = [2 0.05]; %use jack-knife resampling confidence intervals p = 0.05
% (could also try multitaper frequency domain bootstrapping (MFDB) for less noisy CI - but jackknife used by so many people in NVC)

[PowerX, HzX, ErrorX] = mtspectrumc(positionData(:,1),params);
[PowerY, HzY, ErrorY] = mtspectrumc(positionData(:,2),params);

if ~exist('plotTF','var')
    plotTF = false;
end

if plotTF
    h(1) = figure('Color','White');
    semilogy(HzX,PowerX,'k')
    hold on
    f = fill([HzX flip(HzX)],[ErrorX(1,:) flip(ErrorX(2,:))],'r','Linestyle','none');
    set(f,'facea',[.2]);
    hold off
    title(['\fontsize{20pt}\bf{X Position Frequency Domain}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Frequency (Hz)')
    ylabel('Power')
    
    h(2) = figure('Color','White');
    semilogy(HzY,PowerY,'k')
    hold on
    f = fill([HzY flip(HzY)],[ErrorY(1,:) flip(ErrorY(2,:))],'r','Linestyle','none');
    set(f,'facea',[.2]);
    hold off
    title(['\fontsize{20pt}\bf{Y Position Frequency Domain}' 10 '\fontsize{10pt}\rm{' subtitle '}' 10 '\fontsize{10pt}\rm{' movementData.commentString '}'])
    xlabel('Frequency (Hz)')
    ylabel('Power')
end
end