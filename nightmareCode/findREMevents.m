function events = findREMevents(emgL, emgR, fs, minDurSec, maxGapSec)
    % Smooth
    winSamp = round(1*fs); % 1s median window
    sL = movmedian(emgL, winSamp);
    sR = movmedian(emgR, winSamp);

    % Adaptive threshold per channel (2-cluster split)
    thL = otsuThreshold(sL);
    thR = otsuThreshold(sR);

    mask = (sL < thL) & (sR < thR);

    % Close small gaps, then remove short runs
    mask = imclose(mask, ones(1, round(maxGapSec*fs)));
    mask = bwareaopen(mask, round(minDurSec*fs));

    % Extract onset/offset sample indices
    d = diff([0 mask 0]);
    onsets  = find(d==1);
    offsets = find(d==-1) - 1;

    events = table(onsets(:), offsets(:), ...
        onsets(:)/fs, offsets(:)/fs, ...
        'VariableNames', {'onsetSample','offsetSample','onsetSec','offsetSec'});
end

function th = otsuThreshold(x)
    xn = rescale(x);
    th_norm = graythresh(xn);
    th = th_norm * range(x) + min(x);
end