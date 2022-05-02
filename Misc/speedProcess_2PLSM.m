function [Speeddata, SpeedaccBinary, SpeedaveBinary] = speedProcess_2PLSM(Speed, Fs, Fr, numImage)
    % INPUTS:
    %       Fr: frame rate
    %       rawSpeed: un-processed walking speed data
    %       rawImage: Dalsa Image
    %
    % OUTPUTS:
    %       speedata: downsampled data used to plot on the front pannel
    %       speedaccBinary: Binarized accleration
    %       speedaveBinary: averaged binarized acceleration during one frame

%     timespan = numImage/Fr; % time
    timespan = floor(length(Speed)/Fs);
    
    Freal = floor(length(Speed)/timespan);
    
    % 10 Hz low pass filter
    [zeroa, poleb, gain] = butter(5,2/Fs*10, 'low');
    [sos,g] = zp2sos(zeroa,poleb, gain);
    Speedfilt = filtfilt(sos,g, Speed);
    Speeddata = arrayfun(@(x)mean(Speedfilt(round((x-1)*Fs/Fr+1 : x*Fs/Fr))), 1:numImage); % plot this one on the panel

    % --- Continue working on this filtered speed
    % compute accleration, in m/s/point
    Speedacc = [0; diff(Speedfilt)]; % in Bing's program, this step reduced speed data by 1 sample point

    % Binarize the acceleration
    if Fs == 10000
        accThreshold = 3e-6; % 3e-6 m/s/point, i,e, 0.03 m/s^2
    elseif Fs == 1000
        accThreshold = 3e-5; % 3e-5 m/s/point, i,e, 0.03 m/s^2
    else
        error('Wrong sample rate');
    end
    SpeedaccBinary = (abs(Speedacc)>accThreshold);

    % If the mouse running more than 10% time between consecutive image frames,
    % we define this frame as a running frame
    SpeedaveBinary = arrayfun(@(x)mean(SpeedaccBinary(round((x-1)*Fs/Fr+1 : x*Fs/Fr)))>0.1, 1:numImage); % logical array
    tmp = NaN*zeros(length(SpeedaveBinary),1);
    tmp(SpeedaveBinary) = 1;
    SpeedaveBinary = tmp; % numerical array
    clear tmp;
end