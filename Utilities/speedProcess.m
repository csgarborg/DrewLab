function [speeddata, speedaccBinary, speedaveBinary] = speedProcess(Fr, rawSpeed, rawImage)
    % INPUTS:
    %       Fr: frame rate
    %       rawSpeed: un-processed walking speed data
    %       rawImage: Dalsa Image
    %
    % OUTPUTS:
    %       speedata: downsampled data used to plot on the front pannel
    %       speedaccBinary: Binarized accleration
    %       speedaveBinary: averaged binarized acceleration during one frame
    

 


        % --- Match length --- %       
        Fs = 10000; % sample freq, 20 kHz
        
%         % Method 1: Bing's Method
%         p = size(rawImage,3);
%         q = floor(length(rawSpeed)/(Fs/Fr));
%         speedresamp = resample(rawSpeed, p, q);
%         speedtrim = speedresamp(1:Fs*p/Fr); 
        
         % Method 2: Qingguang's Modification
        numImage = rawImage; % how many frames in Dalsa image
        
        timespan = numImage/Fr;  

        Freal = floor(length(rawSpeed)/timespan);        
        % resample raw speed data to match length
        speedresamp = resample(rawSpeed, Fs, Freal);                 
   
        % trim data
        speedtrim = speedresamp(1:Fs*timespan);
        % 10 Hz low pass filter
        [zeroa, poleb, gain] = butter(5,2/Fs*10, 'low');
        [sos,g] = zp2sos(zeroa,poleb, gain);
        speedfilt = filtfilt(sos,g, speedtrim);
        speeddata = arrayfun(@(x)mean(speedfilt(round((x-1)*Fs/Fr+1) : min(round(x*Fs/Fr),length(speedfilt)))), 1:numImage); % plot this one on the panel

        % --- Continue working on this filtered speed
        % compute accleration, in m/s2
        speedacc = [0; diff(speedfilt)]; % in Bing's program, this step reduced speed data by 1 sample point
        
        % Binarize the acceleration
        accThreshold = 1e-6;
        speedaccBinary = (abs(speedacc)>accThreshold);
        
        % If the mouse running more than 10% time between consecutive image frames,
        % we define this frame as a running frame        
        speedaveBinary = arrayfun(@(x)mean(speedaccBinary(round((x-1)*Fs/Fr+1) : min(round(x*Fs/Fr),length(speedfilt))))>0.1, 1:numImage); % logical array
        tmp = NaN*zeros(length(speedaveBinary),1);
        tmp(speedaveBinary) = 1;
        speedaveBinary = tmp; % numerical array
        clear tmp;  
end