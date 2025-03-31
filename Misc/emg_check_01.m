close all;
clear all
x=load('C:\Workspace\210323_004.txt');
%%

emg_raw=x(:,3);
emg=filloutliers(emg_raw, 'center','mean', 'ThresholdFactor', 4);%this takes out any large artifacts

params.Fs=10000;
the_time=(1:length(emg))/10000;
params.tapers=[5 9];
[S,f]=mtspectrumc(emg-mean(emg),params);
df=f(2)-f(2);

% band pass filter emg between 300-3kHz
[b,a]=butter(5,[300]/5000,'high');
[b1,a1]=butter(5,[3000]/5000,'low');
emg_band=filtfilt(b,a,emg);
emg_band=filtfilt(b1,a1,emg_band);

% band pass filter ekg between 70-200Hz
[bk,ak]=butter(5,[70]/5000,'high');
[bk1,ak1]=butter(5,[200]/5000,'low');
ekg_band=filtfilt(bk,ak,emg);
ekg_band=filtfilt(bk1,ak1,ekg_band);
[S1,f1]=mtspectrumc(emg_band,params);

figure(100)
hold off
loglog(f,S)
hold on
loglog(f1,S1)
%%
figure(101)
subplot(211)
hold off
plot(the_time,emg-mean(emg),'k')
hold on;
plot(the_time,ekg_band,'r')

plot(the_time,emg_band,'b')
legend('raw','EKG','bandpassed EMG')

title('filtered EMG')
%% plot the envelops,, 200 Hz low pass filtered version of power
[b2,a2]=butter(5,[300]/5000,'low');
raw_pow=filtfilt(b2,a2,((emg-mean(emg)).^2));
filt_pow=filtfilt(b2,a2,((emg_band).^2));
ekg_pow=filtfilt(b2,a2,((ekg_band).^2));
figure(101)
subplot(212)
hold off
plot(the_time,raw_pow/mean(raw_pow),'k');
hold on
plot(the_time,ekg_pow/mean(ekg_pow),'r')
plot(the_time,filt_pow/mean(filt_pow),'b')

legend('raw','EKG','bandpassed EMG')
xlabel('time, s')
ylabel('norm EMG power')

