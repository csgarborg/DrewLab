function powerSpectrumPlotSimple(data,sampleRate,tapers,fpass,err)
% example: powerSpectrumPlotSimple(data,10000,[2 3],[0.025 100],[2 0.05])
dt = 1/sampleRate;
params.Fs = 1/dt;
params.tapers = tapers;
params.fpass = fpass;
params.err = err; %use jack-knife resampling confidence intervals p = 0.05
% [Power, Hz, Error] = mtspectrumc(rawEKG',params);
[Power, Hz, Error] = mtspectrumc(data,params);
semilogy(Hz,Power,'k')
end