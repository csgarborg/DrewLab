close all
a = {'001','002','003','004','005','006','007','008','009','010'};
for n = 1:size(a,2)
    h(n) = figure('Color','White');
    load(['sgvd04_211014_' a{n} 'a_pArt1_MScanData.mat']);
    subplot(5,1,1)
    t = (1:length(MScanData.data.vesselDiameter))*(1/MScanData.notes.frameRate);
    procData = smoothdata(MScanData.data.vesselDiameter,'gaussian',10);
    plot(t,procData)
    xlabel('Time (s)')
    ylabel('Diameter (\mum)')
    title('Surface Vessel (0 \mum depth)')
    load(['sgvd04_211014_' a{n} 'b_pArt1_MScanData.mat']);
    subplot(5,1,2)
    t = (1:length(MScanData.data.vesselDiameter))*(2/MScanData.notes.frameRate);
    procData = smoothdata(MScanData.data.vesselDiameter,'gaussian',10);
    plot(t,procData)
    xlabel('Time (s)')
    ylabel('Diameter (\mum)')
    title('Penetrating Vessel (75 \mum depth)')
    subplot(5,1,3)
    t = (1:length(MScanData.data.ball))*(1/MScanData.notes.analogSamplingRate);
    procData = smoothBallData([t' MScanData.data.ball],MScanData.notes.analogSamplingRate);
    procData(:,2) = abs(convertBallVoltToMPS(procData(:,2)));
    plot(procData(:,1),procData(:,2));
    xlabel('Time (s)')
    ylabel('Locomotion (m/s)')
    title('Ball Rotation')
    subplot(5,1,4)
    plot((1:length(MScanData.data.leftWhiskerStim))*(1/MScanData.notes.analogSamplingRate),MScanData.data.leftWhiskerStim)
    xlabel('Time (s)')
    ylabel('Voltage (V)')
    title('Left Whisker Stim')
    subplot(5,1,5)
    plot((1:length(MScanData.data.rightWhiskerStim))*(1/MScanData.notes.analogSamplingRate),MScanData.data.rightWhiskerStim)
    xlabel('Time (s)')
    ylabel('Voltage (V)')
    title('Right Whisker Stim')
end

savefig(h,'H:\21-10-11_MouseExp\211011_vesselDiameter.fig');