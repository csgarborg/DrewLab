close all
a = {'005','006'};
dataMatrix = {};
for n = 1:size(a,2)
%     h(n) = figure('Color','White');
    load(['sgvd04_211011_' a{n} 'a_pArt1_MScanData.mat']);
    t1 = (1:length(MScanData.data.vesselDiameter))*(1/MScanData.notes.frameRate);
    procData1 = MScanData.data.vesselDiameter;
    procData1 = medfilt1(procData1,13);
    procData1 = smoothdata(procData1,'gaussian',10);
    
    plot(t1,procData1)
    title('Select start and stop time (t) to use to determine mean baseline')
    tVals = ginput(2);
    close
    i = tVals(1,1) <= t1 & t1 <= tVals(2,1);
    baseline = mean(procData1(i));
    procData1 = 100*((procData1 - baseline)./(baseline));

    
    
    load(['sgvd04_211011_' a{n} 'b_pArt1_MScanData.mat']);
    t2 = (1:length(MScanData.data.vesselDiameter))*(2/MScanData.notes.frameRate);
    procData2 = MScanData.data.vesselDiameter;
    procData2 = medfilt1(procData2,13);
    procData2 = smoothdata(procData2,'gaussian',10);
    
    plot(t2,procData2)
    title('Select start and stop time (t) to use to determine mean baseline')
    tVals = ginput(2);
    close
    i = tVals(1,1) <= t2 & t2 <= tVals(2,1);
    baseline = mean(procData2(i));
    procData2 = 100*((procData2 - baseline)./(baseline));
    
    dataMatrix(end+1:end+4) = {t1;procData1;t2;procData2};
end

for n = 1:size(a,2)
    h(n) = figure('Color','White');
    plot(dataMatrix{(n-1)*4+1},dataMatrix{(n-1)*4+2})
    xlabel('Time (s)')
    ylabel('Percent Change in Diameter')
    if n == 1
        title('Percent Change of Vessel Diameter (0 and 75 \mum depth)')
    else
        title('Percent Change of Vessel Diameter (0 and 50 \mum depth)')
    end
    
    hold on
    
    plot(dataMatrix{(n-1)*4+3},dataMatrix{(n-1)*4+4})
    
    legend('Surface','Penetrating')
    
    hold off
end

savefig(h,'H:\21-10-11_MouseExp\211011_vesselDiameterPercentChange.fig');