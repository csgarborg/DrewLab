close all

set1 = [.098 .091 .086 .078 .071 .065 .061 .058 .054 .05 .047];
set2 = [.094 .087 .079 .072 .067 .064 .061 .057 .054 .051 .048];
set3 = [.096 .09 .083 .076 .071 .067 .061 .056 .054 .051 .049];
set4 = [.22 .2 .187 .175 .167 .158 .15 .146 .139 .133 .129];
set5 = [.145 .125 .114 .103 .093 .084 .079 .072 .068 .065 .061];

mAVals = 0:10:100;

figure(1)
plot(mAVals,set1,'g');
hold on
plot(mAVals,set2,'r');
plot(mAVals,set3,'b');
plot(mAVals,set4,'k');
plot(mAVals,set5,'c');
hold off
xlabel('Current Input (mA)')
ylabel('Focal Length (m)')
title('Focal Length Range of ETL')
ylim([0 (ceil(max([set1 set2 set4])*10)/10)])
legend('Original 1','Original 2','Original 3','No offset lens','-10 D')

figure(2)
plot(mAVals,set1.^-1,'g');
hold on
plot(mAVals,set2.^-1,'r');
plot(mAVals,set3.^-1,'b');
plot(mAVals,set4.^-1,'k');
plot(mAVals,set5.^-1,'c');
hold off
xlabel('Current Input (mA)')
ylabel('Diopters (m^-1)')
title('Diopter Range of ETL')
ylim([0 ceil(max([1./set1 1./set2 1./set4]))])
legend('Original 1','Original 2','Original 3','No offset lens','-10 D')

figure(3)
plot(mAVals,set1,'k');
hold on
plot(mAVals,1./((1./set1)-10),'r');
plot(mAVals,1./((1./set1)-8),'k--');
plot(mAVals,1./((1./set1)-6),'b');
plot(mAVals,1./((1./set1)-4),'m');
plot(mAVals,1./((1./set1)-2),'b--');
hold off
xlabel('Current Input (mA)')
ylabel('Focal Length (m)')
title('Focal Length Range of ETL with Offset Lens')
legend('No Offset','10D','8D','6D','4D','2D')