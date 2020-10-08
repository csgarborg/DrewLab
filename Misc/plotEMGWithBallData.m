function plotEMGWithBallData(filename)

ballData = load(filename);

t = ballData(:,1);
rot = ballData(:,2);
emg = ballData(:,3);

h(1) = figure('Color','White');
x1 = subplot(2,1,1);
plot(t,emg,'k')
title('\fontsize{20pt}\bf{Abdominal EMG Recording}')
xlabel('Time (s)')
ylabel('Volts')
grid on
axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) -1 ceil(max(movementData.ballData(:,2)))])


x2 = subplot(2,1,2);
plot(t,rot,'k')
title('\fontsize{20pt}\bf{Ball Movement}')
xlabel('Time (s)')
ylabel('Movement')
grid on
axis([min(movementData.ballData(:,1)) max(movementData.ballData(:,1)) -1 ceil(max(movementData.ballData(:,2)))])
end