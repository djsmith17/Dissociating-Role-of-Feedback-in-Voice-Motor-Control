% Requires the Data Acquisition Toolbox (Mathworks)

dev = 'Dev2';

s = daq.createSession('ni');
addAnalogOutputChannel(s, dev, 0, 'Voltage');

s.Rate = 8000;

outputSignal1 = [zeros(1,6000) 3*ones(1,12000) zeros(1,6000)]';

plot(outputSignal1);
xlabel('Time'); ylabel('Voltage'); legend('Analog Output 0');

queueOutputData(s, outputSignal1);
s.startForeground;