% clear all; close all; clc;

s = daq.createSession('ni');
addAnalogOutputChannel(s,'Dev3',0,'Voltage');

s.Rate = 8000;
outputSingleValue = 2;
outputSingleScan(s,outputSingleValue);

outputSignal1 = [zeros(1,6000) 3*ones(1,12000) zeros(1,6000)]';

plot(outputSignal1);
xlabel('Time'); ylabel('Voltage'); legend('Analog Output 0');

queueOutputData(s, outputSignal1);
s.startForeground;