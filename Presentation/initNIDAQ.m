function s = initNIDAQ
%This function sets up the daq object for the controlling parts of the
%experiment. The sampling rate is currently hard set at 8000hz, but may
%eventually become an input variable. 
%Check this website if things go wrong: 
%https://www.mathworks.com/hardware-support/nidaqmx.html

dev = 'Dev3';
s = daq.createSession('ni');
addAnalogOutputChannel(s, dev, 0, 'Voltage'); %Output to the Perturbatron
addAnalogOutputChannel(s, dev, 1, 'Voltage'); %Output voltage to Force Sensor
addAnalogInputChannel(s, dev, 0, 'Voltage'); %Input from Force Sensor 1
addAnalogInputChannel(s, dev, 1, 'Voltage'); %Input from Force Sensor 1

s.Rate = 8000;
end