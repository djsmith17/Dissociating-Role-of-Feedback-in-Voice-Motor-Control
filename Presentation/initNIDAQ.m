function s = initNIDAQ
%This function sets up the daq object for the controlling parts of the
%experiment. The sampling rate is currently hard set at 8000hz, but may
%eventually become an input variable.

s = daq.createSession('ni');
addAnalogOutputChannel(s,'Dev3',0,'Voltage'); %Output to the Perturbatron
addAnalogOutputChannel(s,'Dev3',1,'Voltage'); %Output voltage to Force Sensor
addAnalogInputChannel(s,'Dev3',0,'Voltage'); %Input from Force Sensor 1
addAnalogInputChannel(s,'Dev3',1,'Voltage'); %Input from Force Sensor 1

s.Rate = 8000;
end