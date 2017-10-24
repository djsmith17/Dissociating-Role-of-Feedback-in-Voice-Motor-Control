function [s, sd, niCh, nVS] = initNIDAQSess(dev, trialLen)
%This function sets up the daq object for the controlling parts of the
%experiment, most specifically the Perturbatron. 
%The sampling rate is currently hard set at 8000hz, but may eventually
%become an input variable. 
%
%Requires the Data Acquisition Toolbox
%Check this website if things go wrong: 
%https://www.mathworks.com/hardware-support/nidaqmx.html

%inputs:
%dev:      The device to use. Defaults to Dev2
%trialLen: Length of the recording. This is used to create a -1V signal

%outputs:
%s:        The NIDAQ session for recording
%niCh:     NIDAQ Channels names (Input/Output)
%nVS:      Negative Voltage Source. A -1V signal needed for some circuits. 

if isempty(dev)
    dev = 'Dev2';
end

s      = daq.createSession('ni');
s.Rate = 8000;
s.IsContinuous = true;
s.DurationInSeconds = trialLen;

% A0 Channels
addAnalogOutputChannel(s, dev, 0, 'Voltage'); %Output signal to the Perturbatron
addAnalogOutputChannel(s, dev, 1, 'Voltage'); %Output voltage to power the Force Sensor

% AI Channels
addAnalogInputChannel(s, dev, 0, 'Voltage'); %Input signal from NIDAQ Perturbatron Signal
addAnalogInputChannel(s, dev, 1, 'Voltage'); %Input signal from Force Sensor 1: Collar
addAnalogInputChannel(s, dev, 2, 'Voltage'); %Input signal from Force Sensor 2: Neck
addAnalogInputChannel(s, dev, 3, 'Voltage'); %Input signal from Pressure Sensor
addAnalogInputChannel(s, dev, 4, 'Voltage'); %Input signal from Microphone 
addAnalogInputChannel(s, dev, 5, 'Voltage'); %Input signal from Headphones
addAnalogInputChannel(s, dev, 6, 'Voltage'); %Input signal from Optical Triggerbox

%Document which channels are which. Manual hardset at the moment. When more
%versions are needed, this will be converted to a switch case. 
niCh.ao0 = 'Perturbatron';
niCh.ao1 = 'Negative Voltage Source';
niCh.ai0 = 'Perturbatron';
niCh.ai1 = 'Force Sensor: Collar';
niCh.ai2 = 'Force Sensor: Neck';
niCh.ai3 = 'Pressure Sensor';
niCh.ai4 = 'Microphone';
niCh.ai5 = 'Headphones';
niCh.ai6 = 'Optical Triggerbox';

addTriggerConnection(s, 'External', 'Dev2/PFI0', 'StartTrigger')

nVS = zeros(s.Rate*trialLen, 1) - 1;
nVS(1) = 0; nVS(end) = 0;

lrec  = addlistener(s,'DataAvailable', ...
                  @(src,event) updateNIDAQdata(sd, event.TimeStamps, event.Data));

lsend = addlistener(s,'DataRequired', ...
                  @(src,event) src.queueOutputData(NIDAQsig)); % @myFunction(src,evt,NIDAQsig)


sd = perturbDAQClass;

disp('NIDAQ has been initialized!')
end