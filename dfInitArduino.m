function [ard, ardCh] = dfInitArduino(comPort,ardType,trialLen)
% This function sets up the Arduino (Uno by default) for controlling parts of the
% experiment, most specifically the Perturbatron, the pressure and temperature 
% sensing as well as extra trial oputputs useful for online and offline analysis 
% (trial start/end triggers, etc). The sampling rate is currently hardcoded to 
% meet the Arduino clock, at 8000Hz. The compatibility with the flexiforce
% sensors was removed since no negative voltage is input by the Uno.
% 
% Requires the Arduino support for Matlab, file 'arduinoio.mlpkginstall'. 
% For more info go to: 
% https://www.mathworks.com/hardware-support/arduino-matlab.html
% https://www.mathworks.com/help/supportpkg/arduinoio/ug/configure-setup-for-arduino-hardware.html#bvn8vm8
% 
% INPUTS:
% trialLen:  Length of the recording. This is used to create a -1V signal
% comPort:   string. The COM (serial) port where the arduino is connected. Check 'Devices and Printers' 
%            on the 'Control panel' in Windows if not sure. Example: 'COM4'
% ardType:   string. Device name, can be 'MEGA', 'Uno'
% 
% OUTPUTS:
% ard:       The arduino object used for recording (was s)
% ardCh:     arduino channel names (Input/Output) (was niCh)
% 
%% Common commands after configuring arduino (using the 'ard = arduino(comPort,ardType);' command)
% % Read analog voltage
% readVoltage(ard,'A0');
% % Write digital pin out
% writeDigitalPin(ard,'D2',1)
% % Read digital pin in
% readDigitalPin(ard,'D3');
% 
% Valid pin numbers for Arduino Uno are 
% % D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13,  (D0 and D1 are not usefull for regular digital read/write tasks)
% % A0, A1, A2, A3, A4, A5.
%
%%
%
% Dante     :   v1  : init. 2015-2018 Using ad reference code dfInitNIDAQ.m
% Andres    :   v2  : changed code to be compatible with Arduino Uno. 05 April 2018


% AFSG 20180405
% if isempty(dev)
%     dev = 'Dev2';
% end
% 
% s = daq.createSession('ni');
% 
% s.Rate = 8000;
% s.DurationInSeconds = trialLen;
% 
% addAnalogOutputChannel(s, dev, 0, 'Voltage'); %Output signal to the Perturbatron
% addAnalogOutputChannel(s, dev, 1, 'Voltage'); %Output voltage to power the Force Sensor
% 
% addAnalogInputChannel(s, dev, 0, 'Voltage'); %Input signal from NIDAQ Perturbatron Signal
% addAnalogInputChannel(s, dev, 1, 'Voltage'); %Input signal from Force Sensor 1: Collar
% addAnalogInputChannel(s, dev, 2, 'Voltage'); %Input signal from Force Sensor 2: Neck
% addAnalogInputChannel(s, dev, 3, 'Voltage'); %Input signal from Pressure Sensor
% addAnalogInputChannel(s, dev, 4, 'Voltage'); %Input signal from Microphone 
% addAnalogInputChannel(s, dev, 5, 'Voltage'); %Input signal from Headphones
% addAnalogInputChannel(s, dev, 6, 'Voltage'); %Input signal from Optical Triggerbox
% 
% %Document which channels are which. Manual hardset at the moment. When more
% %versions are needed, this will be converted to a switch case. 
% 
% niCh.ao0 = 'Perturbatron';
% niCh.ao1 = 'Negative Voltage Source';
% niCh.ai1 = 'Force Sensor: Collar';
% niCh.ai2 = 'Force Sensor: Neck';
% niCh.ai3 = 'Pressure Sensor';
% niCh.ai4 = 'Microphone';
% niCh.ai5 = 'Headphones';
% niCh.ai6 = 'Optical Triggerbox';

% Default inputs

%if nargin < 1, comPort= 'COM19'; end
if nargin < 1, comPort= 'COM4'; end
if nargin < 2, ardType = 'Uno'; end

% Setup Arduino
ard = arduino(comPort,ardType);         % to close just clear the arduino object: 'clear ard'

% arduino configure info
ardCh.Rate = 8000;
ardCh.DurationInSeconds = trialLen;

% Name of I/O ports
ardCh.do2 = 'Perturbatron_out';         % digital output to trigger the Perturbatron hardware
% ardCh.do3 = 'trialStart';               % digital output --> Code for event = trial start
% ardCh.do4 = 'trialEnd';                 % digital output --> Code for event = trial end
ardCh.do3 = 'trialON';               % digital output --> Code for event = trial start
ardCh.do4 = 'trialInfo';                 % digital output --> Code for event = trial end
ardCh.di5 = 'Perturbatron_in';          % digital input from Perturbatron trigger output
ardCh.di6 = 'Optical Triggerbox';       % digital input from optical photodiode (on screen stim time)
ardCh.ai0 = 'Pressure Sensor';          % analog input pressure sensor
ardCh.ai1 = 'Temperature';              % analog input temperature sensor
ardCh.ai2 = 'Microphone';               % analog input user's microphone
ardCh.ai3 = 'Headphones';               % analog input user's headphone

% Proper pins assignment for arduino
ardCh.pertOut   = 'D2';
ardCh.trialON   = 'D3';
ardCh.trialInfo = 'D4';
% ardCh.trialStart = 'D3';
% ardCh.trialEnd = 'D4';
ardCh.pertIn    = 'D5';
ardCh.optTrig   = 'D6';
ardCh.press     = 'A0';
ardCh.temp      = 'A1';
ardCh.micropho  = 'A2';
ardCh.headpho   = 'A3';
   
% Everything is ready!
fprintf('%s\n%s\n%s\n','==================================',...
    'Arduino Uno has been initialized!','==================================')

end