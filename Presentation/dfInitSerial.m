function [ard,ardCh] = dfInitSerial(comPort,baudRate,trialLen)
% function [ard,ardCh] = dfInitSerial(comPort,baudRate,trialLen)
%
% configures comPort as a serial port for reading/writing data at the
% specified baudRate.
%
% INPUT
% comPort:      string. Port number wihere arduino is connected
% baudRate:     integer. serial port baudrate. Must match the one of the device
%               and of the dfined port
% trialLen:     integer. Trial length in seconds
%
% OUTPUT
% ard:          Serial port object
% ardCh:     arduino channel names (Input/Output) (was niCh)
%
% Andres    :   v1  : init. 11 April 2018

% Create arduino Object
ard = serial(comPort,'baudrate',baudRate);

% Configure
ard.BytesAvailableFcnCount = 40;
ard.BytesAvailableFcnMode = 'byte';
ard.BytesAvailableFcn = @instrcallback;

ard.BytesAvailableFcnMode

% arduino configure info
ardCh.Rate = 8000;
ardCh.DurationInSeconds = trialLen;

% Open ard object port
fopen(ard);

% Verbose
fprintf('%s\nOpenned serial port %s with baudrate %i\n%s\n',...
    '==============================================',...
    comPort,baudRate,...
    '==============================================')

end
