function [dataOUT,timeOUT,nextStep] = dfArduinoGetData(dataOUT,timeOUT,ard,ardCh,trialTIC,nextStep)
% function [dataOUT,timeOUT,nextStep] = dfArduinoGetData(dataOUT,timeOUT,ard,ardCh,trialTIC,nextStep)
%
% This function reads out digital and analog inputs from the Arduino using the 
% 'ard' object, for the specified 'ardCh' channels and ap[pends the data
% and time staps to 'dataOUT' and 'timeOUT' respectively.
%
% INPUT
% dataOUT:      matrix. [nINPUTs,nSamps] Values for digital and analog inputs 
%               for perturbatron in, photodiode, pressure, temperature,
%               microphone, audio signals
% timeOUT:      matrix. [nINPUTs,nSamps] Time values during the trial (in seconds) 
%               using trialTIC = zero = trial start time. 
% ard:          Arduino object
% ardCh:        structure. Has names and channels/pins for the specified (ard)
%               arduino. 
% trialTIC:     unit64 values. Output from tic command established as the
%               trial start time
% nextStep:     integer. Time step number
%
% OUTPUT
% dataOUT:      matrix. [nINPUTs,nSamps] Updated values for digital and analog inputs 
%               for perturbatron in, photodiode, pressure, temperature,
%               microphone, audio signals
% timeOUT:      matrix. [nINPUTs,nSamps] Updated time values during the trial (in seconds) 
%               using trialTIC = zero = trial start time. 
%
% --------------------------------------------------------------------------
% You first need to configure the arduino using the command 'dfInitArduino.m'
% % This is the listof pins currently configure by 'dfInitArduino.m'
%
% ardCh.do2 = 'Perturbatron_out';         % digital output to trigger the Perturbatron hardware
% ardCh.do3 = 'trialON';               % digital output --> Code for event = trial start
% ardCh.do4 = 'trialInfo';                 % digital output --> Code for event = trial end
% ardCh.di5 = 'Perturbatron_in';          % digital input from Perturbatron trigger output
% ardCh.di6 = 'Optical Triggerbox';       % digital input from optical photodiode (on screen stim time)
% ardCh.ai0 = 'Pressure Sensor';          % analog input pressure sensor
% ardCh.ai1 = 'Temperature';              % analog input temperature sensor
% ardCh.ai2 = 'Microphone';               % analog input user's microphone
% ardCh.ai3 = 'Headphones';               % analog input user's headphone
% 
% % Proper pins assignment for arduino
%
% ardCh.pertOut   = 'D2';
% ardCh.trialON = 'D3';
% ardCh.trialInfo = 'D4';
% ardCh.pertIn    = 'D5';
% ardCh.optTrig   = 'D6';
% ardCh.press     = 'A0';
% ardCh.temp      = 'A1';
% ardCh.micropho  = 'A2';
% ardCh.headpho   = 'A3';
% --------------------------------------------------------------------------
%
% Andres    :   v1  : init. 05 April 2018

%% read everything each step, and save time for it
timeOUT.pertIn(nextStep) = toc(trialTIC);
dataOUT.pertIn(nextStep)= readDigitalPin(ard,ardCh.pertIn);        % digital perturbatron out signal back in
timeOUT.optTrig(nextStep) = toc(trialTIC);
dataOUT.optTrig(nextStep) = readDigitalPin(ard,ardCh.optTrig);     % digital optical trigger (screen photoDiode)

% timeOUT.press(nextStep) = toc(trialTIC);
% dataOUT.press(nextStep) = readVoltage(ard,ardCh.press);              % pressure analog

nextStep = nextStep + 1;

% timeOUT.pertIn(end+1) = toc(trialTIC);
% dataOUT.pertIn(end+1)= readDigitalPin(ard,ardCh.pertIn);        % digital perturbatron out signal back in
% 
% timeOUT.optTrig(end+1) = toc(trialTIC);
% dataOUT.optTrig(end+1) = readDigitalPin(ard,ardCh.optTrig);     % digital optical trigger (screen photoDiode)
% 
% timeOUT.press(end+1) = toc(trialTIC);
% dataOUT.press(end+1) = readVoltage(ard,ardCh.press);              % pressure analog

% 
% timeOUT.pertIn = cat(2,timeOUT.pertIn, toc(trialTIC));
% dataOUT.pertIn = [dataOUT.pertIn, readDigitalPin(ard,ardCh.pertIn)];        % digital perturbatron out signal back in
% 
% timeOUT.optTrig = [timeOUT.optTrig, toc(trialTIC)];
% dataOUT.optTrig = [dataOUT.optTrig, readDigitalPin(ard,ardCh.optTrig)];     % digital optical trigger (screen photoDiode)
% 
% timeOUT.press = [timeOUT.press, toc(trialTIC)];
% dataOUT.press = [dataOUT.press, readVoltage(ard,ardCh.press)];              % pressure analog
% 
% timeOUT.temp = [timeOUT.temp, toc(trialTIC)];                       
% dataOUT.temp = [dataOUT.temp, readVoltage(ard,ardCh.temp)];                 % temperature analog
% 
% timeOUT.micropho = [timeOUT.micropho, toc(trialTIC)];
% dataOUT.micropho = [dataOUT.micropho, readVoltage(ard,ardCh.micropho)];     % microphone analog
% 
% timeOUT.headpho = [timeOUT.headpho, toc(trialTIC)];
% dataOUT.headpho = [dataOUT.headpho, readVoltage(ard,ardCh.headpho)];        % headphone analog


end

