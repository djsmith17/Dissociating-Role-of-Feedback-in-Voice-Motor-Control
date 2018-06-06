function [dataOUT,timeOUT,nextStep] = dfArduinoSerial(dataOUT,timeOUT,ard,trialTIC,nextStep)
%
%
%
%
%
%
%
% Andres    :   v1  : init. 11 April 2018

%% read everything each step, and save time for it
timeOUT.press(nextStep) = toc(trialTIC);
dataOUT.press(nextStep)= fscanf(ard,'%f');        % digital perturbatron out signal back in
%fscanf(ard,'%f',10)

nextStep = nextStep + 1;
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS START HERE %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testSerial 
%
%
%
%
%

clear all testTime testData

ard = serial('COM19','baudrate',115200);

% Configure
ard.BytesAvailableFcnCount = 512;
ard.BytesAvailableFcnMode = 'byte';
% ard.BytesAvailableFcn = @instrcallback;
ard.BytesAvailableFcnMode

% Open ard object port
fopen(ard);

% Read data
nTrials = 1000;
testTime = nan(nTrials,1);
testData = nan(nTrials,1);
trialStartTic = tic;

for iTrial = 1:nTrials, iTrial
    testTime(iTrial) = toc(trialStartTic);
    %temp = fscanf(ard,'%i')
    temp = fread(ard,1,'uint8')
    %temp = fread(ard)
    testData(iTrial) = temp; pause(0.001)
end
fclose(ard)



if ~isempty(instrfind)
    fclose(instrfind);
    delete(instrfind);
end

for ii = 1:nTrials
    testData{iTrial}
end



end
