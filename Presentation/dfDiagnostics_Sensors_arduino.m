function dfDiagnostics_Sensors_arduino()
% dfDiagnostics_Sensors() collects data from the NIDAQ in the way that the
% main experimental setup. A quick test of the force sensors before running 
% the actual experiment. This makes sure that the sensors are working the way 
% they should be and we can continue with the experiment. Eventually this will 
% also include the pressure sensor. 
%
% This script calls the following (4) functions:
% -dfDirs.m
% -dfInitArduino.m
% -dfSetTrialOrder.m
% -dfMakePertSignal.m
%
% This script includes the following (2) subfunctions:
% -initLiveResult
% -updateLiveResult
%
% This script requires the following toolboxes:
% -Arduino Support from MATLAB 
% https://www.mathworks.com/hardware-support/arduino-matlab.html
%
% Dante     :   init    : 2015-2018. Created body of the code compatible with NI-card. 
% Andres    :   v2      : added compatibility with Arduino Uno to replace NI/card 
%                        (after installing arduino driver and toolbox in Matlab 2017b. 
%                         See installation doc). 05 April 2018

close all;
ArdOrSerial = 'arduino';

% Main Experimental prompt: Subject/Run Information
prompt = {'Subject ID:',...
          'Session ID:',...
          'Number of Trials:',...
          'Percent Perturbed (Dec)',...
          'Balloon:', ...
          'Collect New Data?:'};
name = 'Subject Information';
numlines = 1;
%defaultanswer = {'null', 'DS1', '10', '0.5', '2.0K_4','yes'};
defaultanswer = {'null', 'DS1', '1', '1', '2.0K_4','yes'};
ExpPrompt = inputdlg(prompt, name, numlines, defaultanswer);

if isempty([ExpPrompt{1},ExpPrompt{2},ExpPrompt{3},ExpPrompt{4},ExpPrompt{5},ExpPrompt{6}])
    msgbox('Error, please input values in the cells!')
    return
end

%Experiment Configurations
%expParam.project       = 'NIDAQSensorDiagnostics';         % Changed AFSG 20180405
expParam.project        = 'ArdUno_SensorDiagn';
expParam.expType        = 'Somatosensory Perturbation_Perceptual';
expParam.subject        = ExpPrompt{1}; %Subject#, Pilot#, null
expParam.run            = ExpPrompt{2};
expParam.curSess        = [expParam.subject ' ' expParam.run];
expParam.gender         = 'N/A';

% expParam.niDev         = 'Dev2';                      % NIDAQ Device Name. For more information, see dfInitNIDAQ
expParam.dev            = 'ArdUno';                     % Arduino Uno Device Name. For more information, see dfInitArdUNO. % Changed AFSG 20180405
expParam.comPort        = 'COM8'; %COM18                      % Port wehre arduino is connected. In Windows, Check 'Devices and Printers' in 'Control Panel'
expParam.ardType        = 'Uno';                        % type of Arduino. Useful if later Arduino MEGA or other is used.
expParam.baudRate       = 115200;                       % Double check baudrate with drive properties
expParam.trialLen       = 4;                            % Seconds
expParam.numTrial       = str2double(ExpPrompt{3});
expParam.curTrial       = [];
expParam.perCatch       = str2double(ExpPrompt{4});
expParam.balloon        = ExpPrompt{5};
expParam.AudFB          = 'Masking Noise';
expParam.AudFBSw        = 2;
expParam.sRateAnal      = 16000;
expParam.resPause       = 3;
expParam.trialLenLong   = expParam.numTrial*(expParam.trialLen + expParam.resPause);
expParam.sigLong        = [];

% Logic for recording data AFSG 20180405
sv2F                    = 1; %Boolean
collectNewData          = ExpPrompt{6};

%expParam.project %={'NIDAQSensorDiagnostics','ArduinoUnoSensorDiagnostics'};  if expParam.dev = 'ArdUno', expParam.niDev = ''; end; if expParam.dev = 'Dev2', expParam.niDev = 'Dev2', end
    
%Set our dirs based on the project
dirs = dfDirs(expParam.project);

dirs.RecFileDir    = fullfile(dirs.RecData, expParam.subject, expParam.run); % Where to save data
dirs.SavResultsDir = fullfile(dirs.RecData, expParam.subject, expParam.run); % Where to save results

if exist(dirs.RecFileDir, 'dir') == 0
    mkdir(dirs.RecFileDir)
end
if exist(dirs.SavResultsDir, 'dir') == 0
    mkdir(dirs.SavResultsDir)
end
dirs.RecFile   = fullfile(dirs.RecFileDir, [expParam.subject expParam.run 'NSD.mat']);

if strcmp(collectNewData, 'yes')

    %Set up Parameters to control NIDAQ and Perturbatron
    switch ArdOrSerial
        case 'arduino'
            [ard,ardCh] = dfInitArduino(expParam.comPort,expParam.ardType,expParam.trialLen);
        case 'serial'
            [ard,ardCh] = dfInitSerial(expParam.comPort,expParam.baudRate,expParam.trialLen);
    end
    
    expParam.sRateQ = ardCh.Rate;   % arduino dummy sampling rate (sampling rate will depend on serial communication times)
    expParam.niCh   = ardCh;        % Structure of Channel Names
    
    % Set up the order of trials (Order of perturbed, control, etc)
    expParam.trialType = dfSetTrialOrder(expParam.numTrial, expParam.perCatch);
    
    % Select the trigger points for perturbation onset and offset and creating
    [expParam.sigs, expParam.trigs] = dfMakePertSignal(expParam.trialLen, expParam.numTrial, expParam.sRateQ, expParam.sRateAnal, expParam.trialType, 1);  
    
    % Convert perturbation/output signal into time start and time end for using Arduino
    trigTimes = squeeze(expParam.trigs(:,:,1));         % in the form [start,end], in seconds
    
    % Plot 'live' pressure values
    presH = initLiveResult(expParam, 2);
    
    % Pre-allocate memory for data collected
    % Experiment Data and Time vals from recording
    %     DAQin = []; DAQtime = [];
    %     pltStr = [];
    
    %     dataOUT = struct('pertIn',[],'optTrig',[],'press',[],'temp',[],'micropho',[],'headpho',[]);     % data matrix (indices match timeOUT indices)
    %     timeOUT = struct('pertIn',[],'optTrig',[],'press',[],'temp',[],'micropho',[],'headpho',[]);     % time matrix (indices match dataOUT indices)
    
    %521 --> all. increase file
    %878 --> no micro, no audio, no temp. increase file 
    
    nSamps = 100*expParam.numTrial;
    dataOUT = struct('pertIn',[1,nSamps],'optTrig',[1,nSamps],'press',[1,nSamps],'temp',[1,nSamps],'micropho',[1,nSamps],'headpho',[1,nSamps]);     % data matrix (indices match timeOUT indices)
    timeOUT = struct('pertIn',[1,nSamps],'optTrig',[1,nSamps],'press',[1,nSamps],'temp',[1,nSamps],'micropho',[1,nSamps],'headpho',[1,nSamps]);     % time matrix (indices match dataOUT indices)

    % Inititalize arduino to zeros. Set to zero all digital pins out
    switch ArdOrSerial
        case 'arduino'
            writeDigitalPin(ard,ardCh.pertOut,0)
            writeDigitalPin(ard,ardCh.trialON,0)
            writeDigitalPin(ard,ardCh.trialInfo,0)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Perform perturbation and data collection per trial
    LgdStruct = []; % Live-Result iterable Legend properties structure
    nextStep = 1;
    valvCount = 0;
    for ii = 1:expParam.numTrial
        
        expParam.curTrial = ii;
        alreadyStim = 0;
        fprintf('Running Trial %d and perturbed %i\n', ii,expParam.trialType(ii))
        
        % Perturb only if it is a perturbation trial
        if expParam.trialType(ii)
            
            if mod(valvCount,2) == 0; valvON = 1; else valvON = 0; end %#ok<*SEPEX>
            disp(valvON)
            valvCount = valvCount + 1;

            % Set trial start pin to high
            switch ArdOrSerial
                % case 'arduino', writeDigitalPin(ard,ardCh.trialON,randValveVal)
                case 'arduino', writeDigitalPin(ard,ardCh.trialON,valvON)
            end
            
            % first start trial at zero
            currTrigTime = 0;
            trialStartTic = tic;
            
            %% read everything each step, and save time for it            
            switch ArdOrSerial
                case 'arduino'
                    [dataOUT,timeOUT,nextStep] = dfArduinoGetData(dataOUT,timeOUT,ard,ardCh,trialStartTic,nextStep);
                case 'serial'
                    [dataOUT,timeOUT,nextStep] = dfArduinoSerial(dataOUT,timeOUT,ard,trialStartTic,nextStep);
            end
            
            % Perform analysis for the expParam.trialLen set time
            while currTrigTime <= expParam.trialLen

                % If current trial time is before trigStart or after trigEnd
                if or(currTrigTime <= trigTimes(ii,1), currTrigTime >= trigTimes(ii,2))
                    
                    % Set perturbatron out to zero
                    switch ArdOrSerial
                        case 'arduino', writeDigitalPin(ard,ardCh.pertOut,0), outVal = 0;
                    end
                    % Read all data
                    switch ArdOrSerial
                        case 'arduino'
                            [dataOUT,timeOUT,nextStep] = dfArduinoGetData(dataOUT,timeOUT,ard,ardCh,trialStartTic,nextStep);
                        case 'serial'
                            [dataOUT,timeOUT,nextStep] = dfArduinoSerial(dataOUT,timeOUT,ard,trialStartTic,nextStep);
                    end
      
                else   % if currTrigTime is between trialStat and TrialEnd, then pertOut = 1
                    
                    % Set perturbatron out to one
                    if ~alreadyStim
                        switch ArdOrSerial
                            case 'arduino', writeDigitalPin(ard,ardCh.pertOut,1),
                        end
                        alreadyStim = 1; outVal = 1;
                    end
                    
                    % Read all data
                    switch ArdOrSerial
                        case 'arduino'
                            [dataOUT,timeOUT,nextStep] = dfArduinoGetData(dataOUT,timeOUT,ard,ardCh,trialStartTic,nextStep);
                        case 'serial'
                            [dataOUT,timeOUT,nextStep] = dfArduinoSerial(dataOUT,timeOUT,ard,trialStartTic,nextStep);
                    end

                end
                
                %                 % Set trial start pin to zero after 80 ms
                %                 if currTrigTime > 0.08, writeDigitalPin(ard,ardCh.trialStart,0), end
                
                % Read all data
                switch ArdOrSerial
                    case 'arduino'
                        [dataOUT,timeOUT,nextStep] = dfArduinoGetData(dataOUT,timeOUT,ard,ardCh,trialStartTic,nextStep);
                    case 'serial'
                        [dataOUT,timeOUT,nextStep] = dfArduinoSerial(dataOUT,timeOUT,ard,trialStartTic,nextStep);
                end
                
                % Update current time
                currTrigTime = toc(trialStartTic);
                
                %fprintf('%0.4f - pin%i - trial%i\n',currTrigTime,outVal,ii)
            end
            
            %             % Set back to zero trigger flag for next experiment
            %             alreadyStim = 0;
            
            %             % Set trial start pin to high for 80 ms
            %             switch ArdOrSerial
            %                 case 'arduino',  writeDigitalPin(ard,ardCh.trialON,0);
            %             end
            
            % Read all data
            switch ArdOrSerial
                case 'arduino'
                    [dataOUT,timeOUT,nextStep] = dfArduinoGetData(dataOUT,timeOUT,ard,ardCh,trialStartTic,nextStep);
                case 'serial'
                    [dataOUT,timeOUT,nextStep] = dfArduinoSerial(dataOUT,timeOUT,ard,trialStartTic,nextStep);
            end

            %             writeDigitalPin(ard,ardCh.trialEnd,1)
            %             pause(0.08)
            %             writeDigitalPin(ard,ardCh.trialEnd,0)
        else
            
            
        end
        
        % Plot 'live' pressure values
        LgdStruct = updateLiveResult(timeOUT, dataOUT, expParam, LgdStruct);
        pause(expParam.resPause)
    end
    
    NSD.ardCh   = ardCh;
    NSD.dataOUT = dataOUT;
    NSD.timeOUT = timeOUT;
    
    save(dirs.RecFile, 'NSD')
else
    load(dirs.RecFile)
end

% Close ard object port
if strcmp(ArdOrSerial,'serial'), fclose(ard); end

% title('20 trials, only pressure, pertOut, pertIn, trialON')
% legend('0.05 sec = 50 milliseconds')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions go here!
%%%%%%%%%%%%%%%%%%%%%%%%%

% initLiveResult
function presH = initLiveResult(expParam, defMon)

curSess  = expParam.curSess;
balloon  = expParam.balloon;
balloon(strfind(balloon, '_')) = '';

monitorSize = get(0, 'Monitor');
numMon = size(monitorSize, 1);
plotDim = [800 600];

if numMon == 2 && defMon == 2
    [~, mon] = max(monitorSize(:,1));
    
    halfW  = monitorSize(mon, 3)/2;
    halfWD = halfW - plotDim(1)/2 + monitorSize(mon, 1) - 1;
    
    figPosition = [halfWD 80 plotDim];
else
    
    halfW = monitorSize(1, 3)/2;
    halfWD = halfW - plotDim(1)/2 + monitorSize(1, 1) - 1;
    
    figPosition = [halfWD 80 plotDim];
end
winPos = figPosition;

presH = figure('NumberTitle', 'off', 'Color', [1 1 1], 'Position', winPos);

mark = plot([1 1], [-10 50], 'k-', 'LineWidth', 2);
axis([0 4 -0.5 6.0])
box off
set(gca,'FontSize', 12,...
        'FontWeight', 'bold')
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold') 
ylabel('Voltage (V)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k') 
title({'Pressure Recording, Live Result';
       curSess;
       ['Balloon: ' balloon]})

hold on
drawnow
end

% updateLiveResult
function LgdStruct = updateLiveResult(timeOUT, dataOUT, expParam, LgdStruct)

time     = timeOUT.press;
sig      = dataOUT.press;
numTrial = expParam.numTrial;
curTrial = expParam.curTrial;
trialColors = distinguishable_colors(numTrial);

tag = ['Trial ' num2str(curTrial)];
trPrs = plot(time, sig, 'LineWidth', 2, 'Color', trialColors(curTrial, :));

if curTrial == 1
    LgdStruct.tag = {tag};
    LgdStruct.curve = trPrs;
else
    LgdStruct.tag   = cat(1, LgdStruct.tag, tag);
    LgdStruct.curve = cat(1, LgdStruct.curve, trPrs);
end

lgd = legend(LgdStruct.curve, LgdStruct.tag);
set(lgd, 'box', 'off',...
         'location', 'NorthWest'); 
end
