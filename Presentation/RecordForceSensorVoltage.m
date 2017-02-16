function RecordForceSensorVoltage(n)
close all; 

pltFolder = 'C:\Users\djsmith\Documents\Pilot Data\Force Sensor Calibration\';
method = 'Pilot1';

s = initPerturb;

numTrial    = n;
trialLen    = 4; %Seconds
trialLenPts = trialLen*s.Rate;
pauseLen    = 0;

trialType = zeros(numTrial,1) + 1; %orderTrialsLaced(numTrial, 0.25); %numTrials, percentCatch

[sigs, spans] = createPerturbs(s, numTrial, trialLenPts, trialType);

svData = [];
for ii = 1:numTrial
    queueOutputData(s, sigs(:,ii));
    fprintf('Running Trial %d\n', ii)
    [data_DAQ, time] = s.startForeground;
    
    svData = cat(3, svData, data_DAQ);
    
    pause(pauseLen)      
end

plot_data_DAQ(s, spans, svData, method, pltFolder)

ForceSensorData.pltFolder   = pltFolder;
ForceSensorData.method      = method;
ForceSensorData.sRate       = s.Rate;
ForceSensorData.numTrial    = numTrial;
ForceSensorData.trialLen    = trialLen;
ForceSensorData.trialLenPts = trialLenPts;
ForceSensorData.pauseLen    = pauseLen;
ForceSensorData.trialType   = trialType;
ForceSensorData.sigs        = sigs;
ForceSensorData.spans       = spans;
ForceSensorData.svData      = svData;

save([pltFolder method '_ForceSensorData.mat'],'ForceSensorData')
end

function s = initPerturb

s = daq.createSession('ni');
addAnalogOutputChannel(s,'Dev3',0,'Voltage'); %Output to the Perturbatron
addAnalogOutputChannel(s,'Dev3',1,'Voltage'); %Output voltage to Force Sensor
addAnalogInputChannel(s,'Dev3',0,'Voltage'); %Input from Force Sensor 1
addAnalogInputChannel(s,'Dev3',1,'Voltage'); %Input from Force Sensor 1

s.Rate = 8000;
end

function trialType = orderTrialsLaced(numTrial, per)
%This function organizes trials in sets based on variable 'per'. 
%For example if per = 0.25, then a set of four trials will have 
%three control trials, and one catch trial. The catch trial will always be 
%the second or third trial in the set. This will generate enough sets for
%all trials in 'numTrial'.

trialType = [];
for ii = 1:(numTrial*per)
    place = round(rand) + 2; %This will break if per!= 0.25 **Fix this**
    set = zeros(1,(1/per));
    set(place) = 1;
    
    trialType = cat(2, trialType, set);
end
end

function [sigs, spans] = createPerturbs(s, numTrial, trialLen, trialType)
%This function creates the digital signal for the NIDAQ needed for each
%trial. It also keeps track of the time points when the perturbation should
%activate and deactivate. This function calculates these points for the 
%sampling rate of the NIDAQ (s.Rate). 'spans' is then converted to the 
%Audapter sampling rate (sRate) in the main function.

%The perturbation starts at a semi-random time between 1.7s and 2.1s after 
%phonation. The pertrubation stops at a semi-random time between 
%0.7 and 1.1s after the start of pertrubation.

%s:         NIDAQ object handle for keeping track of variables and I/O
%numTrial:  The number of trials
%trialLen:  The length of each trial in points
%trialType: Vector (length = numTrial) of order of trials (control/catch)

%sigs:  Per-trial digital signal to be outputted to NIDAQ
%spans: Per-trial pertrubation start and stop points to be aligned with mic data

sigs  = zeros(trialLen, numTrial);
spans = zeros(numTrial,2);
for i = 1:numTrial
    St   = round(s.Rate*(1.7 + (2.1-1.7)*rand)); %Hardset (1.7-2.1 seconds)
    pLen = round(s.Rate*(0.7 + (1.1-0.7)*rand)); %Hardset (0.7-1.1 seconds)   
    Sp   = St+pLen;    
    span = St:Sp; 

    sig  = zeros(trialLen,1);
    if trialType(i) == 1
        sig(span) = 3; %3V
    end
    
    sigs(:,i)  = sig;
    spans(i,:) = [St Sp];
end
end

function plot_data_DAQ(s, spans, svData, method, pltFolder)

[r, c] = size(spans);
pts = length(svData);
time = 0:1/s.Rate:(pts-1)/s.Rate;

for ii = 1:r
    perturb = zeros(1, pts);
    perturb(spans(ii,1):spans(ii,2)) = -0.5;
    
    ForceSensorV(ii) = figure('Color', [1 1 1]);
    
    subplot(1,2,1)
    plot(time, perturb, 'k')
    hold on
    plot(time, svData(:,1,ii), 'b')
    
    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); 
    ylabel('Voltage (V)', 'FontSize', 18, 'FontWeight', 'bold')
    title('Collar Sensor')
    axis([0 4 -5 5]); box off
    
    subplot(1,2,2)
    plot(time, perturb, 'k')
    hold on
    plot(time, svData(:,2,ii), 'b')
    
    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); 
    ylabel('Voltage (V)', 'FontSize', 18, 'FontWeight', 'bold')
    title('Neck Sensor')
    axis([0 4 -5 5]); box off
    
    legend('Perturbation', 'Voltage from Force Sensor', 'location', 'best')
    
    suptitle('Voltage Change in Force Sensors due to Balloon Inflation')
    set(gca, 'FontSize', 10,...
             'FontWeight', 'bold')
   
    plTitle = [method  '_ForceSensor_Test ' num2str(ii)];     
    saveFileName = [pltFolder plTitle '.png'];
    export_fig(saveFileName)   
end
end