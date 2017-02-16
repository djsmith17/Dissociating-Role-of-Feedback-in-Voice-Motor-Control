function RecordForceSensorVoltage(n)
close all; 

pltFolder = 'C:\Users\djsmith\Documents\Pilot Data\Force Sensor Calibration\';
method = 'Pilot1';

s = initNIDAQ;

numTrial    = n;
trialLen    = 4; %Seconds
trialLenPts = trialLen*s.Rate;
pauseLen    = 0;

sig = zeros(s.Rate*trialLen,1) - 1;
sig(1) = 0; sig(end) = 0;

trialType = zeros(numTrial,1) + 1; %orderTrialsLaced(numTrial, 0.25); %numTrials, percentCatch

[sigs, spans] = createPerturbSignal(s, numTrial, trialLenPts, trialType);

svData = [];
for ii = 1:numTrial
    queueOutputData(s, [sigs(:,ii) sig]);
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