function RecordForceSensorVoltage(n)
close all; 

pltFolder = 'C:\Users\djsmith\Documents\Pilot Data\Force Sensor Calibration\';
method = 'null';

s = initNIDAQ;

numTrial    = n;
trialLen    = 4; %Seconds
trialLenPts = trialLen*s.Rate;
pauseLen    = 0;

negVolSrc = zeros(s.Rate*trialLen, 1) - 1;
negVolSrc(1) = 0; negVolSrc(end) = 0;

trialType = ones(numTrial,1);

[sigs, spans] = createPerturbSignal(s, numTrial, trialLen, trialType);

svData = [];
for ii = 1:numTrial
    NIDAQsig = [sigs(:,ii) negVolSrc];
    queueOutputData(s, NIDAQsig);
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

function plot_data_DAQ(s, spans, svData, method, pltFolder)

[r, c] = size(spans);
pts = length(svData);
time = 0:1/s.Rate:(pts-1)/s.Rate;

plotpos = [500 500];
plotdim = [1000 400];

for ii = 1:r
    perturb = zeros(1, pts);
    perturb(spans(ii,1):spans(ii,2)) = -0.5;
    
    ForceSensorV(ii) = figure('Color', [1 1 1]);
    set(ForceSensorV(ii), 'Position',[plotpos plotdim],'PaperPositionMode','auto')
    
    subplot(1,2,1)
    plot(time, perturb, 'k')
    hold on
    plot(time, svData(:,1,ii), 'b')
    
    xlabel('Time (s)', 'FontSize', 10, 'FontWeight', 'bold') 
    ylabel('Voltage (V)', 'FontSize', 10, 'FontWeight', 'bold')
    title('Collar Sensor', 'FontSize', 10, 'FontWeight', 'bold')
    axis([0 4 -5 5]); box off
    
    subplot(1,2,2)
    plot(time, perturb, 'k')
    hold on
    plot(time, svData(:,2,ii), 'b')
    
    xlabel('Time (s)', 'FontSize', 10, 'FontWeight', 'bold')
    ylabel('Voltage (V)', 'FontSize', 10, 'FontWeight', 'bold')
    title('Neck Sensor', 'FontSize', 10, 'FontWeight', 'bold')
    axis([0 4 -5 5]); box off
    
    suptitle('Voltage Change in Force Sensors due to Balloon Inflation')
    
    pltlgd = legend('Perturbation', 'Voltage from Force Sensor');
    set(pltlgd, 'box', 'off',...
                'location', 'best');
   
    set(gca, 'FontSize', 12,...
             'FontWeight', 'bold')
   
    plTitle = [method  '_ForceSensor_Test ' num2str(ii)];     
    saveFileName = [pltFolder plTitle '.png'];
    export_fig(saveFileName)   
end
end