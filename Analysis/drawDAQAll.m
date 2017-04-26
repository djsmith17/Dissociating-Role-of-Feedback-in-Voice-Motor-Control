function drawDAQAll(niAn, plotFlag, saveResultsDir, sv2F)
%This will draw every channel that was collected by the DAQ and then
%analyzed by dfAnalysisNIDAQ. This just gives a top down view of what
%happened in this recording. If multiple trials are given, this will either
%plot each trial, or aligned trials, based on flag

curExp   = niAn.curExp;   %The current experiment detials (Subject/Run)
numTrial = niAn.numTrial; %Number of Trials
trigs    = niAn.pertTrig; %Where the perturbations occur

time     = niAn.time;
pertSig  = niAn.pertSig;
sensorFC = niAn.sensorFC;
sensorFN = niAn.sensorFN;
sensorP  = niAn.sensorP;
audioM   = niAn.audioM;
audioH   = niAn.audioH;
sensorO  = niAn.sensorO;

% sensorP        = niAn.sensorP_DN;
% lagTimeP       = niAn.lagsPres;
% riseTimeP      = niAn.riseTimeP;
% rangePressures = niAn.rangePressures;
% pLimits        = niAn.pLimits;
% 
% sensorFC   = niAn.sensorFC_DN;
% sensorFN   = niAn.sensorFN_DN;
% lagTimeFC  = niAn.lagsFC;
% lagTimeFN  = niAn.lagsFN;
% fLimits    = niAn.fLimits;
% 
% time_audio = niAn.time_audio;
% micf0     = niAn.audioMf0_norm;
% headf0    = niAn.audioHf0_norm;
% aLimits   = niAn.aLimits;

curExp(strfind(curExp, '_')) = ' ';

plotpos = [150 150];
plotdim = [1200 600];
pertColor = [0.8 0.8 0.8];

for ii = 1:numTrial
    DAQOutput(ii) = figure('Color', [1 1 1]);
    set(DAQOutput(ii), 'Position',[plotpos plotdim],'PaperPositionMode','auto')
    
    pertAx  = [trigs(ii,1), trigs(ii,2)];
    pertAy  = [200 200];
    
    ha = tight_subplot(2,4,[0.15 0.05],[0.12 0.15],[0.05 0.01]);

    axes(ha(1))
    pA = area(pertAx, pertAy, -200, 'FaceColor', pertColor, 'EdgeColor', pertColor);
    hold on
    plot(time, pertSig(:,ii))
    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('Voltage (V)', 'FontSize', 18, 'FontWeight', 'bold')
    title('Perturbation Signal', 'FontSize', 18, 'FontWeight', 'bold')
    axis([0 4 0 3.5]); box off
    
    axes(ha(2))
    pA = area(pertAx, pertAy, -200, 'FaceColor', pertColor, 'EdgeColor', pertColor);
    hold on
    plot(time, sensorFC(:,ii))
    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('Voltage (V)', 'FontSize', 18, 'FontWeight', 'bold')
    title('Force Sensor: Collar', 'FontSize', 18, 'FontWeight', 'bold')
    axis([0 4 0 5]); box off
    
    axes(ha(3))
    pA = area(pertAx, pertAy, -200, 'FaceColor', pertColor, 'EdgeColor', pertColor);
    hold on
    plot(time, sensorFN(:,ii))
    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('Voltage (V)', 'FontSize', 18, 'FontWeight', 'bold')
    title('Force Sensor: Neck', 'FontSize', 18, 'FontWeight', 'bold')
    axis([0 4 0 5]); box off
    
    axes(ha(4))
    pA = area(pertAx, pertAy, -200, 'FaceColor', pertColor, 'EdgeColor', pertColor);
    hold on
    plot(time, sensorP(:,ii))
    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('Pressure (psi)', 'FontSize', 18, 'FontWeight', 'bold')
    title('Pressure Sensor', 'FontSize', 18, 'FontWeight', 'bold')
    axis([0 4 0 4]); box off
    
    axes(ha(5))
    pA = area(pertAx, pertAy, -200, 'FaceColor', pertColor, 'EdgeColor', pertColor);
    hold on
    plot(time, audioM(:,ii))
    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('SPL (dB)', 'FontSize', 18, 'FontWeight', 'bold')
    title('Microphone', 'FontSize', 18, 'FontWeight', 'bold')
    axis([0 4 -0.01 0.01]); box off
    
    axes(ha(6))
    pA = area(pertAx, pertAy, -200, 'FaceColor', pertColor, 'EdgeColor', pertColor);
    hold on
    plot(time, audioH(:,ii))
    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('SPL (dB)', 'FontSize', 18, 'FontWeight', 'bold')
    title('Headphones', 'FontSize', 18, 'FontWeight', 'bold')
    axis([0 4 -0.01 0.01]); box off
    
    axes(ha(7))
    pA = area(pertAx, pertAy, -200, 'FaceColor', pertColor, 'EdgeColor', pertColor);
    hold on
    plot(time, sensorO(:,ii))
    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('-', 'FontSize', 18, 'FontWeight', 'bold')
    title('Optical TriggerBox', 'FontSize', 18, 'FontWeight', 'bold')
    axis([0 4 0 2]); box off
    
    suptitle([curExp ': All NIDAQ Channels'])    
    
    plotTitle = 'DAQAll';
    if sv2F == 1
        plTitle = [curExp  plotTitle num2str(ii) '.jpg'];     
        saveFileName = fullfile(saveResultsDir, plTitle);
        export_fig(saveFileName) 
    end
end
end