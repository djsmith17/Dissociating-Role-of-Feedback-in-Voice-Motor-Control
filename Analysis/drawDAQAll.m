function drawDAQAll(res, saveResultsDir, sv2F)
% drawDAQAll draws the recordings from all selected channels for a given
% all trials of a recording that was collected by the DAQ and then
%analyzed by dfAnalysisNIDAQ. This just gives a top down view of what
%happened in this recording. If multiple trials are given, this will either
%plot each trial, or aligned trials, based on flag

curExp   = res.curSess;       % Ccurrent experiment details (Subject/Run)
numTrial = res.numPertTrials; % Number of Trials
trigs    = res.pertTrig;      % Where the perturbation occurs

time     = res.time;
pertSig  = res.pertSig;
sensorFC = res.sensorFC;
sensorFN = res.sensorFN;
sensorP  = res.sensorP;
audioM   = res.audioM;
audioH   = res.audioH;

curExp(strfind(curExp, '_')) = ' ';

plotpos = [150 150];
plotdim = [1200 600];
pertColor = [0.8 0.8 0.8];

for ii = 1:numTrial
    DAQOutput(ii) = figure('Color', [1 1 1]);
    set(DAQOutput(ii), 'Position',[plotpos plotdim],'PaperPositionMode','auto')
    
    pertAx  = [trigs(ii,1), trigs(ii,2)];
    pertAy  = [200 200];
    w = res.ctIdx(ii);
    
    ha = tight_subplot(2,3,[0.15 0.08],[0.12 0.15],[0.05 0.01]);

    axes(ha(1))
    pA = area(pertAx, pertAy, -200, 'FaceColor', pertColor, 'EdgeColor', pertColor);
    hold on
    plot(time, pertSig(:,w), 'LineWidth', 3)
    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('Voltage (V)', 'FontSize', 18, 'FontWeight', 'bold')
    title('Perturbation Signal', 'FontSize', 18, 'FontWeight', 'bold')
    axis([0 4 0 3.5]); box off
    set(gca, 'YTick', 0:1:3)
    
    axes(ha(2))
    pA = area(pertAx, pertAy, -200, 'FaceColor', pertColor, 'EdgeColor', pertColor);
    hold on
    plot(time, sensorFC(:,w), 'LineWidth', 3)
    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('Voltage (V)', 'FontSize', 18, 'FontWeight', 'bold')
    title('Force Sensor: Collar', 'FontSize', 18, 'FontWeight', 'bold')
    axis([0 4 0 5]); box off
    
    axes(ha(3))
    pA = area(pertAx, pertAy, -200, 'FaceColor', pertColor, 'EdgeColor', pertColor);
    hold on
    plot(time, sensorFN(:,w), 'LineWidth', 3)
    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('Voltage (V)', 'FontSize', 18, 'FontWeight', 'bold')
    title('Force Sensor: Neck', 'FontSize', 18, 'FontWeight', 'bold')
    axis([0 4 0 5]); box off
    
    axes(ha(4))
    pA = area(pertAx, pertAy, -200, 'FaceColor', pertColor, 'EdgeColor', pertColor);
    hold on
    plot(time, sensorP(:,w), 'LineWidth', 3)
    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('Pressure (psi)', 'FontSize', 18, 'FontWeight', 'bold')
    title('Pressure Sensor', 'FontSize', 18, 'FontWeight', 'bold')
    axis([0 4 0 5]); box off
%     t = annotation('textbox',[0.05 0.32 0.40 0.1],...
%                     'string', {['Onset Lag: ' num2str(1000*lagTimeP(ii,1)) 'ms'];...
%                                ['Rise Time: ' num2str(1000*riseTimeP(ii,1)) 'ms'];...
%                                ['Onset/Offset Values: ' num2str(rangePressures(ii,1), '%1.2f') 'psi, ' num2str(rangePressures(ii,2), '%1.2f') 'psi']},...
%                     'LineStyle','none',...
%                     'FontWeight','bold',...
%                     'FontSize',8,...
%                     'FontName','Arial');
    
    axes(ha(5))
    pA = area(pertAx, pertAy, -200, 'FaceColor', pertColor, 'EdgeColor', pertColor);
    hold on
    plot(time, audioM(:,w), 'LineWidth', 3)
    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('SPL (dB)', 'FontSize', 18, 'FontWeight', 'bold')
    title('Microphone', 'FontSize', 18, 'FontWeight', 'bold')
    axis([0 4 -0.01 0.01]); box off
    
    axes(ha(6))
    pA = area(pertAx, pertAy, -200, 'FaceColor', pertColor, 'EdgeColor', pertColor);
    hold on
    plot(time, audioH(:,w), 'LineWidth', 3)
    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('SPL (dB)', 'FontSize', 18, 'FontWeight', 'bold')
    title('Headphones', 'FontSize', 18, 'FontWeight', 'bold')
    axis([0 4 -0.01 0.01]); box off
    
%     axes(ha(7))
%     pA = area(pertAx, pertAy, -200, 'FaceColor', pertColor, 'EdgeColor', pertColor);
%     hold on
%     plot(time, sensorO(:,ii))
%     xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('-', 'FontSize', 18, 'FontWeight', 'bold')
%     title('Optical TriggerBox', 'FontSize', 18, 'FontWeight', 'bold')
%     axis([0 4 0 2]); box off
    
    suptitle({[curExp ' Trial ' num2str(w)], 'All NIDAQ Channels'})    
    
    plotTitle = 'DAQAllControl';
    if sv2F == 1
        plTitle = [curExp  plotTitle num2str(w) '.jpg'];     
        saveFileName = fullfile(saveResultsDir, plTitle);
        export_fig(saveFileName) 
    end
end
end