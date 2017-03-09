function drawDAQsignal(sRate, trigs, DAQin, curRecording, saveResultsDir)
%Good for seeing the whole signal

[r, c] = size(trigs);
pts = length(DAQin);
time = 0:1/sRate:(pts-1)/sRate;

plotpos = [500 500];
plotdim = [800 400];

sv2File = 0;

for ii = 1:r
    ForceSensorV(ii) = figure('Color', [1 1 1]);
    set(ForceSensorV(ii), 'Position',[plotpos plotdim],'PaperPositionMode','auto')
    
    ha = tight_subplot(1,2,[0.1 0.05],[0.12 0.15],[0.05 0.03]);
    
    perturb = zeros(1, pts);
    perturb(trigs(ii,1):trigs(ii,2)) = -0.5;
    
    axes(ha(1))
    plot(time, perturb, 'k')
    hold on
    plot(time, DAQin(:,1,ii), 'b')
    
    xlabel('Time (s)', 'FontSize', 10, 'FontWeight', 'bold') 
    ylabel('Voltage (V)', 'FontSize', 10, 'FontWeight', 'bold')
    title('Collar Sensor', 'FontSize', 10, 'FontWeight', 'bold')
    axis([0 4 -5 5]); box off
    
    set(gca, 'FontSize', 12,...
             'FontWeight', 'bold',...
             'YTick', -5:1:5)
    
    axes(ha(2))
    plot(time, perturb, 'k')
    hold on
    plot(time, DAQin(:,2,ii), 'b')
    
    xlabel('Time (s)', 'FontSize', 10, 'FontWeight', 'bold')
    ylabel('Voltage (V)', 'FontSize', 10, 'FontWeight', 'bold')
    title('Neck Sensor', 'FontSize', 10, 'FontWeight', 'bold')
    axis([0 4 -5 5]); box off
    
    set(gca, 'FontSize', 12,...
             'FontWeight', 'bold',...
             'YTick', -5:1:5)
    
    suptitle('Voltage Change in Force Sensors due to Balloon Inflation')
    
    pltlgd = legend('Perturbation', 'Voltage from Force Sensor');
    set(pltlgd, 'box', 'off',...
                'location', 'best'); 
   
    if sv2File == 1
        plTitle = [curRecording  '_ForceSensorTest ' num2str(ii) '.png'];     
        saveFileName = fullfile(saveResultsDir, plTitle);
        export_fig(saveFileName) 
    end
end
end