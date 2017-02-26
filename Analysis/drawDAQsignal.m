function drawDAQsignal(sRate, spans, DAQin, curRecording, saveResultsDir)
spans = spans*8000/16000;

[r, c] = size(spans);
pts = length(DAQin);
time = 0:1/sRate:(pts-1)/sRate;

plotpos = [500 500];
plotdim = [1000 400];

sv2File = 0;

for ii = 1:r
    ForceSensorV(ii) = figure('Color', [1 1 1]);
    set(ForceSensorV(ii), 'Position',[plotpos plotdim],'PaperPositionMode','auto')
    
    perturb = zeros(1, pts);
    perturb(spans(ii,1):spans(ii,2)) = -0.5;
    
    subplot(1,2,1)
    plot(time, perturb, 'k')
    hold on
    plot(time, DAQin(:,1,ii), 'b')
    
    xlabel('Time (s)', 'FontSize', 10, 'FontWeight', 'bold') 
    ylabel('Voltage (V)', 'FontSize', 10, 'FontWeight', 'bold')
    title('Collar Sensor', 'FontSize', 10, 'FontWeight', 'bold')
    axis([0 4 -5 5]); box off
    
    subplot(1,2,2)
    plot(time, perturb, 'k')
    hold on
    plot(time, DAQin(:,2,ii), 'b')
    
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
   
    if sv2File == 1
        plTitle = [curRecording  '_ForceSensor_Test ' num2str(ii)];     
        saveFileName = [saveResultsDir plTitle '.png'];
        export_fig(saveFileName) 
    end
    pause()
    close all
end
end