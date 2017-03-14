function drawDAQsignal(sRate, trigs, DAQin, curRecording, saveResultsDir)
%Good for seeing the whole signal

[r, c] = size(trigs);
pts = length(DAQin);
time = 0:1/sRate:(pts-1)/sRate;

plotpos = [500 500];
plotdim = [800 400];
limits  = [0 4 1 3.5];
pertColor = [0.8 0.8 0.8];

sv2File = 0;

for ii = 1:r
    ForceSensorV(ii) = figure('Color', [1 1 1]);
    set(ForceSensorV(ii), 'Position',[plotpos plotdim],'PaperPositionMode','auto')
    
    pertAx  = [trigs(ii,1),trigs(ii,2)];
    pertAy  = [6 6];
    
    area(pertAx, pertAy, -1, 'FaceColor', pertColor, 'EdgeColor', pertColor)
    hold on
    plot(time, DAQin(:,1,ii), 'b', 'LineWidth', 2) %Collar
    hold on
    plot(time, DAQin(:,2,ii), 'r', 'LineWidth', 2) %Neck
    
    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold') 
    ylabel('Voltage (V)', 'FontSize', 18, 'FontWeight', 'bold')
    title('Voltage Change in Force Sensors due to Balloon Inflation')
    axis(limits); box off
    
    set(gca, 'FontSize', 12,...
             'FontWeight', 'bold')
    
    
    pltlgd = legend('Perturbation Period', 'Collar Sensor', 'NeckSensor');
    set(pltlgd, 'box', 'off',...
                'location', 'NorthWest'); 
   
    if sv2File == 1
        plTitle = [curRecording  '_ForceSensorTest ' num2str(ii) '.png'];     
        saveFileName = fullfile(saveResultsDir, plTitle);
        export_fig(saveFileName) 
    end
end
end