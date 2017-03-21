function drawDAQsignal(time, DAQin, trigs, pLimits, fLimits, curRecording, saveResultsDir)
%Good for seeing the whole signal
[r, ~] = size(trigs);

plotpos = [500 300];
plotdim = [800 600];
pertColor = [0.8 0.8 0.8];

sv2F = 1;
for ii = 1:r
    ForceSensorV(ii) = figure('Color', [1 1 1]);
    set(ForceSensorV(ii), 'Position',[plotpos plotdim],'PaperPositionMode','auto')
    
    pertAx  = [trigs(ii,1), trigs(ii,2)];
    pertAy  = [6 6];
    
    pA = area(pertAx, pertAy, -1, 'FaceColor', pertColor, 'EdgeColor', pertColor);
    hold on
    %Pressure ('4') followed by Force Sensor Collar ('2') and Force Sensor
    %Neck ('3')
    [aX, fS, pS] = plotyy([time' time'], [DAQin(:,2,ii) DAQin(:,3,ii)], time, DAQin(:,4,ii));
    
    set(pS, 'LineWidth', 2,... 
            'Color', 'k');
    
    set(fS(1), 'LineWidth', 2,... 
               'Color', 'b');      
    set(fS(2), 'LineWidth', 2,... 
               'Color', 'r');
    
    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold') 
    ylabel(aX(1), 'Voltage (V)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k')
    ylabel(aX(2), 'Pressure (psi)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k') 
    title('Sensors Measurements due to Balloon Inflation', 'FontSize', 18, 'FontWeight', 'bold')
    axis(aX(1), fLimits); 
    axis(aX(2), pLimits);
    box off
    
    set(aX, {'ycolor'}, {'k';'k'},...
            'FontSize', 12,...
            'FontWeight', 'bold')
        
    pltlgd = legend([pA pS fS(1) fS(2)], 'Perturbation Period', 'Pressure Sensor', 'Force Sensor: Collar', 'Force Sensor: Neck');
    set(pltlgd, 'box', 'off',...
                'location', 'NorthWest'); 
   
    if sv2F == 1
        plTitle = [curRecording  '_ForceSensorTest ' num2str(ii) '.png'];     
        saveFileName = fullfile(saveResultsDir, plTitle);
        export_fig(saveFileName) 
    end
end
end