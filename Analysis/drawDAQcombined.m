function drawDAQcombined(time, pSensor, trigs, niAn, pLimits, curRecording, saveResultsDir, sv2F)
%Good for seeing the whole signal
[r, ~] = size(trigs);

plotpos = [500 300];
plotdim = [800 600];
pertColor = [0.8 0.8 0.8];

CombinedSensor = figure('Color', [1 1 1]);
set(CombinedSensor, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

[aX, H] = plotyy([time' time'], [fSensorC(:,ii) fSensorN(:,ii)], time, pSensor(:,ii));

set(H, 'LineWidth', 2,... 
        'Color', 'k');

xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold') 
ylabel(aX, 'Pressure (psi)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k') 
title('Sensors Measurements due to Balloon Inflation', 'FontSize', 18, 'FontWeight', 'bold')
axis(aX, pLimits);
box off

set(aX, 'ycolor', 'k',...
        'FontSize', 12,...
        'FontWeight', 'bold')

pltlgd = legend([pA H fS(1) fS(2)], 'Perturbation Period', 'Pressure Sensor', 'Force Sensor: Collar', 'Force Sensor: Neck');
set(pltlgd, 'box', 'off',...
            'location', 'NorthWest'); 

t = annotation('textbox',[0.55 0.82 0.35 0.1],...
               'string', {['Pressure Sensor Onset Lag: ' num2str(1000*niAn.PresLags(ii,1)) 'ms'];...
                          ['Onset & Offset Pressures: ' num2str(niAn.rangePressures(ii,1)) 'psi, ' num2str(niAn.rangePressures(ii,2)) 'psi']},...
                'LineStyle','none',...
                'FontWeight','bold',...
                'FontSize',10,...
                'FontName','Arial');

if sv2F == 1
    plTitle = [curRecording  '_CombinedDAQSignalOutput.png'];     
    saveFileName = fullfile(saveResultsDir, plTitle);
    export_fig(saveFileName) 
end
end