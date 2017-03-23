function drawDAQcombined(time, pSensor, trigs, niAn, pLimits, curExp, saveResultsDir, sv2F)
%Good for seeing the whole signal
[r, ~] = size(trigs);

curExp(strfind(curExp, '_')) = ' ';

plotpos = [500 300];
plotdim = [800 600];

CombinedSensor = figure('Color', [1 1 1]);
set(CombinedSensor, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

colors = distinguishable_colors(r);

plot([1 1], [-1 5], 'k-', 'LineWidth', 2)

for ii = 1:niAn.numTrial
    hold on
    h(ii) = plot(time, pSensor(:,ii), 'LineWidth', 2, 'Color', colors(ii,:));
end

xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold') 
ylabel('Pressure (psi)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k') 
title({'Pressure Sensor Measurements due to Balloon Inflation';
        curExp}, 'FontSize', 12, 'FontWeight', 'bold')
axis(pLimits);
box off

set(gca,'FontSize', 12,...
        'XTickLabel', {'-1.0' '-0.5' '0' '0.5' '1.0' '1.5' '2.0' '2.5'},...
        'FontWeight', 'bold')

pltlgd = legend(h, 'Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'Trial 6', 'Trial 7', 'Trial 8', 'Trial 9', 'Trial 10');
set(pltlgd, 'box', 'off',...
            'location', 'NorthWest'); 

t = annotation('textbox',[0.5 0.82 0.45 0.1],...
               'string', {['Mean Pressure Onset Lag: ' num2str(1000*niAn.PresLagVals(1)) 'ms'];...
                          ['Mean Onset & Offset Pressure Values: ' num2str(niAn.meanRangePressure(1)) 'psi, ' num2str(niAn.meanRangePressure(2)) 'psi']},...
                'LineStyle','none',...
                'FontWeight','bold',...
                'FontSize',10,...
                'FontName','Arial');

if sv2F == 1
    plTitle = [curExp  '_CombinedDAQSignalOutput.png'];     
    saveFileName = fullfile(saveResultsDir, plTitle);
    export_fig(saveFileName) 
end
end