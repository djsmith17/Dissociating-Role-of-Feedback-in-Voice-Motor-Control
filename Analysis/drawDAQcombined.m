function drawDAQcombined(niAn, time, sensor, saveResultsDir, sv2F)
%Plots multiple trials on top of each other. Currently only plotting one 
%sensor. Assumes the trials have been aligned.

curExp   = niAn.curSess;  % The current experiment detials (Subject/Run)
numTrial = niAn.ncTrials; % Number of Catch Trials (Only relevant ones)

meanLagTimeP      = niAn.meanLagTimeP;
meanRiseTimeP     = niAn.meanRiseTimeP;
meanRangePressure = niAn.meanRangePressure;
SensLim_Al        = niAn.pLimits_Al;

curExp(strfind(curExp, '_')) = ' ';

plotpos = [500 300];
plotdim = [800 600];
trialColors = distinguishable_colors(numTrial);

CombinedSensor = figure('Color', [1 1 1]);
set(CombinedSensor, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

plot([1 1], [-1 5], 'k-', 'LineWidth', 2)

for ii = 1:numTrial
    hold on
    h(ii) = plot(time, sensor(:,ii), 'LineWidth', 2, 'Color', trialColors(ii,:));
end

xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold') 
ylabel('Pressure (psi)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k') 
title({'Mean Pressure Sensor Measurements aligned at Perturbation Onset';
        curExp}, 'FontSize', 12, 'FontWeight', 'bold')
axis(SensLim_Al);
box off

set(gca,'FontSize', 12,...
        'XTickLabel', {'-1.0' '-0.5' '0' '0.5' '1.0' '1.5' '2.0' '2.5'},...
        'FontWeight', 'bold')

lgdNames = cell(numTrial, 1);
for i = 1:numTrial; lgdNames{i} = ['Trial ' num2str(i)]; end;

pltlgd = legend(h, lgdNames);
set(pltlgd, 'box', 'off',...
            'location', 'NorthWest'); 

t = annotation('textbox',[0.70 0.7 0.45 0.1],...
               'string', {['Onset Lag: ' num2str(1000*meanLagTimeP(1)) 'ms'];...
                          ['Rise Time: ' num2str(1000*meanRiseTimeP) 'ms'];...
                          ['Onset/Offset Val: ' num2str(meanRangePressure(1), '%1.2f') 'psi, ' num2str(meanRangePressure(2), '%1.2f') 'psi']},...
                'LineStyle','none',...
                'FontWeight','bold',...
                'FontSize',10,...
                'FontName','Arial');

if sv2F == 1
    plTitle = [curExp  '_CombinedDAQSignalOutput.jpg'];     
    saveFileName = fullfile(saveResultsDir, plTitle);
    export_fig(saveFileName) 
end
end