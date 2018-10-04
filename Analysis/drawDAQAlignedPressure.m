function drawDAQAlignedPressure(res, saveResultsDir, sv2F)
%Plots multiple trials on top of each other. Currently only plotting one 
%sensor. Assumes the trials have been aligned.

curSess  = res.curSess;         % The current experiment details (Subject/Run)
numTrial = res.numPertTrialsNi; % Number of Catch Trials (Only relevant ones)
AudFB    = res.AudFB;

time     = res.timeSAl;
sensor   = res.sensorPAl;
limits   = res.limitsPAl;

LagTimeM  = 1000*res.lagTimePm(1,1);
LagTimeSE = round(1000*res.lagTimePm(2,1),1);
LagNote = ['Mean Onset Lag: ' num2str(LagTimeM) ' ms ± ' num2str(LagTimeSE) ' ms'];

RiseTimeM  = 1000*res.riseTimePm;
RiseTimeSE = round(1000*res.riseTimePSE,1);
RiseNote = ['Mean Rise Time: ' num2str(RiseTimeM) ' ms ± ' num2str(RiseTimeSE) ' ms'];

PresLossM  = round(res.pTrialLossPm, 3);
PresLossSE = round(res.pTrialLossPSE, 3);
PresLossNote = ['Mean Pressure Loss: ' num2str(PresLossM) ' psi ± ' num2str(PresLossSE) ' psi'];

if isfield(res, 'balloon')
    balloon = res.balloon;
else
    balloon = 'N/A';
end

curSess(strfind(curSess, '_')) = ' ';
balloon(strfind(balloon, '_')) = '';

plotpos = [500 300];
plotdim = [800 600];
trialColors = distinguishable_colors(numTrial);

CombinedSensor = figure('Color', [1 1 1]);
set(CombinedSensor, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

plot([0 0], [-1 5], 'k-', 'LineWidth', 2)

for ii = 1:numTrial
    hold on
    h(ii) = plot(time, sensor(:,ii), 'LineWidth', 2, 'Color', trialColors(ii,:));
end

xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold') 
ylabel('Pressure (psi)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k') 
title({'Mean Pressure Sensor Measurements aligned at Perturbation Onset';
        curSess;
        ['AudFB: ' AudFB]}, 'FontSize', 12, 'FontWeight', 'bold')
axis(limits);
box off

set(gca,'FontSize', 12,...
        'FontWeight', 'bold')

lgdNames = cell(numTrial, 1);
for i = 1:numTrial; lgdNames{i} = ['Trial ' num2str(i)]; end

pltlgd = legend(h, lgdNames);
set(pltlgd, 'box', 'off',...
            'location', 'NorthWest'); 

t = annotation('textbox',[0.65 0.76 0.45 0.1],...
               'string', {LagNote;...
                          RiseNote;...
                          PresLossNote;...
                          ['Balloon: ' balloon ]},...
                'LineStyle','none',...
                'FontWeight','bold',...
                'FontSize',10,...
                'FontName','Arial');

if sv2F == 1
    plTitle = [curSess  '_AlignedPressureRecordings.jpg'];     
    saveFileName = fullfile(saveResultsDir, plTitle);
    export_fig(saveFileName) 
end
end