function drawDAQAlignedPressure(res, saveResultsDir, sv2F, intFlag)
%Plots multiple trials on top of each other. Currently only plotting one 
%sensor. Assumes the trials have been aligned.

curSess  = res.curSess;         % The current experiment details (Subject/Run)
numTrial = res.numPertTrialsNi; % Number of Catch Trials (Only relevant ones)
AudFB    = res.AudFB;
balloon  = res.balloon;

if intFlag == 1
    presSD = res.fSNSD;
    suffix = '_Int';
else
    presSD = res.presSD;
    suffix = '';
end

time     = presSD.timeAl;
sensor   = presSD.sensorAl;
limits   = res.limitsPAl;

LagTimeM  = presSD.lagTimeM(1); % Onset
LagTimeSE = presSD.lagTimeSE(1); % Onset
LagNote = ['Mean Onset Lag: ' num2str(LagTimeM) ' ms ± ' num2str(LagTimeSE) ' ms'];

RiseTimeM  = presSD.riseTimeM(1); % Onset
RiseTimeSE = presSD.riseTimeSE(1); % Onset
RiseNote = ['Mean Rise Time: ' num2str(RiseTimeM) ' ms ± ' num2str(RiseTimeSE) ' ms'];

OnOfValM  = presSD.OnOffValM;
OnOfValSE = presSD.OnOffValSE;
PlatStNote = ['Mean Plataeu (Start) Val: ' num2str(OnOfValM(1)) ' psi ± ' num2str(OnOfValSE(1)) ' psi'];

PresLossM  = presSD.pTrialLossM;
PresLossSE = presSD.pTrialLossSE;
PresLossNote = ['Mean Pressure Loss: ' num2str(PresLossM) ' psi ± ' num2str(PresLossSE) ' psi'];

curSess(strfind(curSess, '_')) = ' ';
balloon(strfind(balloon, '_')) = '';

plotpos = [500 300];
plotdim = [850 600];
trialColors = distinguishable_colors(numTrial);

CombinedSensor = figure('Color', [1 1 1]);
set(CombinedSensor, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

plot([0 0], [-1 20], 'k-', 'LineWidth', 2)

for ii = 1:numTrial
    hold on
    h(ii) = plot(time, sensor(:,ii), 'LineWidth', 2, 'Color', trialColors(ii,:));
end

xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold') 
ylabel('Pressure (psi)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k') 
title({'Mean Pressure Sensor Measurements aligned at Perturbation Onset';
        curSess}, 'FontSize', 12, 'FontWeight', 'bold')
axis(limits);
box off

set(gca,'FontSize', 12,...
        'FontWeight', 'bold')

lgdNames = cell(numTrial, 1);
for i = 1:numTrial; lgdNames{i} = ['Trial ' num2str(i)]; end

pltlgd = legend(h, lgdNames);
set(pltlgd, 'box', 'off',...
            'location', 'NorthWest'); 

t = annotation('textbox',[0.64 0.83 0.9 0.1],...
               'string', {LagNote;...
                          RiseNote;...
                          PlatStNote;...
                          PresLossNote;...
                          ['Balloon: ' balloon ]},...
                'LineStyle','none',...
                'FontWeight','bold',...
                'FontSize',9,...
                'FontName','Arial');

if sv2F == 1
    plTitle = [curSess  '_AlignedPressureRecordings' suffix '.jpg'];     
    saveFileName = fullfile(saveResultsDir, plTitle);
    export_fig(saveFileName) 
end
end