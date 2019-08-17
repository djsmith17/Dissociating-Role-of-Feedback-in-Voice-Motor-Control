function drawMeanSubjf0Resp_Onset(poolRes, targPixDim, plotFolder, fLabel, fStat, varargin)

if isempty(varargin)
    presFlag = 0;
else
    presFlag = varargin{1};
end

curSess          = poolRes.curSess;
cond             = poolRes.pubCond;
numCond          = length(cond);
numControl       = poolRes.numContTrialsFin;
numPerturb       = poolRes.numPertTrialsFin;

time             = poolRes.secTime;
contf0           = poolRes.audioMf0MeanCont;  
pertf0M          = poolRes.audioMf0MeanPert;
pertf0H          = poolRes.audioHf0MeanPert;

limits  = poolRes.limitsAmean;
pltName = poolRes.pltName;

timeP       = poolRes.secTimeP;
sensorPAdj  = poolRes.sensorPMean;
InflDeflT   = poolRes.InflDeflT;

legLines = [];
legNames = {};

pValueThresh = 0.05;

% Plotting Variables
plotpos        = [10 40];
plotdim        = targPixDim;
MeanSubjf0Resp = figure('Color', [1 1 1]);
set(MeanSubjf0Resp, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

dottedStartx = [0 0];
dottedy      = [-500 500];
condColors   = {'b', 'r', 'g', 'm'};
fontN        = 'Arial';
legAnnoFSize = 25;
titleFSize   = 35;
axisLSize    = 30;
lineThick    = 4;
pertColor = [0.6 0.6 0.6];
pressureC = [255 87 51]/255; %Auburn

ax2 = axes('Position',[0.1300 0.1100 0.7750 0.8150]);
ax1 = axes('Position',[0.1300 0.1100 0.7750 0.8150]);

if presFlag == 1
    pertAx  = [InflDeflT(1), InflDeflT(2)];
    pertAy  = [600 600];
    
%     pA = area(pertAx, pertAy, -600, 'FaceColor', pertColor, 'EdgeColor', pertColor);
%     hold on 
    axes(ax2)
    pL = plot(timeP, sensorPAdj(:,1), 'color', pressureC, 'LineStyle', ':', 'LineWidth', lineThick);
    
    switch poolRes.expType
        case 'Somatosensory Perturbation_Perceptual'
            measUnits = 'Pressure (psi)';
            measUnitTicks = [0 1 2 3 4];
            measLimits = [-0.5 1.0 0 4.5];
            legendLabel = 'Pressure Trace';
        case 'Auditory Perturbation_Perceptual'
            measUnits = 'Artificial f0 Shift (cents)';
            measUnitTicks = [-100 0 100];
            legendLabel = 'Artificial f0 Shift';
            limits = [-0.5 1.0 -105 105];
            measLimits = [-0.5 1.0 -105 0];
    end    
    ylabel(ax2, measUnits)
    axis(measLimits)    
    
    box off
    set(gca,'FontName', fontN,...
            'FontSize', axisLSize,...
            'FontWeight','bold',...
            'LineWidth', 2,...
            'YTick', measUnitTicks,...
            'YColor', pressureC)
end

axes(ax1)
plot(dottedStartx, dottedy, 'color', [0.3 0.3 0.3],'LineWidth',lineThick)
hold on

nC = shadedErrorBar(time, contf0(:,1), contf0(:,2), 'lineprops', 'k', 'transparent', 1);
set(nC.mainLine, 'LineWidth', lineThick)
legLines = cat(2, legLines, nC.mainLine);
legNames = cat(2, legNames, 'Control Trials');
hold on

for ii = 1:numCond
    nM = shadedErrorBar(time, pertf0M{ii}(:,1), pertf0M{ii}(:,2), 'lineprops', condColors(ii), 'transparent', 1);
    set(nM.mainLine, 'LineWidth', lineThick, 'LineStyle', '--')
    legLines = cat(2, legLines, nM.mainLine);
    legNames = cat(2, legNames, {[cond{ii} ' Trials-Produced']});
    hold on
    
    nM = shadedErrorBar(time, pertf0H{ii}(:,1), pertf0H{ii}(:,2), 'lineprops', condColors(ii), 'transparent', 1);
    set(nM.mainLine, 'LineWidth', lineThick, 'LineStyle', '-')
    legLines = cat(2, legLines, nM.mainLine);
    legNames = cat(2, legNames, {[cond{ii} ' Trials-Heard']});
    hold on
end

xlabel('Time (s)', 'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold'); 
ylabel(ax1, 'f0 (cents)', 'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold')
% title({'Mean Participant Response'; 'Onset of Perturbation'}, 'FontName', fontN, 'FontSize', titleFSize, 'FontWeight', 'bold')
axis(limits); box off

set(gca,'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold',...
        'LineWidth', 2)
    
ax2.YAxisLocation = 'right';
ax2.YDir = 'reverse';
set(gca, 'Color', 'None')
if presFlag == 1
    legLines = cat(2, legLines, pL);
    legNames = cat(2, legNames, legendLabel);
end

if fStat == 1
    % Done plotting, now to add some annotations
    annoStim = ['SM (M/NM): ' num2str(statSMM) ' cents / ' num2str(statSMV) ' cents'];
    annoResp = ['RM (M/NM): ' num2str(statRMM) ' cents / ' num2str(statRMV) ' cents'];
    annoPerc = ['RP (M/NM): ' num2str(statRPM) '% / ' num2str(statRPV) '%'];

    annoStim = checkSig(statSP, pValueThresh, annoStim);
    annoResp = checkSig(statRP, pValueThresh, annoResp);
    annoPerc = checkSig(statPP, pValueThresh, annoPerc);
    
    statBox = annotation('textbox',[.30 .75 0.45 0.1],...
                         'string', {annoStim;
                                    annoResp
                                    annoPerc},...
                          'LineStyle','none',...
                          'FontName', fontN,...
                          'FontSize', legAnnoFSize,...
                          'FontWeight','bold');
end

if fLabel == 1
    figureL = pltName(end);
    figureMark = annotation('textbox', [0.01 0.88 0.05 0.1],...
                            'string', figureL,...
                            'LineStyle', 'none',...
                            'FontName', fontN,...
                            'FontSize', titleFSize,...
                            'FontWeight','bold');
end

legend(legLines, legNames,...
       'Position', [0.72 0.15 0.1 0.1],...
       'FontName', fontN,...
       'FontSize', legAnnoFSize,...
       'FontWeight', 'bold');

plots = {'Figure'};
for i = 1:length(plots)
    plTitle = [pltName '_Onset.jpg'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end
end

function anno = checkSig(stat, thresh, anno)

if stat < thresh
    anno = [anno '*'];
end
end