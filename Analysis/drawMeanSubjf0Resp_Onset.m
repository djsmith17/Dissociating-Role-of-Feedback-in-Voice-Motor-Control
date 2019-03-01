function drawMeanSubjf0Resp_Onset(poolRes, targPixDim, plotFolder, fLabel, fStat, varargin)

if isempty(varargin)
    presFlag = 0;
else
    presFlag = varargin{1};
end

curSess          = poolRes.curSess;
cond             = poolRes.cond;
numCond          = length(cond);
numControl       = poolRes.numContTrialsFin;
numPerturb       = poolRes.numPertTrialsFin;

time             = poolRes.secTime;
contf0           = poolRes.audioMf0MeanCont;  
pertf0           = poolRes.audioMf0MeanPert;

statLib = poolRes.statLib;      % Stats Library
statSMM = round(statLib(1), 1); % Mean Stimulus Magnitude (Masked)
statSMV = round(statLib(2), 1); % Mean Stimulus Magnitude (Voice)
statRMM = round(statLib(3), 1); % Mean Response Magnitude (Masked)
statRMV = round(statLib(4), 1); % Mean Response Magnitude (Voice)
statRPM = round(statLib(5));    % Mean Response Percentage (Masked)
statRPV = round(statLib(6));    % Mean Response Percentage (Voice)
statSP  = statLib(7);           % p-value of Stimulus Magnitude t-test
statRP  = statLib(8);           % p-value of Response Magnitude t-test
statPP  = statLib(9);           % p-value of Response Percentage t-test

limits  = poolRes.limitsAmean;
pltName = poolRes.pltName;

timeP       = poolRes.secTimeP;
sensorPAdj  = poolRes.sensorPAdjust;
InflDeflT   = poolRes.InflDeflT;

legLines = [];
legNames = {};

pValueThresh = 0.05;

% Plotting Variables
plotpos        = [10 10];
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

%Onset of Perturbation
if presFlag == 1
    pertAx  = [InflDeflT(1), InflDeflT(2)];
    pertAy  = [600 600];
    
    pA = area(pertAx, pertAy, -600, 'FaceColor', pertColor, 'EdgeColor', pertColor);
    hold on 
    
    plot(timeP, sensorPAdj(:,1), '--m', 'LineWidth', 3)
    hold on
end  

plot(dottedStartx, dottedy,'color',[0.3 0.3 0.3],'LineWidth',lineThick)
hold on

nC = shadedErrorBar(time, contf0(:,1), contf0(:,2), 'lineprops', 'k', 'transparent', 1);
set(nC.mainLine, 'LineWidth', lineThick)
legLines = cat(2, legLines, nC.mainLine);
legNames = cat(2, legNames, 'Control Trials');
hold on

for ii = 1:numCond
    nM = shadedErrorBar(time, pertf0{ii}(:,1), pertf0{ii}(:,2), 'lineprops', condColors(ii), 'transparent', 1);
    set(nM.mainLine, 'LineWidth', lineThick)
    legLines = cat(2, legLines, nM.mainLine);
    legNames = cat(2, legNames, {[cond{ii} ' Trials']});
    hold on
end

xlabel('Time (s)',   'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold'); 
ylabel('f0 (cents)', 'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold')
title({'Mean Participant Response'; 'Onset of Perturbation'}, 'FontName', fontN, 'FontSize', titleFSize, 'FontWeight', 'bold')
axis(limits); box off

set(gca,'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold')

% Done plotting, now to add some annotations
annoStim = ['SM (M/NM): ' num2str(statSMM) ' cents / ' num2str(statSMV) ' cents'];
annoResp = ['RM (M/NM): ' num2str(statRMM) ' cents / ' num2str(statRMV) ' cents'];
annoPerc = ['RP (M/NM): ' num2str(statRPM) '% / ' num2str(statRPV) '%'];

annoStim = checkSig(statSP, pValueThresh, annoStim);
annoResp = checkSig(statRP, pValueThresh, annoResp);
annoPerc = checkSig(statPP, pValueThresh, annoPerc);

if fStat == 1
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
       'Position', [0.73 0.72 0.1 0.1],...
       'box','off',...
       'FontName', fontN,...
       'FontSize', legAnnoFSize,...
       'FontWeight', 'bold');

plots = {'Figure'};
for i = 1:length(plots)
    plTitle = [pltName '_Onset.jpg'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName, '-nocrop')
end
end

function anno = checkSig(stat, thresh, anno)

if stat < thresh
    anno = [anno '*'];
end
end