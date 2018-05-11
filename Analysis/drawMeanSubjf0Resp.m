function drawMeanSubjf0Resp(poolRes, targPixDim, plotFolder, fLabel, fStat, varargin)

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

statLib = poolRes.statLib;
statSMM = round(statLib(1), 1);
statSMV = round(statLib(2), 1);
statRMM = round(statLib(3), 1);
statRMV = round(statLib(4), 1);
statRPM = round(statLib(5));
statRPV = round(statLib(6));
statSP  = statLib(7);
statRP  = statLib(8);
statPP  = statLib(9);

limits  = poolRes.limitsAmean;
pltName = poolRes.pltName;

sensorP     = poolRes.sensorPMean;
pressureLim = poolRes.limitsPmean;

legLines = [];
legNames = {};

pValueThresh = 0.05;

% Plotting Variables
plotpos        = [10 100];
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

ha = tight_subplot(1,2,[0.1 0.03],[0.12 0.15],[0.05 0.05]);

%Onset of Perturbation
axes(ha(1))
plot(dottedStartx, dottedy,'color',[0.3 0.3 0.3],'LineWidth',lineThick)
hold on

if presFlag == 1
    yyaxis right
    plot(time, sensorP(:,1), '--k', 'LineWidth', 1.5)

    ylabel('Pressure (psi)')
    axis(pressureLim);
    set(gca,'FontSize', 14,...
            'FontWeight','bold')
    yyaxis left
end  

nC = shadedErrorBar(time, contf0(:,1), contf0(:,2), 'lineprops', 'k', 'transparent', 1);
set(nC.mainLine, 'LineWidth', lineThick)
hold on

for ii = 1:numCond
    nM = shadedErrorBar(time, pertf0{ii}(:,1), pertf0{ii}(:,2), 'lineprops', condColors(ii), 'transparent', 1);
    set(nM.mainLine, 'LineWidth', lineThick)
    hold on
end

xlabel('Time (s)',   'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold'); 
ylabel('f0 (cents)', 'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold')
title('Onset of Perturbation', 'FontName', fontN, 'FontSize', titleFSize, 'FontWeight', 'bold')
axis(limits); box off

set(gca,'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold')

%Offset of Perturbation
axes(ha(2))
plot(dottedStartx, dottedy,'color',[0.3 0.3 0.3],'LineWidth',lineThick)
hold on

fC = shadedErrorBar(time, contf0(:,3), contf0(:,4), 'lineprops', 'k', 'transparent', 1);
set(fC.mainLine, 'LineWidth', lineThick)
legLines = cat(2, legLines, fC.mainLine);
legNames = cat(2, legNames, {[num2str(numControl) ' Control Trials']});
hold on

for ii = 1:numCond
    fM = shadedErrorBar(time, pertf0{ii}(:,3), pertf0{ii}(:,4), 'lineprops', condColors(ii), 'transparent', 1);
    set(fM.mainLine, 'LineWidth', lineThick)
    legLines = cat(2, legLines, fM.mainLine);
    legNames = cat(2, legNames, {[num2str(numPerturb(ii)) ' ' cond{ii} ' Trials']});
    hold on
end

xlabel('Time (s)',   'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold'); 
ylabel('f0 (cents)', 'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold')
title('Offset of Perturbation', 'FontName', fontN, 'FontSize', titleFSize, 'FontWeight', 'bold')
axis(limits); box off

set(gca,'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold',...
        'YAxisLocation', 'right');

sup = suptitle(curSess);
set(sup, 'FontName', fontN,...
         'FontSize', titleFSize,...
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
       'Position', [0.83 0.75 0.1 0.1],...
       'Box', 'off',...
       'Edgecolor', [1 1 1],...
       'FontName', fontN,...
       'FontSize', legAnnoFSize,...
       'FontWeight', 'bold');

plots = {'Figure'};
for i = 1:length(plots)
    plTitle = [pltName '.jpg'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end
end

function anno = checkSig(stat, thresh, anno)

if stat < thresh
    anno = [anno '*'];
end
end