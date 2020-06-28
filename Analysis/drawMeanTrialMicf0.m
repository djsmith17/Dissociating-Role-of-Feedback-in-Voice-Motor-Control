function drawMeanTrialMicf0(res, plotFolder, varargin)
% drawDAQMeanTrialMicf0(res, plotFolder) plots differences in microphone 
% recordings between perturbed and control trials. 

if isempty(varargin)
    presFlag = 0;
else
    presFlag = varargin{1};
end

curSess          = res.curSess;
f0b              = round(res.f0b, 1); % Baseline f0 rounded to 0.1 Hz
AudFB            = res.AudFB;
numCT            = res.numContTrialsFin;
numPT            = res.numPertTrialsFin;

time             = res.secTime;
meanf0PertOnset  = res.audioMf0MeanPert(:,1);
CIf0PertOnset    = res.audioMf0MeanPert(:,2);
meanf0PertOffset = res.audioMf0MeanPert(:,3);
CIf0PertOffset   = res.audioMf0MeanPert(:,4);

meanf0ContOnset  = res.audioMf0MeanCont(:,1);
CIf0ContOnset    = res.audioMf0MeanCont(:,2);
meanf0ContOffset = res.audioMf0MeanCont(:,3);
CIf0ContOffset   = res.audioMf0MeanCont(:,4);
limits           = res.limitsAmean;

timeP   = res.presSDsv.timeSec;
sensorP = res.presSDsv.sensorSecM;

plotpos = [10 100];
plotdim = [1600 600];
InterTrialf0 = figure('Color', [1 1 1]);
set(InterTrialf0, 'Position', [plotpos plotdim],'PaperPositionMode','auto')

AD = res.audioDynamics;
if ~isempty(AD)
    statSM = round(AD.respVarM(2), 1);
    statRM = round(AD.respVarM(3), 1);
    statRP = round(AD.respVarM(4));
else
    statSM = '';
    statRM = '';
    statRP = '';
end
curSess(strfind(curSess, '_')) = ' ';

lgdCurv = [];
lgdLabl = {};

dottedStartx   = [0 0];
dottedy        = [-300 300];
controlColor   = 'k';
perturbedColor = 'b';
pertLineC      = [0.3 0.3 0.3];
fontN          = 'Arial';
legAnnoFSize   = 12;
titleFSize     = 14;
axisLSize      = 14;
lineThick      = 4;

ha = tight_subplot(1,2,[0.1 0.03],[0.12 0.15],[0.05 0.05]);

%Onset of Perturbation
axes(ha(1))

plot(dottedStartx, dottedy, 'color', pertLineC, 'LineWidth', lineThick)
hold on
if ~isempty(meanf0ContOnset)
    cH = shadedErrorBar(time, meanf0ContOnset, CIf0ContOnset, 'lineprops', {'color', controlColor}, 'transparent', 1); %Unperturbed
    lgdCurv = [lgdCurv cH.mainLine];
    lgdLabl = [lgdLabl, [num2str(numCT) ' Control Trials']];
    set(cH.mainLine, 'LineWidth', lineThick)
    hold on
end
pH = shadedErrorBar(time, meanf0PertOnset, CIf0PertOnset, 'lineprops', {'color', perturbedColor}, 'transparent', 1); %Perturbed
lgdCurv = [lgdCurv pH.mainLine];
lgdLabl = [lgdLabl, [num2str(numPT) ' Perturbed Trials']];

set(pH.mainLine, 'LineWidth', lineThick)
xlabel('Time (s)', 'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold'); 
ylabel('f0 (cents)', 'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold')
title('Onset of Perturbation', 'FontName', fontN, 'FontSize', titleFSize, 'FontWeight', 'bold')
axis(limits); box off

set(gca,'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold')

if presFlag == 1
%     pertAx  = [InflDeflT(1), InflDeflT(2)];
%     pertAy  = [600 600];
%     
%     pA = area(pertAx, pertAy, -600, 'FaceColor', pertColor, 'EdgeColor', pertColor);
%     hold on 
    yyaxis right
    plot(timeP, sensorP(:,1), '--m', 'LineWidth', 3)
end  

%Offset of Perturbation
axes(ha(2))

plot(dottedStartx, dottedy, 'color', pertLineC, 'LineWidth', lineThick)
hold on
if ~isempty(meanf0ContOffset)
    cH2 = shadedErrorBar(time, meanf0ContOffset, CIf0ContOffset, 'lineprops', {'color', controlColor}, 'transparent', 1); %Unperturbed
    set(cH2.mainLine, 'LineWidth', lineThick)
    hold on
end
pH2 = shadedErrorBar(time, meanf0PertOffset, CIf0PertOffset, 'lineprops', {'color', perturbedColor}, 'transparent', 1); %Perturbed
hold on

set(pH2.mainLine, 'LineWidth', lineThick)
xlabel('Time (s)', 'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold'); 
ylabel('f0 (cents)', 'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold')
title('Offset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
axis(limits); box off

set(gca,'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold',...
        'YAxisLocation', 'right');

if presFlag == 1    
%     pertAx  = [InflDeflT(3), InflDeflT(4)];
%     pertAy  = [600 600];
%     
%     pA = area(pertAx, pertAy, -600, 'FaceColor', pertColor, 'EdgeColor', pertColor);
%     hold on 
    yyaxis right
    plot(timeP, sensorP(:,3), '--m', 'LineWidth', 3)
end 
    
sup = suptitle({curSess; ['AudFB: ' AudFB]; ['f0: ' num2str(f0b) 'Hz']});
set(sup, 'FontName', fontN,...
         'FontSize', titleFSize,...
         'FontWeight','bold')

annoStim = ['SM: ' num2str(statSM) ' cents'];
annoResp = ['RM: ' num2str(statRM) ' cents'];
annoPerc = ['RP: ' num2str(statRP) ' %'];

statBox = annotation('textbox',[.38 .75 0.45 0.1],...
                     'string', {annoStim;
                                annoResp
                                annoPerc},...
                        'LineStyle','none',...
                        'FontWeight','bold',...
                        'FontSize',12,...
                        'FontName','Arial');
                    
legend(lgdCurv, lgdLabl,...
            'Box', 'off',...
            'Edgecolor', [1 1 1],...
            'FontSize', 12,...
            'FontWeight', 'bold',...
            'Position', [0.8 0.93 0.05 0.05]);         
                    
plots = {'InterTrialf0'};
for i = 1:length(plots)
    plTitle = [curSess '_' plots{i} '.jpg'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end
end