function drawMeanTrialMicf0(res, plotFolder)
%drawDAQMeanTrialMicf0(niRes, plotFolder) is plotting function for
%displaying differences in microphone channel recordings between perturbed
%and control trials. If your data set does not have any control trials,
%this will give an error

curSess          = res.curSess;
f0b              = round(10*res.f0b)/10;
f0Type           = res.f0Type;
etMH             = res.etMH;
AudFB            = res.AudFB;
numCT            = res.numContTrialsPP;
numPT            = res.numPertTrialsPP;

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

statSM = round(10*res.respVarM(2))/10;
statRM = round(10*res.respVarM(3))/10;
statRP = round(res.respVarM(4));

lgdCurv = [];
lgdLabl = {};

plotpos = [10 100];
plotdim = [1600 600];
InterTrialf0 = figure('Color', [1 1 1]);
set(InterTrialf0, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

dottedStartx = [0 0];
dottedy      = [-300 300];

ha = tight_subplot(1,2,[0.1 0.05],[0.12 0.15],[0.08 0.08]);

%Onset of Perturbation
axes(ha(1))
if ~isempty(meanf0ContOnset)
    uH = shadedErrorBar(time, meanf0ContOnset, CIf0ContOnset, 'lineprops', 'b', 'transparent', 1); %Unperturbed
    lgdCurv = [lgdCurv uH.mainLine];
    lgdLabl = [lgdLabl, [num2str(numCT) ' Control Trials']];
    hold on
end
pH = shadedErrorBar(time, meanf0PertOnset, CIf0PertOnset, 'lineprops', 'r', 'transparent', 1); %Perturbed
lgdCurv = [lgdCurv pH.mainLine];
lgdLabl = [lgdLabl, [num2str(numPT) ' Perturbed Trials']];
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')

title('Onset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
axis(limits); box off

set(gca,'XTickLabel',{'-0.5' '0' '0.5' '1.0'},...
        'FontSize', 16,...
        'FontWeight','bold')

%Offset of Perturbation
axes(ha(2))
if ~isempty(meanf0ContOffset)
    shadedErrorBar(time, meanf0ContOffset, CIf0ContOffset, 'lineprops', 'b', 'transparent', 1)  %Unperturbed
    hold on
end
shadedErrorBar(time, meanf0PertOffset, CIf0PertOffset, 'lineprops', 'r', 'transparent', 1) %Perturbed
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')

title('Offset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
axis(limits); box off
set(gca,'XTickLabel', {'-0.5' '0' '0.5' '1.0'},...
        'FontSize', 16,...
        'FontWeight','bold',...
        'YAxisLocation', 'right');

suptitle({curSess; ['AudFB: ' AudFB]; ['f0: ' num2str(f0b) 'Hz']})

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
            'Position', [0.35 0.2 0.1 0.1]);
        
timeBox = annotation('textbox',[.80 .88 0.45 0.1],...
                     'string', {f0Type;
                            ['Analysis Time: ' num2str(etMH) ' min']},...
                        'LineStyle','none',...
                        'FontWeight','bold',...
                        'FontSize',8,...
                        'FontName','Arial');

plots = {'InterTrialf0'};
for i = 1:length(plots)
    plTitle = [curSess '_' plots{i} f0Type '.jpg'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end
end