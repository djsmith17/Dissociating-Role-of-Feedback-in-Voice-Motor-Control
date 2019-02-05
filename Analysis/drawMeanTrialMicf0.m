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
% f0Type           = res.f0Type;
% etMH             = res.etMH;
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

if ~isempty(res.respVarM)
    statSM = round(res.respVarM(2), 1);
    statRM = round(res.respVarM(3), 1);
    statRP = round(res.respVarM(4));
else
    statSM = '';
    statRM = '';
    statRP = '';
end
curSess(strfind(curSess, '_')) = ' ';

lgdCurv = [];
lgdLabl = {};

controlColor   = [0 0 0];
perturbedColor = [0 0 1];
lineThick      = 4;

plotpos = [10 100];
plotdim = [1600 600];
InterTrialf0 = figure('Color', [1 1 1]);
set(InterTrialf0, 'Position', [plotpos plotdim],'PaperPositionMode','auto')

dottedStartx = [0 0];
dottedy      = [-300 300];

ha = tight_subplot(1,2,[0.1 0.05],[0.12 0.15],[0.08 0.08]);
%Onset of Perturbation
axes(ha(1))
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
hold on
plot(dottedStartx, dottedy, 'k', 'LineWidth', 4)

set(pH.mainLine, 'LineWidth', lineThick)
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')

title('Onset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
axis(limits); box off

set(gca,'YTickLabelMode', 'auto',...
        'XTickLabel',{'-0.5' '0' '0.5' '1.0'},...
        'FontSize', 16,...
        'FontWeight','bold')

%Offset of Perturbation
axes(ha(2))
if ~isempty(meanf0ContOffset)
    cH2 = shadedErrorBar(time, meanf0ContOffset, CIf0ContOffset, 'lineprops', {'color', controlColor}, 'transparent', 1); %Unperturbed
    set(cH2.mainLine, 'LineWidth', lineThick)
    hold on
end
pH2 = shadedErrorBar(time, meanf0PertOffset, CIf0PertOffset, 'lineprops', {'color', perturbedColor}, 'transparent', 1); %Perturbed
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
set(pH2.mainLine, 'LineWidth', lineThick)
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')

title('Offset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
axis(limits); box off
set(gca,'YTickLabelMode', 'auto',...
        'XTickLabel', {'-0.5' '0' '0.5' '1.0'},...
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
            'Position', [0.8 0.93 0.05 0.05]);
        
% timeBox = annotation('textbox',[.80 .88 0.45 0.1],...
%                      'string', {f0Type;
%                             ['Analysis Time: ' num2str(etMH) ' min']},...
%                         'LineStyle','none',...
%                         'FontWeight','bold',...
%                         'FontSize',8,...
%                         'FontName','Arial');
         
                    
plots = {'InterTrialf0'};
for i = 1:length(plots)
    plTitle = [curSess '_' plots{i} '.jpg'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end
end