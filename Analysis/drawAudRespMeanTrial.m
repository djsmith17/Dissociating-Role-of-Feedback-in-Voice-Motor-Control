function drawAudRespMeanTrial(res, plotFolder)
% drawAudRespMeanTrial(res, plotFolder) displays results of pitch-shift
% experiments 

curSess          = res.curSess;
f0b              = round(res.f0b, 1); % Baseline f0 rounded to 0.1 Hz
AudFB            = res.AudFB;
numPT            = res.numPertTrialsFin;

time             = res.secTime;
meanf0MicOnset   = res.audioMf0MeanPert(:,1);
CIf0MicOnset     = res.audioMf0MeanPert(:,2);
meanf0MicOffset  = res.audioMf0MeanPert(:,3);
CIf0MicOffset    = res.audioMf0MeanPert(:,4);

meanf0HeadOnset  = res.audioHf0MeanPert(:,1);
CIf0HeadOnset    = res.audioHf0MeanPert(:,2);
meanf0HeadOffset = res.audioHf0MeanPert(:,3);
CIf0HeadOffset   = res.audioHf0MeanPert(:,4);
limits           = res.limitsAMH;

plotpos = [10 100];
plotdim = [1600 600];
MeanTrialAudResp = figure('Color', [1 1 1]);
set(MeanTrialAudResp, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

dottedStartx = [0 0];
dottedy      = [-500 500];
micColor     = 'b';
headColor    = 'r';
pertLineC    = [0.3 0.3 0.3];
fontN        = 'Arial';
legAnnoFSize = 10;
titleFSize   = 14;
axisLSize    = 12;
lineThick    = 1;

ha = tight_subplot(1,2,[0.1 0.03],[0.12 0.15],[0.05 0.05]);

%Onset of Perturbation
axes(ha(1))
mH = shadedErrorBar(time, meanf0MicOnset, CIf0MicOnset, 'lineprops', micColor, 'transparent', 1); %Pertrubed Microphone
hold on
hH = shadedErrorBar(time, meanf0HeadOnset, CIf0HeadOnset, 'lineprops', headColor, 'transparent', 1); %Perturbed Headphones
hold on
plot(dottedStartx, dottedy, 'color', pertLineC, 'LineWidth', 4)

set(mH.mainLine, 'LineWidth', lineThick)
set(hH.mainLine, 'LineWidth', lineThick)
xlabel('Time (s)', 'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold')
ylabel('f0 (cents)', 'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold')
title('Onset of Perturbation', 'FontName', fontN, 'FontSize', titleFSize, 'FontWeight', 'bold')
axis(limits); box off

set(gca,'XTickLabel',{'-0.5' '0' '0.5' '1.0'},...
        'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold')

%Offset of Perturbation
axes(ha(2))
plot(dottedStartx, dottedy, 'color', pertLineC, 'LineWidth', lineThick)
hold on
mH2 = shadedErrorBar(time, meanf0MicOffset, CIf0MicOffset, 'lineprops', micColor, 'transparent', 1);  %Perturbed Microphone
hold on
hH2 = shadedErrorBar(time, meanf0HeadOffset, CIf0HeadOffset, 'lineprops', headColor, 'transparent', 1); %Perturbed Headphones

set(mH2.mainLine, 'LineWidth', lineThick)
set(hH2.mainLine, 'LineWidth', lineThick)
xlabel('Time (s)', 'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold')
ylabel('f0 (cents)', 'FontName', fontN, 'FontSize', axisLSize, 'FontWeight', 'bold')
title('Offset of Perturbation', 'FontName', fontN, 'FontSize', titleFSize, 'FontWeight', 'bold')
axis(limits); box off

set(gca,'XTickLabel', {'-0.5' '0' '0.5' '1.0'},...
        'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold',...
        'YAxisLocation', 'right');

sup = suptitle({curSess; ['AudFB: ' AudFB]; ['f0: ' num2str(f0b) 'Hz, ' num2str(numPT) ' trial(s)']});
set(sup, 'FontName', fontN,...
         'FontSize', titleFSize,...
         'FontWeight','bold')
     
legend([mH.mainLine hH.mainLine],{'Microphone', 'Headphones'},...
            'Position', [0.8 0.30 0.1 0.1],...
            'Box', 'off',...
            'Edgecolor', [1 1 1],...
            'FontName', fontN,...
            'FontSize', legAnnoFSize,...
            'FontWeight', 'bold');


plots = {'AudResp_MeanTrial'};
for i = 1:length(plots)
    plTitle = [curSess '_' plots{i} '.jpg'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end
end