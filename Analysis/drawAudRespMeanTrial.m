function drawAudRespMeanTrial(res, plotFolder)
%Showing the mean 

curSess          = res.curSess;
f0b              = round(10*res.f0b)/10;
AudFB            = res.AudFB;
numPT            = res.numPertTrialsPP;

time             = niRes.secTime;
meanf0MicOnset  = niRes.audioMf0MeanPert(:,1);
CIf0MicOnset    = niRes.audioMf0MeanPert(:,2);
meanf0MicOffset = niRes.audioMf0MeanPert(:,3);
CIf0MicOffset   = niRes.audioMf0MeanPert(:,4);

meanf0HeadOnset  = niRes.audioHf0MeanPert(:,1);
CIf0HeadOnset    = niRes.audioHf0MeanPert(:,2);
meanf0HeadOffset = niRes.audioHf0MeanPert(:,3);
CIf0HeadOffset   = niRes.audioHf0MeanPert(:,4);
limits           = niRes.limitsAmean;

plotpos = [10 100];
plotdim = [1600 600];
MeanTrialAudResp = figure('Color', [1 1 1]);
set(MeanTrialAudResp, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

dottedStartx = [0 0];
dottedy      = [-300 300];

ha = tight_subplot(1,2,[0.1 0.05],[0.12 0.15],[0.08 0.08]);

%Onset of Perturbation
axes(ha(1))
mH = shadedErrorBar(time, meanf0MicOnset, CIf0MicOnset, 'b', 1); %Pertrubed Microphone
hold on
hH = shadedErrorBar(time, meanf0HeadOnset, CIf0HeadOnset, 'r', 1); %Perturbed Headphones
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')

title('Onset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
axis(limits); box off

set(gca,'XTickLabel',{'-0.5' '0' '0.5' '1.0'},...
        'FontSize', 16,...
        'FontWeight','bold')

l0 = legend([mH.mainLine hH.mainLine], 'Microphone', 'Headphones'); 
set(l0,'box', 'off','FontSize', 14, 'FontWeight', 'bold');

%Offset of Perturbation
axes(ha(2))
shadedErrorBar(time, meanf0MicOffset, CIf0MicOffset, 'b', 1)  %Perturbed Microphone
hold on
shadedErrorBar(time, meanf0HeadOffset, CIf0HeadOffset, 'r', 1) %Perturbed Headphones
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

plots = {'InterTrialf0AudResp'};
for i = 1:length(plots)
    plTitle = [curSess '_' plots{i} '.jpg'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end
end