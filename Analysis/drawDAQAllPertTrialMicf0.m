function drawDAQAllPertTrialMicf0(niRes, plotFolder)

curSess          = niRes.curSess;
f0b              = niRes.f0b;
numPT            = niRes.numPertTrials;

time             = niRes.timeA;
sigs             = niRes.audioMf0TrialPert;
limits           = [0 4 -150 150];

% meanf0PertOnset  = niRes.audioMf0MeanPert(:,1);
% CIf0PertOnset    = niRes.audioMf0MeanPert(:,2);
% meanf0PertOffset = niRes.audioMf0MeanPert(:,3);
% CIf0PertOffset   = niRes.audioMf0MeanPert(:,4);
% 
% meanf0ContOnset  = niRes.audioMf0MeanCont(:,1);
% CIf0ContOnset    = niRes.audioMf0MeanCont(:,2);
% meanf0ContOffset = niRes.audioMf0MeanCont(:,3);
% CIf0ContOffset   = niRes.audioMf0MeanCont(:,4);
% limits           = niRes.limitsAmean;

plotpos = [200 100];
plotdim = [1300 500];
InterTrialf0 = figure('Color', [1 1 1]);
set(InterTrialf0, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

curSess(strfind(curSess, '_')) = ' ';

dottedStartx = [0 0];
dottedy      = [-300 300];

ha = tight_subplot(1,2,[0.1 0.05],[0.12 0.15],[0.08 0.08]);

for i = 1:numPT
    subplot(2,5,i)
    plot(time, sigs(:,i))
    axis(limits)
end

% %Onset of Perturbation
% axes(ha(1))
% uH = shadedErrorBar(time, meanf0ContOnset, CIf0ContOnset, 'b', 1); %Unperturbed
% hold on
% pH = shadedErrorBar(time, meanf0PertOnset, CIf0PertOnset, 'r', 1); %Perturbed
% hold on
% plot(dottedStartx, dottedy,'k','LineWidth',4)
% xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')
% 
% title('Onset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
% axis(limits); box off
% 
% set(gca,'XTickLabel',{'-0.5' '0' '0.5' '1.0'},...
%         'FontSize', 16,...
%         'FontWeight','bold')
% 
% l0 = legend([uH.mainLine pH.mainLine],[num2str(numCT) ' Control Trials'], [num2str(numPT) ' Perturb Trials']); 
% set(l0,'box', 'off','FontSize', 14, 'FontWeight', 'bold');
% 
% %Offset of Perturbation
% axes(ha(2))
% shadedErrorBar(time, meanf0ContOffset, CIf0ContOffset, 'b', 1)  %Unperturbed
% hold on
% shadedErrorBar(time, meanf0PertOffset, CIf0PertOffset, 'r', 1) %Perturbed
% hold on
% plot(dottedStartx, dottedy,'k','LineWidth',4)
% xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')
% 
% title('Offset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
% axis(limits); box off
% set(gca,'XTickLabel', {'-0.5' '0' '0.5' '1.0'},...
%         'FontSize', 16,...
%         'FontWeight','bold',...
%         'YAxisLocation', 'right');
% 
% suptitle({[curSess ': Mic Recording']; ['   f0: ' num2str(f0b) 'Hz']})
% 
plots = {'IntraTrialf0DAQ'};
for i = 1:length(plots)
    plTitle = [curSess '_' plots{i} '.jpg'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end
end