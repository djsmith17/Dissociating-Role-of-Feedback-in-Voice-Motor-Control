function drawAudResp_InterTrial(time, meanTrialf0_St, meanTrialf0_Sp, limits, counts, meanTrialf0b, curSess, curRec, plotFolder)
plotpos = [200 100];
plotdim = [1300 500];
InterTrialAudResp = figure('Color', [1 1 1]);
set(InterTrialAudResp, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

curSess(strfind(curSess, '_')) = ' ';
curRec(strfind(curRec, '_')) = ' ';

dottedStartx = [0.5 0.5];
dottedy      = [-300 300];

ha = tight_subplot(1,2,[0.1 0.05],[0.12 0.15],[0.08 0.08]);

axes(ha(1))
mH = shadedErrorBar(time, meanTrialf0_St(:,1), meanTrialf0_St(:,2), 'b', 1); %Pertrubed Microphone
hold on
hH = shadedErrorBar(time, meanTrialf0_St(:,3), meanTrialf0_St(:,4), 'r', 1); %Perturbed Headphones
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

axes(ha(2))
shadedErrorBar(time, meanTrialf0_Sp(:,1), meanTrialf0_Sp(:,2), 'b', 1)  %Pertrubed Microphone
hold on
shadedErrorBar(time, meanTrialf0_Sp(:,3), meanTrialf0_Sp(:,4), 'r', 1) %Perturbed Headphones
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')

title('Offset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
axis(limits); box off
set(gca,'XTickLabel', {'-0.5' '0' '0.5' '1.0'},...
        'FontSize', 16,...
        'FontWeight','bold',...
        'YAxisLocation', 'right');

suptitle({[curSess ': ' num2str(counts) ' Perturbed Trials']; [curRec '   f0: ' num2str(meanTrialf0b) 'Hz']})

plots = {'InterTrialf0_AudResp'};
for i = 1:length(plots)
    plTitle = [curSess '_' plots{i} '_' curRec '.jpg'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end
end