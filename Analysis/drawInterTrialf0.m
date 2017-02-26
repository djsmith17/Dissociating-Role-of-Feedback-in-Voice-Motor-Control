function drawInterTrialf0(time, meanTrialf0_St, meanTrialf0_Sp, limits, counts, meanTrialf0b, curExp, curRecording, plotFolder)
plotpos = [200 100];
plotdim = [1300 500];
AveInterTrialNHR = figure('Color', [1 1 1]);
set(AveInterTrialNHR, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

curExp(strfind(curExp, '_')) = ' ';

dottedStartx = [0.5 0.5];
dottedy = [-300 300];

ha = tight_subplot(1,2,[0.1 0.05],[0.12 0.15],[0.05 0.03]);

axes(ha(1))
errorbar(time, meanTrialf0_St(:,1,1), meanTrialf0_St(:,2,1), 'blue', 'LineWidth',2) %Unperturbed
hold on
errorbar(time, meanTrialf0_St(:,1,2), meanTrialf0_St(:,2,2), 'black', 'LineWidth',2) %Perturbed
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')

title('Onset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
axis(limits); box off

set(gca,'XTickLabel',{'-0.5' '0' '0.5' '1.0'},...
        'FontSize', 16,...
        'FontWeight','bold')

l0 = legend([num2str(counts(1)) ' Control Trials'], [num2str(counts(2)) ' Perturb Trials']); set(l0,'box', 'off','FontSize', 14, 'FontWeight', 'bold');

axes(ha(2))
errorbar(time, meanTrialf0_Sp(:,1,1), meanTrialf0_Sp(:,2,1), 'blue', 'LineWidth',2)  %Unperturbed
hold on
errorbar(time, meanTrialf0_Sp(:,1,2), meanTrialf0_Sp(:,2,2), 'black', 'LineWidth',2) %Perturbed
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')

title('Offset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
axis(limits); box off
set(gca,'XTickLabel', {'-0.5' '0' '0.5' '1.0'},...
        'FontSize', 16,...
        'FontWeight','bold');

suptitle({[curExp ': Mic Recording']; [curRecording '   f0: ' num2str(meanTrialf0b) 'Hz']})

plots = {'InterTrialf0'};
for i = 1:length(plots)
    plTitle = [curRecording '_' plots{i} '.png'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end
end