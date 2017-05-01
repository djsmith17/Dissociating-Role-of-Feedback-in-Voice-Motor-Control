function drawInterTrialf0(time, meanTrialf0_St, meanTrialf0_Sp, limits, counts, meanTrialf0b, curExp, curRecording, plotFolder)
plotpos = [200 100];
plotdim = [1300 500];
InterTrialf0 = figure('Color', [1 1 1]);
set(InterTrialf0, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

curExp(strfind(curExp, '_')) = ' ';
curRecording(strfind(curRecording, '_')) = ' ';

dottedStartx = [0.5 0.5];
dottedy      = [-300 300];

ha = tight_subplot(1,2,[0.1 0.05],[0.12 0.15],[0.05 0.03]);

axes(ha(1))
uH = shadedErrorBar(time, meanTrialf0_St(:,1,1), meanTrialf0_St(:,2,1), 'b', 1); %Unperturbed
hold on
pH = shadedErrorBar(time, meanTrialf0_St(:,1,2), meanTrialf0_St(:,2,2), 'r', 1); %Perturbed
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')

title('Onset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
axis(limits); box off

set(gca,'XTickLabel',{'-0.5' '0' '0.5' '1.0'},...
        'FontSize', 16,...
        'FontWeight','bold')

l0 = legend([uH.mainLine pH.mainLine],[num2str(counts(1)) ' Control Trials'], [num2str(counts(2)) ' Perturb Trials']); 
set(l0,'box', 'off','FontSize', 14, 'FontWeight', 'bold');

axes(ha(2))
shadedErrorBar(time, meanTrialf0_Sp(:,1,1), meanTrialf0_Sp(:,2,1), 'b', 1)  %Unperturbed
hold on
shadedErrorBar(time, meanTrialf0_Sp(:,1,2), meanTrialf0_Sp(:,2,2), 'r', 1) %Perturbed
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')

title('Offset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
axis(limits); box off
set(gca,'XTickLabel', {'-0.5' '0' '0.5' '1.0'},...
        'FontSize', 16,...
        'FontWeight','bold',...
        'YAxisLocation', 'right');

suptitle({[curExp ': Mic Recording']; [curRecording '   f0: ' num2str(meanTrialf0b) 'Hz']})

plots = {'InterTrialf0'};
for i = 1:length(plots)
    plTitle = [curRecording '_' plots{i} '.jpg'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end
end