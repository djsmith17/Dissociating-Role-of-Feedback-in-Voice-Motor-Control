function drawIntraTrialf0(time, Trialf0_St, Trialf0_Sp,  limits, trialType, meanTrialf0b, curExp, curRecording, plotFolder)
plotpos = [200 100];
plotdim = [1300 500];
IntraTrialf0 = figure('Color', [1 1 1]);
set(IntraTrialf0, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

curExp(strfind(curExp, '_')) = ' ';

dottedStartx = [0.5 0.5];
dottedy      = [-300 300];

ha = tight_subplot(1,2,[0.1 0.05],[0.1 0.03],[0.05 0.03]);

axes(ha(1))
plot(time, Trialf0_St(:,1),'b')
hold on
plot(time, Trialf0_St(:,2),'r')
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')

title('Onset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
axis(limits); box off

set(gca,'XTickLabel',{'-0.5' '0' '0.5' '1.0'},...
        'FontSize', 16,...
        'FontWeight','bold')

l0 = legend('Microphone', 'Headphones'); 
set(l0,'box', 'off','FontSize', 14, 'FontWeight', 'bold');

axes(ha(2))
plot(time, Trialf0_Sp(:,1),'b')
hold on
plot(time, Trialf0_Sp(:,2),'r')
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')

title('Offset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
axis(limits); box off

set(gca,'XTickLabel',{'-0.5' '0' '0.5' '1.0'},...
        'FontSize', 16,...
        'FontWeight','bold')

Types ={'Unperturbed'; 'Perturbed'};

suptitle({[curExp ': ' Types{trialType+1} ' Trial']; [curRecording '   f0: ' num2str(meanTrialf0b) 'Hz']})

plots = {'IntraTrialf0_AudResp'};
for i = 1:length(plots)
    plTitle = [curRecording '_' plots{i} '.png'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end
end