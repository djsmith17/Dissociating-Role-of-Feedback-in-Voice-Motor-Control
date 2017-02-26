function drawIntraTrialf0(time, plotf0pts_St, plotf0pts_Sp, pert, limits, curRecording, k, plotFolder)
plotpos = [200 100];
plotdim = [1300 500];
InterTrialNHR = figure('Color', [1 1 1]);
set(InterTrialNHR, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

dottedStartx = [0.5 0.5];
dottedy = [-300 300];

ha = tight_subplot(1,2,[0.1 0.05],[0.1 0.03],[0.05 0.03]);

axes(ha(1))
plot(time, plotf0pts_St(:,1))
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)'); ylabel('f0 (Hz)')

title('Onset of Perturbation', 'FontSize', 10, 'FontWeight', 'bold')
axis(limits); box off

axes(ha(2))
plot(time, plotf0pts_Sp(:,1))
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)'); ylabel('f0 (Hz)')

title('Offset of Perturbation', 'FontSize', 10, 'FontWeight', 'bold')
axis(limits); box off

suptitle([curRecording ' Trial #' num2str(k)])

if pert == 0
    legend('Unperturbed')
else
    legend('Perturbed')
end

plots = {'IntraTrialf0'};
for i = 1:length(plots)
    plTitle = [curRecording '_' plots{i} '.png'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end
end