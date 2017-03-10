function drawForceSensorSignal(time, meanTrialForce_St, meanTrialForce_Sp, limits, counts, curExp, curRecording, saveResultsDir)
%Good for seeing the whole signal
plotpos = [200 100];
plotdim = [1300 500];
ForceSensorS = figure('Color', [1 1 1]);
set(ForceSensorS, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

curExp(strfind(curExp, '_')) = ' ';

dottedStartx = [0.5 0.5];
dottedy      = [-300 300];

ha = tight_subplot(1,2,[0.1 0.05],[0.12 0.15],[0.05 0.03]);

axes(ha(1))
shadedErrorBar(time, meanTrialForce_St(:,1,2), meanTrialForce_St(:,2,2), 'b', 1) %Collar Sensor, Perturbed Trials
hold on
shadedErrorBar(time, meanTrialForce_St(:,3,2), meanTrialForce_St(:,4,2), 'r', 1) %Neck Sensor, Perturbed Trials
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('Voltage (V)', 'FontSize', 18, 'FontWeight', 'bold')

title('Onset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
axis(limits); box off

set(gca, 'YTick', 0:1:5,...
         'FontSize', 16,...
         'FontWeight', 'bold')
     
pltlgd = legend('Collar Sensor', 'Neck Sensor');
set(pltlgd, 'box', 'off',...
            'location', 'best'); 

axes(ha(2))
shadedErrorBar(time, meanTrialForce_Sp(:,1,2), meanTrialForce_Sp(:,2,2), 'b', 1) %Collar Sensor, Perturbed Trials
hold on
shadedErrorBar(time, meanTrialForce_Sp(:,3,2), meanTrialForce_Sp(:,4,2), 'r', 1)%Neck Sensor, Perturbed Trials
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('Voltage (V)', 'FontSize', 18, 'FontWeight', 'bold')

title('Offset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
axis(limits); box off

set(gca, 'YTick', 0:1:5,...
         'FontSize', 16,...
         'FontWeight', 'bold')

suptitle({[curExp ': Inflated Balloon-Collar Force Beta']; [curRecording ': ' num2str(counts(2)) ' Perturbed Trials']})

plTitle = [curRecording  '_ForceSensorAve.png'];     
saveFileName = fullfile(saveResultsDir, plTitle);
export_fig(saveFileName) 
end