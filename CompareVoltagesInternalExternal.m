plotpos = [10 100];
plotdim = [1600 600];

compV = figure('Color', [1 1 1]);
set(compV, 'Position', [plotpos plotdim],'PaperPositionMode','auto');

numTrial = niAn.numTrial;

ha = tight_subplot(1,numTrial,[0.05 0.05],[0.1 0.01],[0.05 0.05]);
for i = 1:numTrial
    axes(ha(i))
    
    plot(niAn.time, niAn.sensorFN(:,i), 'b')
    hold on
    plot(niAn.time, niAn.sensorP(:,i), 'r')
    
    if i == 1
        xlabel('Time (s)')
        ylabel('Voltage (V)')
    end
    
    diff = niAn.sensorP(:,i) - niAn.sensorFN(:,i);
    plot(niAn.time, diff, 'g')
    
    axis([0 4 -0.05 5])
    box off
end

suptitle('SpringNewBoxCalib DS5')

legend('Internal Voltage', 'External Voltage', 'Ext - Int');