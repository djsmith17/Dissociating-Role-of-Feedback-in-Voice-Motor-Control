function drawInflationResultMetrics(onset, ir)

figure('Color', [1 1 1])
plot([ir.tAtOnset ir.tAtOnset], [-300 300], 'k--')
hold on
plot([-0.6 1.1], [ir.vAtOnset, ir.vAtOnset], 'r--')

plot(ir.time, onset, 'k', 'LineWidth', 2)
xlabel('Time (s)')
ylabel('f0 (Cents)')
hold on

plot(ir.tAtOnset, ir.vAtOnset, 'r*', 'LineWidth', 2)
plot(ir.tAtMin, ir.vAtMin, 'm*', 'LineWidth', 2);
plot(ir.tAtResp, ir.vAtResp, 'g*', 'LineWidth', 2);

axis([ir.time(1) ir.time(end) (min(onset) - 20) (max(onset) + 20)])
box off

SMArrowX = [(ir.tAtOnset -0.2), (ir.tAtOnset -0.2)];
SMArrowY = [ir.vAtOnset ir.vAtMin];

% plot(SMArrowX, SMArrowY, 'k--o')

stimMag = round(ir.stimMag, 1);
respMag = round(ir.respMag, 1);
respPer = round(ir.respPer);

annoSM = ['SM: ' num2str(stimMag) ' cents'];
annoRM = ['RM: ' num2str(respMag) ' cents'];
annoRP = ['RP: ' num2str(respPer) '%'];

respVarAnno = annotation('textbox', [0.8 0.8 0.21 0.15],...
                         'string',{annoSM
                                   annoRM
                                   annoRP},...
                                   'LineStyle', 'none',...
                         'FontWeight','bold');
end