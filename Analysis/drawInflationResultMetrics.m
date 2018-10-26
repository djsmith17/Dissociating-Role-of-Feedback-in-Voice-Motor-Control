function drawInflationResultMetrics(onset, ir)

dim     = 50;
plotPos = [500 200];
plotDim = [dim*16 dim*9];

fontN        = 'Arial';
lineThick    = 2;
axisLSize    = 14;
legAnnoFSize = 14;

onsetC = 'r';
minC   = 'm';
respC  = 'g';

respFig = figure('Color', [1 1 1]);
set(respFig, 'Position', [plotPos plotDim]);

plot([ir.tAtOnset ir.tAtOnset], [-300 300], 'k--')
hold on
plot([-0.6 1.1], [ir.vAtOnset, ir.vAtOnset], 'r--')

% Draw the pitch trace
plot(ir.time, onset, 'k', 'LineWidth', lineThick)
xlabel('Time (s)')
ylabel('f0 (Cents)')
hold on

plot(ir.tAtOnset, ir.vAtOnset, 'Color', onsetC, 'Marker', '*', 'LineWidth', 2)
plot(ir.tAtMin, ir.vAtMin, 'Color', minC, 'Marker', '*', 'LineWidth', 2);
plot(ir.tAtResp, ir.vAtResp, 'Color', respC, 'Marker', '*', 'LineWidth', 2);

axis([ir.time(1) ir.time(end) (min(onset) - 20) (max(onset) + 20)])
box off

set(gca,'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold')

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
                         'FontName', fontN,...
                         'FontSize', legAnnoFSize,...
                         'LineStyle', 'none',...
                         'FontWeight','bold');
                     
path = 'E:\Desktop';
fileName = 'ExperimentalOutcomes.jpg';
fullpath = fullfile(path, fileName);
export_fig(fullpath)
end