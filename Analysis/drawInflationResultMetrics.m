function drawInflationResultMetrics(irAll, trial)

ir = irAll(trial);

dim     = 50;
plotPos = [500 200];
plotDim = [dim*16 dim*9];

fontN        = 'Arial';
lineThick    = 3;
axisLSize    = 20;
legAnnoFSize = 20;

onsetC = [252,141,89]/255;
minC   = [145,191,219]/255;
respC  = [255,255,191]/255;
bColor = [0.8 0.8 0.8];

DVFig = figure('Color', bColor);
set(DVFig, 'Position', [plotPos plotDim]);

plot([ir.tAtOnset ir.tAtOnset], [-300 300], 'k--', 'LineWidth', 1.5)
hold on
plot([-0.6 1.1], [ir.vAtOnset, ir.vAtOnset], 'r--', 'LineWidth', 1.5)

% Draw the pitch trace
plot(ir.time, ir.onset, 'k', 'LineWidth', lineThick)
xlabel('Time (s)')
ylabel('f0 (Cents)')
hold on

plot(ir.tAtOnset, ir.vAtOnset, 'Color', onsetC, 'Marker', 'o', 'MarkerFaceColor',onsetC, 'MarkerSize', 10)
plot(ir.tAtMin, ir.vAtMin, 'Color', minC, 'Marker', 'o','MarkerFaceColor', minC,  'MarkerSize', 10);
plot(ir.tAtResp, ir.vAtResp, 'Color', respC, 'Marker', 'o', 'MarkerFaceColor', respC,'MarkerSize', 10);

axis([ir.time(1) ir.time(end) (min(ir.onset) - 20) (max(ir.onset) + 20)])
box off

set(gca,'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold',...
        'Color', bColor)

SMArrowX = [0.33 0.33];
SMArrowY = [0.8 0.24];
annotation('arrow', SMArrowX, SMArrowY, 'LineWidth', 2)
annotation('textbox', [0.175 0.38 0.1 0.1],...
                     'string', 'StimMag',...
                     'LineStyle', 'none',...
                     'FontName', fontN,...
                     'FontWeight','bold',...
                     'FontSize',legAnnoFSize)

RMArrowX = [0.905 0.905];
RMArrowY = [0.24 0.65];
annotation('arrow', RMArrowX, RMArrowY, 'LineWidth', 2)
annotation('textbox', [0.74 0.3 0.1 0.1],...
                     'string', 'RespMag',...
                     'LineStyle', 'none',...
                     'FontWeight','bold',...
                     'FontName', fontN,...
                     'FontSize',legAnnoFSize)

stimMag = round(ir.stimMag, 1);
respMag = round(ir.respMag, 1);
respPer = round(ir.respPer);

annoSM = ['SM: ' num2str(stimMag) ' cents'];
annoRM = ['RM: ' num2str(respMag) ' cents'];
annoRP = ['RP: ' num2str(respPer) '%'];

% respVarAnno = annotation('textbox', [0.8 0.8 0.21 0.15],...
%                          'string',{annoSM
%                                    annoRM
%                                    annoRP},...
%                          'FontName', fontN,...
%                          'FontSize', legAnnoFSize,...
%                          'LineStyle', 'none',...
%                          'FontWeight','bold');
                     
path = 'E:\Desktop';
% path  = 'C:\Users\djsmith\Desktop';
fileName = 'DependentVariablesManuscript.jpg';
fullpath = fullfile(path, fileName);
export_fig(fullpath, '-nocrop')
end