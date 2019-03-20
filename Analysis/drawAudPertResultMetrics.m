function drawAudPertResultMetrics(ir, arrows, vals)
% This has been optimized to diagram DRF4 mean trace results.
% Hopefully in the future, the arrows will 

dim     = 50;
plotPos = [500 200];
plotDim = [dim*16 dim*9];

fontN        = 'Arial';
lineThick    = 3;
axisLSize    = 20;
legAnnoFSize = 20;

onsetC = [252,141,89]/255;
respC  = [255,255,191]/255;
figColor = [1 1 1];
bColor = [0.8 0.8 0.8];

DVFig = figure('Color', figColor);
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
plot(ir.tAtResp, ir.vAtResp, 'Color', respC, 'Marker', 'o', 'MarkerFaceColor', respC,'MarkerSize', 10);

axis([ir.time(1) ir.time(end) (min(ir.onset) - 20) (max(ir.onset) + 20)])
box off

set(gca,'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold',...
        'Color', figColor)

if arrows == 1
    annotation('textbox', [0.415 0.8 0.1 0.1],...
                         'string', 'StimMag = 100 cents',...
                         'LineStyle', 'none',...
                         'FontName', fontN,...
                         'FontWeight','bold',...
                         'FontSize',legAnnoFSize)

    RMArrowX = [0.905 0.905];
    RMArrowY = [0.30 0.80];
    annotation('arrow', RMArrowX, RMArrowY, 'LineWidth', 2)
    annotation('textbox', [0.74 0.35 0.1 0.1],...
                         'string', 'RespMag',...
                         'LineStyle', 'none',...
                         'FontWeight','bold',...
                         'FontName', fontN,...
                         'FontSize',legAnnoFSize)
end

if vals == 1
    stimMag = round(ir.stimMag, 2);
    respMag = round(ir.respMag, 2);
    respPer = round(ir.respPer, 2);

    annoSM = ['SM: ' num2str(stimMag) ' cents'];
    annoRM = ['RM: ' num2str(respMag) ' cents'];
    annoRP = ['RP: ' num2str(respPer) '%'];

    respVarAnno = annotation('textbox', [0.8 0.8 0.21 0.15],...
                             'string',{annoSM
                                       annoRM
                                       annoRP},...
                             'FontName', fontN,...
                             'FontSize', 12,...
                             'LineStyle', 'none',...
                             'FontWeight','bold');
end
                     
path = 'E:\Desktop';
% path  = 'C:\Users\djsmith\Desktop';
fileName = 'DependentVariablesManuscript.jpg';
fullpath = fullfile(path, fileName);
export_fig(fullpath, '-nocrop')
end