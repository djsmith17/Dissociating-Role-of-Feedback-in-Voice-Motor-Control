function drawAudPertResultMetrics(ir, arrows, vals)
% This has been optimized to diagram DRF19 AF2 mean trace results.
% Hopefully in the future, the arrows will 

plotPos = [10 40];
plotDim = [1485 1005];

fontN        = 'Times New Roman';
lineThick    = 3;
axisLSize    = 30;
legAnnoFSize = 30;

onsetC = [217,95,2]/255;
respC  = [27,158,119]/255;
figColor = [1 1 1];
bColor = [0.8 0.8 0.8];

DVFig = figure('Color', figColor);
set(DVFig, 'Position', [plotPos plotDim]);

plot([ir.tAtOnset ir.tAtOnset], [-300 300], 'k--', 'LineWidth', 2.5)
hold on
plot([-0.6 1.1], [ir.vAtOnset, ir.vAtOnset], 'r--', 'LineWidth', 2.5)

area([ir.tAtRespRange(1) ir.tAtRespRange(2)], [min(ir.vAtRespRange)-2 min(ir.vAtRespRange)-2], max(ir.vAtRespRange)+2, 'FaceColor', respC, 'FaceAlpha', 0.25, 'EdgeAlpha', 0, 'ShowBaseline', 'off')

% Draw the pitch trace
plot(ir.time, ir.onset, 'k', 'LineWidth', lineThick)
xlabel('Time (s)')
ylabel('{\it f}_o (Cents)')
hold on

% Draw the markers
plot(ir.tAtOnset, ir.vAtOnset, 'Color', onsetC, 'Marker', 's', 'MarkerFaceColor',onsetC, 'MarkerSize', 20)
plot(ir.tAtResp, ir.vAtResp, 'Color', respC, 'Marker', 'o', 'MarkerFaceColor', respC,'MarkerSize', 20);

% set the plot characteristics
axis([ir.time(1) ir.time(end) (min(ir.onset) - 10) (max(ir.onset) + 10)])
box off

set(gca,'LineWidth', 2,...
        'YTick', 0:50:100,...
        'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold',...
        'Color', figColor)

% Toggle to draw arrows
if arrows == 1

    RMArrowX = [0.853 0.853];
    RMArrowY = [0.25 0.765];
    annotation('arrow', RMArrowX, RMArrowY, 'LineWidth', 5, 'HeadWidth', 20)
    annotation('textbox', [(RMArrowX(1)-0.140) 0.33 0.1 0.1],...
                         'string', {'Response'; 'Magnitude'},...
                         'LineStyle', 'none',...
                         'HorizontalAlignment', 'center',...
                         'FontWeight','bold',...
                         'FontName', fontN,...
                         'FontSize',legAnnoFSize)
end

% Toggle to draw the measure values
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
  
% save the figure
path = 'E:\Desktop';
% path  = 'C:\Users\djsmith\Desktop';
fileName = 'AudDependentVariablesManuscript.png';
fullpath = fullfile(path, fileName);
export_fig(fullpath)
end