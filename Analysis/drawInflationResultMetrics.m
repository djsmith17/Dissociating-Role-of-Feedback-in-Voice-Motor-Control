function drawInflationResultMetrics(ir, arrows, vals, toggle)
% Optimized for DRF18 SF3 mean trace results (default).
% Can also be toggled for DRF13 mean auditory trace (toggle = 1).

plotPos = [0 40];
plotDim = [1000 667];

fontN        = 'Arial';
lineThick    = 3;
axisLSize    = 20;
arrowTextSize = 20;

onsetC = [217,95,2]/255;
minC   = [117,112,179]/255;
respC  = [27,158,119]/255;
figColor = [1 1 1];
bColor = [0.8 0.8 0.8];

if toggle == 1
    % Auditory
    lowf0PixelY  = 0.22;
    lastf0PixelY = 0.455;
    lowf0PixelYLow = lowf0PixelY - 0.06;
    
    SMArrowX = [0.31 0.31];
    SMArrowY = [0.80 lowf0PixelY];
    SMAnnoBox = [0.165 0.32 0.1 0.1];
    
    RMArrowX = [0.854 0.854];
    RMArrowY = [lowf0PixelY lastf0PixelY];
    RMAnnoBox = [(RMArrowX(1)-0.130) 0.26 0.1 0.1];
    
    arrowXEnd = 0.47;
    tMArrowX = [0.39 arrowXEnd];
    tMArrowY = [lowf0PixelYLow lowf0PixelYLow];
    tMAnnoBox = [(arrowXEnd + .09) (lowf0PixelYLow-0.076) 0.1 0.1];
    expType = 'Auditory';
else
    % Default
    lowf0PixelY  = 0.34;
    lastf0PixelY = 0.685;
    lowf0PixelYLow = lowf0PixelY - 0.1;
    
    SMArrowX = [0.31 0.31];
    SMArrowY = [0.813 lowf0PixelY];
    SMAnnoBox = [0.17 0.37 0.1 0.1];
    
    RMArrowX = [0.854 0.854];
    RMArrowY = [lowf0PixelY lastf0PixelY];
    RMAnnoBox = [(RMArrowX(1)-0.13) 0.41 0.1 0.1];
    
    arrowXEnd = 0.45;
    tMArrowX = [0.39 arrowXEnd];
    tMArrowY = [lowf0PixelYLow lowf0PixelYLow];
    tMAnnoBox = [(arrowXEnd + 0.03) (lowf0PixelYLow-0.05) 0.1 0.1];
    expType = 'Laryngeal';
end

DVFig = figure('Color', figColor);
set(DVFig, 'Position', [plotPos plotDim]);

plot([ir.tAtOnset ir.tAtOnset], [-300 300], 'k--', 'LineWidth', 2.5)
hold on
plot([-0.6 1.1], [ir.vAtOnset, ir.vAtOnset], 'r--', 'LineWidth', 2.5)

area([ir.tAtRespRange(1) ir.tAtRespRange(2)], [min(ir.vAtRespRange)-3 min(ir.vAtRespRange)-3], max(ir.vAtRespRange)+3, 'FaceColor', respC, 'FaceAlpha', 0.25, 'EdgeAlpha', 0, 'ShowBaseline', 'off')

% Draw the pitch trace
plot(ir.time, ir.onset, 'k', 'LineWidth', lineThick)
xlabel('Time (s)')
ylabel('{\it f}_o (cents)')
hold on

% Draw the markers
plot(ir.tAtOnset, ir.vAtOnset, 'Color', onsetC, 'Marker', 's', 'MarkerFaceColor',onsetC, 'MarkerSize', 20)
plot(ir.tAtMin, ir.vAtMin, 'Color', minC, 'Marker', '^','MarkerFaceColor', minC, 'MarkerSize', 20);
plot(ir.tAtResp, ir.vAtResp, 'Color', respC, 'Marker', 'o', 'MarkerFaceColor', respC,'MarkerSize', 20);

% set the plot characteristics
axis([ir.time(1) ir.time(end) -120 20])
box off

set(gca,'LineWidth', 2,...
        'YTick', -150:50:0,...
        'FontName', fontN,...
        'FontSize', axisLSize,...
        'FontWeight','bold',...
        'Color', figColor)

% Toggle to draw arrows
if arrows == 1
 
    annotation('arrow', SMArrowX, SMArrowY, 'LineWidth', 5, 'HeadWidth', 20)
    annotation('textbox', SMAnnoBox,...
                         'string', {'Stimulus', 'Magnitude'},...
                         'LineStyle', 'none',...
                         'HorizontalAlignment', 'center',...
                         'FontName', fontN,...
                         'FontWeight','bold',...
                         'FontSize',arrowTextSize)

    annotation('arrow', RMArrowX, RMArrowY, 'LineWidth', 5, 'HeadWidth', 20)
    annotation('textbox', RMAnnoBox,...
                         'string', {'Response', 'Magnitude'},...
                         'LineStyle', 'none',...
                         'HorizontalAlignment', 'center',...
                         'FontWeight','bold',...
                         'FontName', fontN,...
                         'FontSize',arrowTextSize)
    
    annotation('arrow', tMArrowX, tMArrowY, 'LineWidth', 5, 'HeadWidth', 20)
    annotation('textbox', tMAnnoBox,...
                         'string', {'Response Latency'},...
                         'LineStyle', 'none',...
                         'HorizontalAlignment', 'center',...
                         'FontWeight','bold',...
                         'FontName', fontN,...
                         'FontSize',arrowTextSize)
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
% path = 'E:\Desktop';
path  = 'C:\Users\djsmith\Desktop';
fileName = [expType 'DependentVariablesManuscript.png'];
fullpath = fullfile(path, fileName);
export_fig(fullpath, '-r300')
end