function drawInflationResultMetrics(ir, arrows, vals, toggle)
% Optimized for DRF18 SF3 mean trace results (default).
% Can also be toggled for DRF13 mean auditory trace (toggle = 1).

plotPos = [-1485 40];
plotDim = [1000 667];

fontN        = 'Times New Roman';
lineThick    = 3;
axisLSize    = 30;
arrowTextSize = 26;

onsetC = [217,95,2]/255;
minC   = [117,112,179]/255;
respC  = [27,158,119]/255;
figColor = [1 1 1];
bColor = [0.8 0.8 0.8];

if toggle == 1
    % Auditory
    lowf0PixelY  = 0.26;
    lastf0PixelY = 0.48;
    lowf0PixelYLow = lowf0PixelY - 0.04;
    
    SMArrowX = [0.37 0.37];
    SMArrowY = [0.81 lowf0PixelY];
    SMAnnoBox = [0.21 0.35 0.1 0.1];
    
    RMArrowX = [0.856 0.856];
    RMArrowY = [lowf0PixelY lastf0PixelY];
    RMAnnoBox = [(RMArrowX(1)-0.140) 0.34 0.1 0.1];
    
    tMArrowX = [0.41 0.495];
    tMArrowY = [lowf0PixelYLow lowf0PixelYLow];
    tMAnnoBox = [0.59 (lowf0PixelYLow-0.07) 0.1 0.1];
    expType = 'Auditory';
else
    % Default
    lowf0PixelY  = 0.37;
    lastf0PixelY = 0.695;
    lowf0PixelYLow = lowf0PixelY - 0.06;
    
    SMArrowX = [0.37 0.37];
    SMArrowY = [0.82 lowf0PixelY];
    SMAnnoBox = [0.21 0.39 0.1 0.1];
    
    RMArrowX = [0.855 0.855];
    RMArrowY = [lowf0PixelY lastf0PixelY];
    RMAnnoBox = [(RMArrowX(1)-0.142) 0.42 0.1 0.1];
    
    arrowXEnd = 0.475;
    tMArrowX = [0.405 arrowXEnd];
    tMArrowY = [lowf0PixelYLow lowf0PixelYLow];
    tMAnnoBox = [(arrowXEnd + 0.03) (lowf0PixelYLow - 0.065) 0.1 0.1];
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
ylabel('{\it f}_o (Cents)')
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
path = 'E:\Desktop';
% path  = 'C:\Users\djsmith\Desktop';
fileName = [expType 'DependentVariablesManuscript.png'];
fullpath = fullfile(path, fileName);
export_fig(fullpath)
end