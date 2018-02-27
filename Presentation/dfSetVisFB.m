function [anMsr, H1, H2, H3, fbLines, LoudRec, visTrig] = dfSetVisFB(targRMS, bounds)
%Overlays for the experiment.
%anMsr = annotation Measurements

monitorSize = get(0, 'Monitor');
if size(monitorSize,1) == 1
    figPosition = [1 200 monitorSize(3) monitorSize(4)-200];
elseif size(monitorSize,1) == 2
    figPosition = [monitorSize(2,1) monitorSize(2,2) monitorSize(1,3) monitorSize(2,4)];
end
anMsr.winPos = figPosition;

%Assuming that perfect RMS lies at 50% of the screen, and the acceptable
%upper and lower bound are set at 45% and 55% of the screen, the edge of
%the screen will correspond to these variables accordingly. This will shift
%with differences in participants varying RMS values
anMsr.targRMS   = targRMS; %Ideally this is between 40-70 dB
anMsr.hRMSrange = bounds*10; %half the possible range of RMS, where each 5% mark represents 'bounds' dB
anMsr.vTopAmp   = anMsr.targRMS + anMsr.hRMSrange; %dB 
anMsr.vBotAmp   = anMsr.targRMS - anMsr.hRMSrange; %dB 

%Starts with some measurements
anMsr.bMar      = 0.05; %Bottom Margin
anMsr.minH      = 0.4;  %From the margin. This will need to change
anMsr.maxH      = 0.5;  %From the margin. This will need to change
anMsr.drawMinH  = anMsr.bMar + anMsr.minH; %45%
anMsr.drawMaxH  = anMsr.bMar + anMsr.maxH; %55% Goal value is 50%

anMsr.lMar      = 0.15; %Left Margin
anMsr.lineWidth = 0.1;  %Arbitrary-Looks good at the moment
anMsr.drawMinW  = anMsr.lMar + anMsr.lineWidth;

anMsr.recWidth  = 0.05; %Arbitrary-Looks good at the moment
anMsr.recHeight = 0.05; %Start at this default
anMsr.recXSt    = anMsr.lMar + (anMsr.lineWidth/2 - anMsr.recWidth/2);
anMsr.recYSt    = anMsr.bMar; %Start at the margin

anMsr.minLx = [anMsr.lMar anMsr.drawMinW];     %First X and Last X
anMsr.minLy = [anMsr.drawMinH anMsr.drawMinH]; %First Y and Last Y
anMsr.maxLx = [anMsr.lMar anMsr.drawMinW];     %First X and Last X
anMsr.maxLy = [anMsr.drawMaxH anMsr.drawMaxH]; %First Y and Last Y

anMsr.recPos = [anMsr.recXSt anMsr.recYSt anMsr.recWidth anMsr.recHeight];

anMsr.r        = 60;
anMsr.visTrigX = 0;   
anMsr.visTrigY = 0.02;
anMsr.visTrigW = round(100*anMsr.r/anMsr.winPos(3))/100;
anMsr.visTrigH = round(100*anMsr.r/anMsr.winPos(4))/100;

anMsr.visTrigPos = [anMsr.visTrigX anMsr.visTrigY anMsr.visTrigW anMsr.visTrigH];

%%%%%%
figure1 = figure('NumberTitle','off','Color',[0 0 0],'Position', anMsr.winPos,'MenuBar','none');

%Plus Sign
H1 = annotation(figure1,'textbox',[0.425 0.425 0.15 0.15],...
                        'Color',[1 1 1],...
                        'String',{'+'},...
                        'LineStyle','none',...
                        'HorizontalAlignment','center',...
                        'VerticalAlignment','middle',...
                        'FontWeight','bold',...
                        'FontSize',160,...
                        'FontName','Arial',...
                        'FitBoxToText','off',...
                        'EdgeColor','none',...
                        'BackgroundColor',[0 0 0],...
                        'Visible','off');

%EEE Directions
H2 = annotation(figure1,'textbox',[0.39 0.42 0.22 0.16],...
                        'Color',[1 1 1],...
                        'String',{'eee'},...
                        'LineStyle','none',...
                        'HorizontalAlignment','center',...
                        'VerticalAlignment','middle',...
                        'FontWeight','bold',...
                        'FontSize',160,...
                        'FontName','Arial',...
                        'FitBoxToText','off',...
                        'EdgeColor','none',...
                        'BackgroundColor',[0 0 0],...
                        'visible','off');

%Ready Note                    
H3 = annotation(figure1,'textbox', [0.3 0.35 0.4 0.3],...
                        'Color',[1 1 1],...
                        'String',{'READY'},...
                        'LineStyle','none',...
                        'HorizontalAlignment','center',...
                        'VerticalAlignment','middle',...
                        'FontSize',130,...
                        'FontName','Arial',...
                        'FitBoxToText','off',...
                        'EdgeColor','none',...
                        'BackgroundColor',[0 0 0],...
                        'Visible','on');
                    
visTrig = annotation(figure1, 'rectangle', anMsr.visTrigPos,...
                              'Color',[1 1 1],...
                              'FaceColor', [1 1 1],...
                              'visible', 'off');
                    
LoudRec = annotation(figure1,'rectangle', anMsr.recPos,...
                             'Color',[0 1 0],...
                             'LineStyle','none',...
                             'FaceColor',[0 1 0],...
                             'visible','off');
                     
LoudminLine = annotation(figure1,'line', anMsr.minLx, anMsr.minLy,...
                                 'Color',[1 1 1],...
                                 'LineStyle', '--',...
                                 'LineWidth', 1, ...
                                 'visible','off');
                         
LoudmaxLine = annotation(figure1,'line', anMsr.maxLx, anMsr.maxLy,...
                                 'Color',[1 1 1],...
                                 'LineStyle', '--',...
                                 'LineWidth', 1, ...
                                 'visible','off');
% ticks = 0:0.05:1;
% for i = ticks
%     tickLine = annotation(figure1,'line', [0.05 0.07], [i i],...
%                              'Color',[1 1 1],...
%                              'LineWidth', 1, ...
%                              'visible','on');
%     tickName = annotation(figure1, 'textbox', [0.07 i 0.02 0.02],...
%                                     'string', num2str(i),...
%                                     'color', 'white',...
%                                     'visible','on');
% end
                
fbLines = [LoudminLine LoudmaxLine];   
end