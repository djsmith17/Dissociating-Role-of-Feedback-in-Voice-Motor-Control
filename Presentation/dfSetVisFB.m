function [anMsr, H1, H2, H3, fbLines, LoudRec, visTrig] = dfSetVisFB(curSess, targRMS, bounds)
% [anMsr, H1, H2, H3, fbLines, LoudRec, visTrig] = dfSetVisFb(targRMS,  bounds) 
% creates a figure, which becomes the means for trial progression
% for the SFPeturb or AFPerturb experiment.
%
% targRMS: Baseline RMS (dB) of the subject and the target for trials
%          where visual feedback about their loudness is presented. 
% bounds:  The acceptable RMS bounds (dB) where the participant is considered
%          phonating at a reasonable level.
%
% anMsr:   Struc of annotation measurements. The collection of the 
%          positions for all the annotations on the figure. Annotations are
%          just objects on the figure that can be moved around.
% H1:      '+' annotation
% H2:      'EEE' annotation
% H3:      'Ready' annotation
% fBLines: The lines representing the bounds of the acceptable RMS range
% LoudRec: Bar annotation representing the amount of loudness the
%          participant had on the last trial
% visTrig: The trigger annotation that is used to trigger the triggerbox

monitorSize = get(0, 'Monitor');
numMon = size(monitorSize, 1);

if numMon == 1
    W = monitorSize(3);
    H = monitorSize(4);
    W2 = W/2;
    H2 = H/2;
    XPos = W2;
    YPos = 50;
    
    figPosition = [XPos YPos W2 H2];
elseif numMon == 2
    [~, mon] = max(monitorSize(:,1));
    
    figPosition = [monitorSize(mon,1) monitorSize(mon,2) monitorSize(1,3) monitorSize(1,4)];
end
anMsr.winPos = figPosition;

% Ready Annotation Dim
anMsr.rdAnoD = [700 300];
anMsr.rdAnoW = round(anMsr.rdAnoD(1)/anMsr.winPos(3), 2); 
anMsr.rdAnoH = round(anMsr.rdAnoD(2)/anMsr.winPos(4), 2);
anMsr.rdAnoX = 0.5 - anMsr.rdAnoW/2;
anMsr.rdAnoY = 0.5 - anMsr.rdAnoH/2;
anMsr.rdAnoPos = [anMsr.rdAnoX anMsr.rdAnoY anMsr.rdAnoW anMsr.rdAnoH];

% Cue Annotation Dim
anMsr.cuAnoD = [250 150];
anMsr.cuAnoW = round(anMsr.cuAnoD(1)/anMsr.winPos(3), 2); 
anMsr.cuAnoH = round(anMsr.cuAnoD(2)/anMsr.winPos(4), 2);
anMsr.cuAnoX = 0.5 - anMsr.cuAnoW/2;
anMsr.cuAnoY = 0.5 - anMsr.cuAnoH/2;
anMsr.cuAnoPos = [anMsr.cuAnoX anMsr.cuAnoY anMsr.cuAnoW anMsr.cuAnoH];

% EEE Annotation Dim
anMsr.eeAnoD = [370 170];
anMsr.eeAnoW = round(anMsr.eeAnoD(1)/anMsr.winPos(3), 2); 
anMsr.eeAnoH = round(anMsr.eeAnoD(2)/anMsr.winPos(4), 2);
anMsr.eeAnoX = 0.5 - anMsr.eeAnoW/2;
anMsr.eeAnoY = 0.5 - anMsr.eeAnoH/2;
anMsr.eeAnoPos = [anMsr.eeAnoX anMsr.eeAnoY anMsr.eeAnoW anMsr.eeAnoH];

%Assuming that perfect RMS lies at 50% of the screen, and the acceptable
%upper and lower bound are set at 45% and 55% of the screen, the edge of
%the screen will correspond to these variables accordingly. This will shift
%with differences in participants varying RMS values
anMsr.targRMS   = targRMS;              % Ideally this is between 60-80 dB
anMsr.hRMSrange = bounds*10;            % half the possible range of RMS, where each 5% mark represents 'bounds' dB
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

% Feedback Rectangle dim measurements 
anMsr.recWidth  = 0.05; % Arbitrary-Looks good at the moment
anMsr.recHeight = 0.05; % Start at this default
anMsr.recXSt    = anMsr.lMar + (anMsr.lineWidth/2 - anMsr.recWidth/2);
anMsr.recYSt    = anMsr.bMar; % Start at the margin

% Feedback Rectangle: Position and Size
anMsr.recPos = [anMsr.recXSt anMsr.recYSt anMsr.recWidth anMsr.recHeight];

% Feedback RMS Bounds Dim
anMsr.minLx = [anMsr.lMar anMsr.drawMinW];     %First X and Last X
anMsr.minLy = [anMsr.drawMinH anMsr.drawMinH]; %First Y and Last Y
anMsr.maxLx = [anMsr.lMar anMsr.drawMinW];     %First X and Last X
anMsr.maxLy = [anMsr.drawMaxH anMsr.drawMaxH]; %First Y and Last Y

% Trigger dim measurements
anMsr.r        = 60;                                % Trigger side L (pixels)
anMsr.visTrigX = 0;                                 % Trigger X Pos
anMsr.visTrigY = 0.02;                              % Trigger Y Pos
anMsr.visTrigW = round(anMsr.r/anMsr.winPos(3), 2); % Trigger X Len
anMsr.visTrigH = round(anMsr.r/anMsr.winPos(4), 2); % Trigger Y Len
anMsr.visTrigPos = [anMsr.visTrigX anMsr.visTrigY anMsr.visTrigW anMsr.visTrigH];

% Current Subject Note Dim
anMsr.subjND = [170 30];
anMsr.subjNX = anMsr.visTrigX;
anMsr.subjNY = anMsr.visTrigY + anMsr.visTrigH + 0.02;
anMsr.subjNW = round(anMsr.subjND(1)/anMsr.winPos(3), 2); 
anMsr.subjNH = round(anMsr.subjND(2)/anMsr.winPos(4), 2); 
anMsr.subjNPos = [anMsr.subjNX anMsr.subjNY anMsr.subjNW anMsr.subjNH];

%%%%%%
VBFig = figure('NumberTitle', 'off', 'Color', [0 0 0], 'Position', anMsr.winPos, 'MenuBar', 'none');

%Plus Sign
H1 = annotation(VBFig,'textbox', anMsr.cuAnoPos,...
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
H2 = annotation(VBFig,'textbox', anMsr.eeAnoPos,...
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
H3 = annotation(VBFig,'textbox', anMsr.rdAnoPos,...
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

%Trigger for the triggerbox                    
visTrig = annotation(VBFig, 'rectangle', anMsr.visTrigPos,...
                              'Color',[1 1 1],...
                              'FaceColor', [1 1 1],...
                              'visible', 'off');

% Visual feedback bar                    
LoudRec = annotation(VBFig, 'rectangle', anMsr.recPos,...
                              'Color',[0 1 0],...
                              'LineStyle','none',...
                              'FaceColor',[0 1 0],...
                              'visible','off');

% Lower limit of acceptable RMS                     
LoudminLine = annotation(VBFig, 'line', anMsr.minLx, anMsr.minLy,...
                                  'Color',[1 1 1],...
                                  'LineStyle', '--',...
                                  'LineWidth', 1, ...
                                  'visible','off');

% Upper limit of acceptable RMS                              
LoudmaxLine = annotation(VBFig, 'line', anMsr.maxLx, anMsr.maxLy,...
                                  'Color',[1 1 1],...
                                  'LineStyle', '--',...
                                  'LineWidth', 1, ...
                                  'visible','off');
% Subject Note
curSessNote = annotation(VBFig,'textbox', anMsr.subjNPos,...
                                'Color',[1 1 1],...
                                'String', curSess,...
                                'LineStyle','none',...
                                'HorizontalAlignment','left',...
                                'VerticalAlignment','middle',...
                                'FontSize', 10,...
                                'FontName','Arial',...
                                'FitBoxToText','off',...
                                'EdgeColor','none',...
                                'BackgroundColor',[0 0 0],...
                                'Visible','on');
                             
fbLines = [LoudminLine LoudmaxLine];   
end