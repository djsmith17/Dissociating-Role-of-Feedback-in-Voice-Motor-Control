function [msrStr, annoStr] = dfSetVisFB(defMon, curSess, targRMS, bounds)
% [msrStr, annoStr] = dfSetVisFb(defMon, curSess, targRMS, bounds) 
% creates a figure, which becomes the visualization involved for trial progression
% in the behavioral experiments.
%
% INPUTS 
% defMon:  Monitor selection. A value of 1 uses the defualt monitor, A value of
%          2 selects the secondary monitor. There is currently not support for 
%          a number of monitors greater than 2.
% curSess: A string input that displays simple information about the current 
%          session being run. This is displayed as an annotation, so this can 
%          display any trial information that you want
% targRMS: Baseline RMS (dB) of the subject and the target for trials
%          where visual feedback about their intensity is presented. 
% bounds:  The acceptable RMS bounds (dB) where the participant is considered
%          vocalizing at the target level.
%
% OUTPUTS
% msrStr:  Struc of annotation measurements. The collection of the 
%          positions for all the annotations on the figure. Annotations are
%          just objects on the figure that can be moved around.
% annoStr: Struc of annotation elememts. The actual words, shapes and colors of 
%          items that are displauyed on the screen. Each of these elements can 
%          be made visible or invisible by toggling the property: 'visible'
%
% Author: Dante J Smith
% Update: 10/28/2021
curSess(strfind(curSess, '_')) = ' ';

monitorSize = get(0, 'Monitor');
numMon = size(monitorSize, 1);

if numMon == 2 && defMon == 2
    [~, mon] = max(monitorSize(:,1));
    
    figPosition = [monitorSize(mon,1) monitorSize(mon,2) monitorSize(mon,3) monitorSize(mon,4)]; 
else
    W = monitorSize(1, 3);
    H = monitorSize(1, 4);
    W2 = W/2;
    H2 = H/2;
    XPos = W2;
    YPos = 50;
    
    figPosition = [XPos YPos W2 H2];
end
msrStr.winPos = figPosition;

% Ready Annotation Dim
msrStr.rdAnoD = [700 300];
msrStr.rdAnoW = round(msrStr.rdAnoD(1)/msrStr.winPos(3), 2); 
msrStr.rdAnoH = round(msrStr.rdAnoD(2)/msrStr.winPos(4), 2);
msrStr.rdAnoX = 0.5 - msrStr.rdAnoW/2;
msrStr.rdAnoY = 0.5 - msrStr.rdAnoH/2;
msrStr.rdAnoPos = [msrStr.rdAnoX msrStr.rdAnoY msrStr.rdAnoW msrStr.rdAnoH];

% Cue Annotation Dim
msrStr.cuAnoD = [250 150];
msrStr.cuAnoW = round(msrStr.cuAnoD(1)/msrStr.winPos(3), 2); 
msrStr.cuAnoH = round(msrStr.cuAnoD(2)/msrStr.winPos(4), 2);
msrStr.cuAnoX = 0.5 - msrStr.cuAnoW/2;
msrStr.cuAnoY = 0.5 - msrStr.cuAnoH/2;
msrStr.cuAnoPos = [msrStr.cuAnoX msrStr.cuAnoY msrStr.cuAnoW msrStr.cuAnoH];

% EEE Annotation Dim
msrStr.eeAnoD = [370 170];
msrStr.eeAnoW = round(msrStr.eeAnoD(1)/msrStr.winPos(3), 2); 
msrStr.eeAnoH = round(msrStr.eeAnoD(2)/msrStr.winPos(4), 2);
msrStr.eeAnoX = 0.5 - msrStr.eeAnoW/2;
msrStr.eeAnoY = 0.5 - msrStr.eeAnoH/2;
msrStr.eeAnoPos = [msrStr.eeAnoX msrStr.eeAnoY msrStr.eeAnoW msrStr.eeAnoH];

%Assuming that perfect RMS lies at 50% of the screen, and the acceptable
%upper and lower bound are set at 45% and 55% of the screen, the edge of
%the screen will correspond to these variables accordingly. This will shift
%with differences in participants varying RMS values
msrStr.targRMS   = targRMS;              % Ideally this is between 60-80 dB
msrStr.hRMSrange = bounds*10;            % half the possible range of RMS, where each 5% mark represents 'bounds' dB
msrStr.vTopAmp   = msrStr.targRMS + msrStr.hRMSrange; %dB 
msrStr.vBotAmp   = msrStr.targRMS - msrStr.hRMSrange; %dB 

%Starts with some measurements
msrStr.bMar      = 0.05; %Bottom Margin
msrStr.minH      = 0.4;  %From the margin. This will need to change
msrStr.maxH      = 0.5;  %From the margin. This will need to change
msrStr.drawMinH  = msrStr.bMar + msrStr.minH; %45%
msrStr.drawMaxH  = msrStr.bMar + msrStr.maxH; %55% Goal value is 50%

msrStr.lMar      = 0.15; %Left Margin
msrStr.lineWidth = 0.1;  %Arbitrary-Looks good at the moment
msrStr.drawMinW  = msrStr.lMar + msrStr.lineWidth;

% Feedback Rectangle dim measurements 
msrStr.recWidth  = 0.05; % Arbitrary-Looks good at the moment
msrStr.recHeight = 0.05; % Start at this default
msrStr.recXSt    = msrStr.lMar + (msrStr.lineWidth/2 - msrStr.recWidth/2);
msrStr.recYSt    = msrStr.bMar; % Start at the margin

% Feedback Rectangle: Position and Size
msrStr.recPos = [msrStr.recXSt msrStr.recYSt msrStr.recWidth msrStr.recHeight];

% Feedback RMS Bounds Dim
msrStr.minLx = [msrStr.lMar msrStr.drawMinW];     %First X and Last X
msrStr.minLy = [msrStr.drawMinH msrStr.drawMinH]; %First Y and Last Y
msrStr.maxLx = [msrStr.lMar msrStr.drawMinW];     %First X and Last X
msrStr.maxLy = [msrStr.drawMaxH msrStr.drawMaxH]; %First Y and Last Y

% Trigger dim measurements
msrStr.r        = 60;                                % Trigger side L (pixels)
msrStr.visTrigX = 0;                                 % Trigger X Pos
msrStr.visTrigY = 0.02;                              % Trigger Y Pos
msrStr.visTrigW = round(msrStr.r/msrStr.winPos(3), 2); % Trigger X Len
msrStr.visTrigH = round(msrStr.r/msrStr.winPos(4), 2); % Trigger Y Len
msrStr.visTrigPos = [msrStr.visTrigX msrStr.visTrigY msrStr.visTrigW msrStr.visTrigH];

% Current Subject Note Dim
msrStr.subjND = [220 25];
msrStr.subjNX = msrStr.visTrigX;
msrStr.subjNY = msrStr.visTrigY + msrStr.visTrigH;
msrStr.subjNW = round(msrStr.subjND(1)/msrStr.winPos(3), 2); 
msrStr.subjNH = round(msrStr.subjND(2)/msrStr.winPos(4), 2); 
msrStr.subjNPos = [msrStr.subjNX msrStr.subjNY msrStr.subjNW msrStr.subjNH];

% Upcoming Trial 
msrStr.trialD = [400 150];
msrStr.trialW = round(msrStr.trialD(1)/msrStr.winPos(3), 2);
msrStr.trialH = round(msrStr.trialD(2)/msrStr.winPos(4), 2);
msrStr.trialX = 0.5 - msrStr.trialW/2;
msrStr.trialY = 0.8 - msrStr.trialH/2;
msrStr.trialPos = [msrStr.trialX msrStr.trialY msrStr.trialW msrStr.trialH];

% Upcoming Auditory Feedback 
msrStr.FBCueD = [600 150];
msrStr.FBCueW = round(msrStr.FBCueD(1)/msrStr.winPos(3), 2);
msrStr.FBCueH = round(msrStr.FBCueD(2)/msrStr.winPos(4), 2);
msrStr.FBCueX = 0.5 - msrStr.FBCueW/2;
msrStr.FBCueY = msrStr.trialY - msrStr.FBCueH;
msrStr.FBCuePos = [msrStr.FBCueX msrStr.FBCueY msrStr.FBCueW msrStr.FBCueH];

%%%%%%
VBFig = figure('NumberTitle', 'off', 'Color', [0 0 0], 'Position', msrStr.winPos, 'MenuBar', 'none');

%Plus Sign
annoStr.plus = annotation(VBFig,'textbox', msrStr.cuAnoPos,...
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
annoStr.EEE = annotation(VBFig,'textbox', msrStr.eeAnoPos,...
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
annoStr.Ready = annotation(VBFig,'textbox', msrStr.rdAnoPos,...
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
annoStr.visTrig = annotation(VBFig, 'rectangle', msrStr.visTrigPos,...
                              'Color',[1 1 1],...
                              'FaceColor', [1 1 1],...
                              'visible', 'off');

% Visual feedback bar                    
annoStr.LoudRec = annotation(VBFig, 'rectangle', msrStr.recPos,...
                              'Color',[0 1 0],...
                              'LineStyle','none',...
                              'FaceColor',[0 1 0],...
                              'visible','off');

% Lower limit of acceptable RMS                     
LoudminLine = annotation(VBFig, 'line', msrStr.minLx, msrStr.minLy,...
                                  'Color',[1 1 1],...
                                  'LineStyle', '--',...
                                  'LineWidth', 1, ...
                                  'visible','off');

% Upper limit of acceptable RMS                              
LoudmaxLine = annotation(VBFig, 'line', msrStr.maxLx, msrStr.maxLy,...
                                  'Color',[1 1 1],...
                                  'LineStyle', '--',...
                                  'LineWidth', 1, ...
                                  'visible','off');
% Combine the limit handles
annoStr.fbLines = [LoudminLine LoudmaxLine]; 
                              
% Subject Note
annoStr.curSessNote = annotation(VBFig,'textbox', msrStr.subjNPos,...
                                'Color',[1 1 1],...
                                'String', curSess,...
                                'FitBoxToText','on',...
                                'LineStyle','none',...
                                'VerticalAlignment','bottom',...
                                'HorizontalAlignment','left',...
                                'FontSize', 10,...
                                'FontName','Arial',...
                                'EdgeColor','none',...
                                'BackgroundColor',[0 0 0],...
                                'Visible','on');
                               
% Trial Title
annoStr.trialT = annotation(VBFig, 'textbox', msrStr.trialPos,...
                           'Color', [1 1 1],...
                           'String', 'Trial 1',...
                           'FitBoxToText','on',...
                           'LineStyle','none',...
                           'VerticalAlignment','bottom',...
                           'HorizontalAlignment','center',...
                           'FontSize', 130,...
                           'FontName','Arial',...
                           'EdgeColor','none',...
                           'BackgroundColor',[0 0 0],...
                           'Visible','off');
                       
% Feedback Cue
annoStr.FBCue = annotation(VBFig, 'textbox', msrStr.FBCuePos,...
                           'Color', [1 1 1],...
                           'String', 'Trial 10',...
                           'FitBoxToText','on',...
                           'LineStyle','none',...
                           'VerticalAlignment','bottom',...
                           'HorizontalAlignment','center',...
                           'FontSize', 100,...
                           'FontName','Arial',...
                           'EdgeColor','none',...
                           'BackgroundColor',[0 0 0],...
                           'Visible','off');
end
