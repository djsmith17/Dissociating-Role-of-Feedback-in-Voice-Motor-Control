classdef LiveSensorResult
   properties
       curSess
       numTrial
       fs
       trigs
       balloon
       defMon
       trialColors
       tag
       curve
       resH
       resA
   end
   methods
       function obj = LiveSensorResult(expParam, defMon)
           if nargin > 0
               
               curSess  = expParam.curSess;
               curSess(strfind(curSess, '_')) = ' ';
               numTrial = expParam.numTrial;
               fs       = expParam.sRateQ;
               trigs    = expParam.trigs(:,:,2);
               balloon  = expParam.balloon;
               balloon(strfind(balloon, '_')) = '';
               
               obj.curSess  = curSess;
               obj.numTrial = numTrial;
               obj.fs       = fs;
               obj.trigs    = trigs;
               obj.balloon  = balloon;
               
               obj.defMon      = defMon;
               obj.trialColors = distinguishable_colors(numTrial);
               
               obj.tag   = {};
               obj.curve = [];
               
               [obj.resH, obj.resA] = initFigure(obj);
           end       
       end
       function [resH, resA] = initFigure(obj)
           % resH = initFigure(obj) creates a figure by which to display 
           % recorded signal during a experimental session. This is 
           % currently configured to display the pressure in the 
           % perturbatron balloon.
           %
           % resH    : Figure handle for the generated figure.
           % resA    : Axes handle for the generated plots
           monitorSize = get(0, 'Monitor');
           numMon = size(monitorSize, 1);
           if numMon == 2 && obj.defMon == 2
               plotDim = [800 600];
               [~, mon] = max(monitorSize(:,1));
               halfW  = monitorSize(mon, 3)/2;
               halfWD = halfW - plotDim(1)/2 + monitorSize(mon, 1) - 1;
               
               figPosition = [halfWD 80 plotDim];
           else
               plotDim = [800 600];
               monW = monitorSize(1, 3);
               halfWD = monW - plotDim(1) - 10;
               
               figPosition = [halfWD 50 plotDim];
           end
           winPos = figPosition;
           
           resH = figure('Name', 'Pressure Live Result', 'Color', [1 1 1], 'Position', winPos);
           resA = axes();
           plot(resA, [1 1], [-1 5], 'k-', 'LineWidth', 2);
           axis([0 3.5 -0.5 5.0])
           box off
           set(gca,'FontSize', 12,...
                   'XTickLabel', {'-1.0' '-0.5' '0' '0.5' '1.0' '1.5' '2.0' '2.5'},...
                   'FontWeight', 'bold')
           xlabel('Time (s)', 'FontSize', 14, 'FontWeight', 'bold')
           ylabel('Voltage (V)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k')
           title(obj.curSess)
              
           annotation('textbox', [0.75 0.8 0.45 0.1],...
                      'string', {['Balloon: ' obj.balloon ]},...
                      'LineStyle','none',...
                      'FontWeight','bold',...
                      'FontSize',10,...
                      'FontName','Arial');
           hold on
       end
       
       function obj = updateLiveResult(obj, daqIn, curTrial)
           sig      = daqIn(:,4);

           St = obj.trigs(curTrial,1) - obj.fs*1 + 1;
           Sp = obj.trigs(curTrial,1) + obj.fs*2.5;
           sigSnip = sig(St:Sp);
           time    = (0:1/obj.fs :(length(sigSnip)-1)/obj.fs)';
           trialNmList = ['Trial ' num2str(curTrial)];
           
           trPrs = plot(obj.resA, time, sigSnip, 'LineWidth', 2, 'Color', obj.trialColors(curTrial, :));

           obj.tag   = cat(1, obj.tag, trialNmList);
           obj.curve = cat(1, obj.curve, trPrs);

           lgd = legend(obj.resA, obj.curve, obj.tag);
           set(lgd, 'box', 'off',...
                    'location', 'NorthWest'); 
       end
   end
end