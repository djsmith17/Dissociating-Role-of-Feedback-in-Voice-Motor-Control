function drawDAQsignal(niAn, time, fSensorC, fSensorN, pSensor, saveResultsDir, sv2F)
%Plots multiple signals against each other. Creates a new plot for each 
%trial. Displays relevant information

curExp   = niAn.curExp;   %The current experiment detials (Subject/Run)
numTrial = niAn.numTrial; %Number of Trials
trigs    = niAn.pertTrig; %Where the perturbations occur

lagTimeP       = niAn.lagsPres;
riseTimeP      = niAn.riseTimeP;
rangePressures = niAn.rangePressures;
pLimits        = niAn.pLimits;

lagTimeFC  = niAn.lagsFC;
lagTimeFN  = niAn.lagsFN;
fLimits    = niAn.fLimits;

curExp(strfind(curExp, '_')) = ' ';

plotpos = [500 300];
plotdim = [800 600];
pertColor = [0.8 0.8 0.8];

for ii = 1:numTrial
    DAQOutput(ii) = figure('Color', [1 1 1]);
    set(DAQOutput(ii), 'Position',[plotpos plotdim],'PaperPositionMode','auto')
    
    pertAx  = [trigs(ii,1), trigs(ii,2)];
    pertAy  = [6 6];
    
    pA = area(pertAx, pertAy, -1, 'FaceColor', pertColor, 'EdgeColor', pertColor);
    hold on

    [aX, fS, pS] = plotyy([time time], [fSensorC(:,ii) fSensorN(:,ii)], time, pSensor(:,ii));
    
    set(pS, 'LineWidth', 2,... 
            'Color', 'k');
    
    set(fS(1), 'LineWidth', 2,... 
               'Color', 'b');      
    set(fS(2), 'LineWidth', 2,... 
               'Color', 'r');
    
    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold') 
    ylabel(aX(1), 'Voltage (V)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k')
    ylabel(aX(2), 'Pressure (psi)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k') 
    title({'Sensors Measurements due to Balloon Inflation';
           curExp}, 'FontSize', 12, 'FontWeight', 'bold')
    axis(aX(1), fLimits); 
    axis(aX(2), pLimits);
    box off
    
    set(aX, {'ycolor'}, {'k';'k'},...
            'FontSize', 12,...
            'FontWeight', 'bold')
        
    pltlgd = legend([pA pS fS(1) fS(2)], 'Perturbation Period', 'Pressure Sensor', 'Force Sensor: Collar', 'Force Sensor: Neck');
    set(pltlgd, 'box', 'off',...
                'location', 'NorthWest'); 
            
    t = annotation('textbox',[0.60 0.82 0.40 0.1],...
                   'string', {['Pressure Onset Lag: ' num2str(1000*lagTimeP(ii,1)) 'ms'];...
                              ['Pressure Rise Time: ' num2str(1000*riseTimeP(ii,1)) 'ms'];...
                              ['Collar Force Onset Lag: ' num2str(1000*lagTimeFC(ii,1)) 'ms'];...
                              ['Neck Force Onset Lag: ' num2str(1000*lagTimeFN(ii,1)) 'ms'];...
                              ['Onset/Offset Pressure Values: ' num2str(rangePressures(ii,1), '%1.2f') 'psi, ' num2str(rangePressures(ii,2), '%1.2f') 'psi']},...
                    'LineStyle','none',...
                    'FontWeight','bold',...
                    'FontSize',8,...
                    'FontName','Arial');
   
    if sv2F == 1
        plTitle = [curExp  '_DAQSignalOutput ' num2str(ii) '.png'];     
        saveFileName = fullfile(saveResultsDir, plTitle);
        export_fig(saveFileName) 
    end
end
end