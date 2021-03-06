function drawDAQsignal_Long(niAn, plotFlag, saveResultsDir, sv2F)
%Plots multiple signals against each other. Includes all trials in the 
%plot. Displays relevant information

curExp   = niAn.curExp;   %The current experiment detials (Subject/Run)
numTrial = niAn.numTrial; %Number of Trials
trigs    = niAn.pertTrig; %Where the perturbations occur

time     = niAn.time_DN;

sensorP        = niAn.sensorP_DN;
lagTimeP       = niAn.lagsPres;
riseTimeP      = niAn.riseTimeP;
rangePressures = niAn.rangePressures;
pLimits        = niAn.pLimits;

sensorFC   = niAn.sensorFC_DN;
sensorFN   = niAn.sensorFN_DN;
lagTimeFC  = niAn.lagsFC;
lagTimeFN  = niAn.lagsFN;
fLimits    = niAn.fLimits;

timef0   = niAn.timef0;
micf0    = niAn.audioMf0_norm;
headf0   = niAn.audioHf0_norm;
aLimits  = niAn.aLimits;

curExp(strfind(curExp, '_')) = ' ';

plotpos = [100 100];
plotdim = [1600 600];
pertColor = [0.8 0.8 0.8];

for ii = 1:numTrial
    DAQOutput(ii) = figure('Color', [1 1 1]);
    set(DAQOutput(ii), 'Position',[plotpos plotdim],'PaperPositionMode','auto')
    
    pertAx  = [trigs(ii,1), trigs(ii,2)];
    pertAy  = [200 200];
    
    pA = area(pertAx, pertAy, -200, 'FaceColor', pertColor, 'EdgeColor', pertColor);
    hold on

    if plotFlag == 1
        plot(time, sensorP(:,ii), 'k', 'LineWidth', 2);   
      
        xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold') 
        ylabel('Pressure (psi)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k') 
        title({'Sensors Measurements due to Balloon Inflation';
                curExp}, 'FontSize', 12, 'FontWeight', 'bold')
 
        axis(pLimits);
        box off
    
        set(gca, 'LineWidth', 2,...
                 'FontSize', 12,...
                 'FontWeight', 'bold')
        
        pltlgd = legend('Perturbation Period', 'Pressure Sensor');
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
        plotTitle =  '_DAQSignalOutput ';
    else
        aM = plot(timef0, micf0(:,ii), 'r');
        hold on
        aH = plot(timef0, headf0(:,ii), 'b');
        
        xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold')
        ylabel('f0 (Cents)', 'FontSize', 18, 'FontWeight', 'bold')
        title({'Audio Measurements due to Audapter Pertrubations';
                curExp}, 'FontSize', 12, 'FontWeight', 'bold')
        axis(aLimits);
        box off
        
        set(gca, 'FontSize', 12,...
                 'FontWeight', 'bold');
        
        pltlgd = legend([pA aM aH], 'Perturbation Period', 'Microphone', 'Headphones');
        set(pltlgd, 'box', 'off',...
                'location', 'NorthWest'); 
            
        plotTitle =  '_DAQAudioInput ';
    end
   
    if sv2F == 1
        plTitle = [curExp  plotTitle num2str(ii) '.jpg'];     
        saveFileName = fullfile(saveResultsDir, plTitle);
        export_fig(saveFileName) 
    end
end
end