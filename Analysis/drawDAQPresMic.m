function drawDAQPresMic(res, saveResultsDir)
%Plots multiple signals against each other. Creates a new plot for each 
%trial. Displays relevant information

curSess  = res.curSess;  % The current experiment details (Subject/Run)
f0b      = res.f0b;
numPT    = res.numPertTrials;
trialIdx = res.pertIdx;
trigs    = res.pertTrig; % Where the perturbations occur

timeS        = res.timeS;
sensorP      = res.sensorPsv;
lagTimeP     = res.lagTimePsv;
riseTimeP    = res.riseTimePsv;
OnOfValP     = res.OnOfValPsv;
limitsP      = res.limitsP;

timef0       = res.timef0;
audioM       = res.audioMf0TrialPert;
limitsf0     = res.limitsA;

curSess(strfind(curSess, '_')) = ' ';

plotpos = [500 300];
plotdim = [800 600];
pertColor = [0.8 0.8 0.8];
micColor  = [0.1 0.5 0.1];

for ii = 1:numPT
    DAQPresMic(ii) = figure('Color', [1 1 1]);
    set(DAQPresMic(ii), 'Position',[plotpos plotdim],'PaperPositionMode','auto') 
    
    pertAx  = [trigs(ii,1), trigs(ii,2)];
    pertAy  = [200 200];
    
    pA = area(pertAx, pertAy, -200, 'FaceColor', pertColor, 'EdgeColor', pertColor);
    
    hold on
    [aX, mS, pS] = plotyy(timef0, audioM(:,ii), timeS, sensorP(:,ii));

    set(mS, 'LineWidth', 2,... 
            'Color', micColor);    

    set(pS, 'LineWidth', 2,... 
            'Color', 'k');

    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold') 
    ylabel(aX(1), 'f0 (Cents)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', micColor)
    ylabel(aX(2), 'Pressure (psi)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k') 
    title({'Pressure and f0 changes due to Balloon Inflation';
            [curSess ' Trial ' num2str(trialIdx(ii))]}, 'FontSize', 12, 'FontWeight', 'bold')
    axis(aX(1), limitsf0); 
    axis(aX(2), limitsP);
    box off

    set(aX, {'ycolor'}, {micColor;'k'},...
             'FontSize', 12,...
             'FontWeight', 'bold')

    pltlgd = legend([pA mS pS], 'Perturbation Period', 'Microphone f0', 'Pressure Sensor');
    set(pltlgd, 'box', 'off',...
            'location', 'NorthWest'); 

    t = annotation('textbox',[0.60 0.82 0.40 0.1],...
                   'string', {['Pressure Onset Lag: ' num2str(1000*lagTimeP(ii,1)) 'ms'];...
                              ['Pressure Rise Time: ' num2str(1000*riseTimeP(ii,1)) 'ms'];...
                              ['Onset/Offset Pressure: ' num2str(OnOfValP(ii,1), '%1.2f') 'psi, ' num2str(OnOfValP(ii,2), '%1.2f') 'psi']},...
                    'LineStyle','none',...
                    'FontWeight','bold',...
                    'FontSize',8,...
                    'FontName','Arial');
    plotTitle =  '_DAQMicPres ';


    plTitle = [curSess  plotTitle num2str(ii) '.jpg'];     
    saveFileName = fullfile(saveResultsDir, plTitle);
    export_fig(saveFileName) 

end
end