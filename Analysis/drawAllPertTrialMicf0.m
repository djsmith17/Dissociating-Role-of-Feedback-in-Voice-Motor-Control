function drawAllPertTrialMicf0(res, plotFolder, varargin)
% drawAllPertTrialMicf0(res, plotFolder) plots the pitch traces 

if isempty(varargin)
    presFlag = 0;
else
    presFlag = varargin{1};
end

curSess          = res.curSess;
f0b              = round(res.f0b, 1); % Baseline f0 rounded to 0.1 Hz
AudFB            = res.AudFB;
numPT            = res.numPertTrialsFin;
pertTrig         = res.pertTrigsFin;

time             = res.timef0;
sigs             = res.audioMf0TrialPert;
limits           = res.limitsA;

timePres         = res.timeS;
sensorP          = res.sensorPsv;
pressureLim      = res.limitsP;

curSess(strfind(curSess, '_')) = ' ';

plotpos = [10 100];
plotdim = [1600 800];
InterTrialf0 = figure('Color', [1 1 1]);
set(InterTrialf0, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

pertColor = [0.8 0.8 0.8];

ha = tight_subplot(2,5,[0.15 0.08],[0.12 0.15],[0.05 0.05]);

for ii = 1:numPT
    axes(ha(ii))
    
    pertAx  = [pertTrig(ii,1), pertTrig(ii,2)];
    pertAy  = [600 600];
    
    pA = area(pertAx, pertAy, -600, 'FaceColor', pertColor, 'EdgeColor', pertColor);
    hold on    
    
    if presFlag == 1
        yyaxis right
        plot(timePres, sensorP(:,ii), '--k', 'LineWidth', 1.5)
        
        ylabel('Pressure (psi)')
        axis(pressureLim);
        set(gca,'FontSize', 14,...
                'FontWeight','bold')
        yyaxis left
    end  
       
    plot(time, sigs(:,ii), 'b', 'LineWidth', 2)
    xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold')
    ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')
    title(['Perturbed Trial ' num2str(ii)], 'FontSize', 18, 'FontWeight', 'bold')
    axis(limits); box off
    
    set(gca,'FontSize', 14,...
            'FontWeight','bold')
    
 
end

sup = suptitle({curSess; [ 'AudFB: ' AudFB]; ['f0: ' num2str(f0b) ' Hz']});
set(sup, 'FontSize', 18, 'FontWeight', 'bold')

plots = {'IntraTrialf0'};
for i = 1:length(plots)
    plTitle = [curSess '_' plots{i} '.jpg'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end
end