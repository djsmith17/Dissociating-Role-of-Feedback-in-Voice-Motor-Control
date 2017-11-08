function drawMaskvVoiceMeanf0(niResM, niResV, statLib, pltName, plotFolder)

curSess          = niResM.subject;
f0b              = round(niResM.f0b);
numMasked        = niResM.numPertTrials;
numVoiced        = niResV.numPertTrials;

time              = niResM.secTime;
meanf0PertOnsetM  = niResM.audioMf0MeanPert(:,1);
CIf0PertOnsetM    = niResM.audioMf0MeanPert(:,2);
meanf0PertOffsetM = niResM.audioMf0MeanPert(:,3);
CIf0PertOffsetM   = niResM.audioMf0MeanPert(:,4);

meanf0PertOnsetV  = niResV.audioMf0MeanPert(:,1);
CIf0PertOnsetV    = niResV.audioMf0MeanPert(:,2);
meanf0PertOffsetV = niResV.audioMf0MeanPert(:,3);
CIf0PertOffsetV   = niResV.audioMf0MeanPert(:,4);
limits            = niResM.limitsAmean;

statSMM = round(10*statLib(1))/10;
statSMV = round(10*statLib(2))/10;
statRMM = round(10*statLib(3))/10;
statRMV = round(10*statLib(4))/10;
statRPM = round(statLib(5));
statRPV = round(statLib(6));
statSP  = statLib(7);
statRP  = statLib(8);
statPP  = statLib(9);

pValueThresh = 0.05;

plotpos = [10 100];
plotdim = [1600 600];
InterTrialf0 = figure('Color', [1 1 1]);
set(InterTrialf0, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

curSess(strfind(curSess, '_')) = ' ';

dottedStartx = [0 0];
dottedy      = [-300 300];

ha = tight_subplot(1,2,[0.1 0.05],[0.12 0.15],[0.05 0.05]);

%Onset of Perturbation
axes(ha(1))
shadedErrorBar(time, meanf0PertOnsetM, CIf0PertOnsetM, 'b', 1); %Masked
hold on
shadedErrorBar(time, meanf0PertOnsetV, CIf0PertOnsetV, 'r', 1); %Voice
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')

title('Onset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
axis(limits); box off

set(gca,'XTickLabel',{'-0.5' '0' '0.5' '1.0'},...
        'FontSize', 16,...
        'FontWeight','bold')

%Offset of Perturbation
axes(ha(2))
uH = shadedErrorBar(time, meanf0PertOffsetM, CIf0PertOffsetM, 'b', 1); %Masked
hold on
pH = shadedErrorBar(time, meanf0PertOffsetV, CIf0PertOffsetV, 'r', 1); %Voice
hold on
plot(dottedStartx, dottedy,'k','LineWidth',4)
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); ylabel('f0 (cents)', 'FontSize', 18, 'FontWeight', 'bold')

title('Offset of Perturbation', 'FontSize', 18, 'FontWeight', 'bold')
axis(limits); box off
set(gca,'XTickLabel', {'-0.5' '0' '0.5' '1.0'},...
        'FontSize', 16,...
        'FontWeight','bold',...
        'YAxisLocation', 'right');

sup = suptitle(curSess);
set(sup, 'FontSize', 20,...
         'FontWeight','bold')
     
annoStim = ['SM (M/NM): ' num2str(statSMM) ' cents / ' num2str(statSMV) ' cents'];
annoResp = ['RM (M/NM): ' num2str(statRMM) ' cents / ' num2str(statRMV) ' cents'];
annoPerc = ['RP (M/NM): ' num2str(statRPM) '% / ' num2str(statRPV) '%'];

annoStim = checkSig(statSP, pValueThresh, annoStim);
annoResp = checkSig(statRP, pValueThresh, annoResp);
annoPerc = checkSig(statPP, pValueThresh, annoPerc);
 
t = annotation('textbox',[.25 .75 0.45 0.1],...
               'string', {annoStim;
                          annoResp
                          annoPerc},...
                'LineStyle','none',...
                'FontWeight','bold',...
                'FontSize',12,...
                'FontName','Arial');

legend([uH.mainLine pH.mainLine],{[num2str(numMasked) ' M Trials'], [num2str(numVoiced) ' NM Trials']},...
            'box', 'off',...
            'FontSize', 14,...
            'FontWeight', 'bold',...
            'Position', [0.80 0.75 0.1 0.1]);

plots = {'Figure'};
for i = 1:length(plots)
    plTitle = [pltName '.jpg'];

    saveFileName = fullfile(plotFolder, plTitle);
    export_fig(saveFileName)
end
end

function anno = checkSig(stat, thresh, anno)

if stat < thresh
    anno = [anno '*'];
end
end