function drawSubjRespVarVoicingCorr(dirs, pA, pooledRunStr)

numSubj = length(pooledRunStr);

fontN = 'Arial';
axisLSize = 14;

plotpos = [10 40];
plotdim = [1800 900];
IndivTrialAudResp = figure('Color', [1 1 1]);
set(IndivTrialAudResp, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

ha = tight_subplot(3, 6, [0.08 0.03],[0.06 0.15],[0.05 0.03]);

pA.pAnalysis(strfind(pA.pAnalysis, '_')) = '';
scattColors = {'b', 'r'};

for ii = 1:numSubj
    axes(ha(ii))
     
    curRes = pooledRunStr(ii);
    curSubj  = curRes.subject;
    
    allRespPer = [];
    allTimes   = [];
    allStatSentence = {};
    for jj = 1:pA.numCond
        respPer      = curRes.respVarSingle{jj}.respPer;
        voicingTimes = curRes.prePertVoicingTimePert{jj};

%         wayWrong = respPer > 300 | respPer < -300;
%         respPer      = respPer(~wayWrong);
%         voicingTimes = voicingTimes(~wayWrong);
%         
        numTrial = length(respPer);

        varToCorr = [voicingTimes, respPer];
        [corrR, corrP] = corrcoef(varToCorr);

        plot(voicingTimes, respPer, 'o', 'color', scattColors{jj}, 'MarkerSize', 5)
        hold on
        
        allRespPer = [allRespPer; respPer];
        allTimes   = [allTimes; voicingTimes];
        
        statSentence = sprintf('r = %0.2f, p = %0.3f, n = %d\n', corrR(1,2), corrP(1,2), numTrial);
        allStatSentence = cat(1, allStatSentence, statSentence);
    end
    respPerMin = min(allRespPer);
    respPerMax = max(allRespPer);
    lowB = respPerMin - 10;
    upB  = respPerMax + 50;
    
    timesMin = min(allTimes);
    timesMax = max(allTimes);
    leftB = timesMin - 0.1;
    rightB = timesMax + 0.1;
    
    box off
    axis([leftB rightB lowB upB])

    if ii == 13
        xlabel('Pre Pert Voicing Time (s)')
        ylabel('RespPer (%)')
    end
     
    set(gca, 'FontName', fontN,...
             'FontSize', axisLSize,...
             'FontWeight','bold')

    title(curSubj)
    
    legend(gca, allStatSentence,...
        'EdgeColor', 'none',...
        'FontSize', 8,...
        'Location', 'northwest')
    
end

% legend(IndivTrialAudResp, pA.pubCond,...
%     'Position', [0.78 0.93, 0.1 0.05],...
%     'Orientation', 'Horizontal',...
%     'EdgeColor', 'none');

suptitle(pA.pAnalysis)

pltTitle = [pA.pAnalysis 'IndividualParticipantVoicingRespPerCorr.jpg'];
plotFileName = fullfile(dirs.SavResultsDir, pltTitle);
export_fig(plotFileName)
end