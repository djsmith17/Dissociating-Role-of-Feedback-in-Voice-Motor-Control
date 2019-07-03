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

fullRespPer = cell(pA.numCond, 1);
fullVTimes  = cell(pA.numCond, 1);
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

        wayWrong = respPer > 300 | respPer < -300;
        respPer      = respPer(~wayWrong);
        voicingTimes = voicingTimes(~wayWrong);

        respPerNorm = NormalizeByMax(respPer);
    
        numTrial = length(respPer);

        varToCorr = [voicingTimes, respPerNorm];
        [corrR, corrP] = corrcoef(varToCorr);

        plot(voicingTimes, respPerNorm, 'o', 'color', scattColors{jj}, 'MarkerSize', 5)
        hold on
        
        allRespPer = [allRespPer; respPerNorm];
        allTimes   = [allTimes; voicingTimes];
        
        fullRespPer{jj} = [fullRespPer{jj}; respPerNorm];
        fullVTimes{jj}  = [fullVTimes{jj}; voicingTimes];
        
        statSentence = sprintf('r = %0.2f, p = %0.3f, n = %d\n', corrR(1,2), corrP(1,2), numTrial);
        allStatSentence = cat(1, allStatSentence, statSentence);
    end
    respPerMin = min(allRespPer);
    respPerMax = max(allRespPer);
    lowB = respPerMin - 0.25;
    upB  = respPerMax + 1;
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%

plotpos = [10 40];
plotdim = [1800 900];
AllTrialVTimingCorr = figure('Color', [1 1 1]);
set(AllTrialVTimingCorr, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

allRespPerBounds = [];
allVTimesBounds = [];
lgdNames = {};
for jj = 1:pA.numCond
    respPer = fullRespPer{jj};
    vTimes  = fullVTimes{jj};
    
    allRespPerBounds = [allRespPerBounds; max(respPer)];
    allRespPerBounds = [allRespPerBounds; min(respPer)];
        
    allVTimesBounds = [allVTimesBounds max(vTimes)];
    allVTimesBounds = [allVTimesBounds min(vTimes)];
    
    numTrial = length(respPer);
    
    varToCorr = [vTimes, respPer];
    [corrR, corrP] = corrcoef(varToCorr);
    
    plot(vTimes, respPer, 'o', 'color', scattColors{jj}, 'MarkerSize', 5)
    hold on
    
    statSentence = sprintf('r = %0.2f, p = %0.3f, n = %d\n', corrR(1,2), corrP(1,2), numTrial);
    lgdNames = cat(1, lgdNames, statSentence);
end
respPerMin = min(allRespPerBounds);
respPerMax = max(allRespPerBounds);
lowB = respPerMin - 0.25;
upB  = respPerMax + 0.25;

timesMin = min(allVTimesBounds);
timesMax = max(allVTimesBounds);
leftB = timesMin - 0.1;
rightB = timesMax + 0.1;

xlabel('Pre Pert Voicing Time (s)')
ylabel('Normalized RespPer (%)')
title(pA.pAnalysis)

box off
axis([leftB rightB lowB upB])

set(gca, 'FontName', fontN,...
         'FontSize', axisLSize,...
         'FontWeight','bold')


legend(lgdNames,...
    'EdgeColor', 'none',...
    'FontSize', 15,...
    'Location', 'northwest')

pltTitle = [pA.pAnalysis 'AllParticipantVoicingRespPerCorr.jpg'];
plotFileName = fullfile(dirs.SavResultsDir, pltTitle);
export_fig(plotFileName)
end

function NormVec = NormalizeByMax(vector)

maxVal = max(vector);

NormVec = vector/maxVal;

end