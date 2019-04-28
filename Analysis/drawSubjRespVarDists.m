function drawSubjRespVarDists(dirs, pooledRunStr)

numSubj = length(pooledRunStr);

fontN = 'Arial';
axisLSize = 14;

plotpos = [10 10];
plotdim = [1600 800];
IndivTrialAudResp = figure('Color', [1 1 1]);
set(IndivTrialAudResp, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

ha = tight_subplot(3, 6, [0.08 0.03],[0.12 0.15],[0.05 0.03]);

allTrials      = 0;
allTrialsLess0 = 0;
for ii = 1:numSubj
     axes(ha(ii))
     
     curRes = pooledRunStr(ii);
     respPer  = curRes.respVarSingle{1,1}.respPer;
     respPerM = round(curRes.respVarSingle{1,1}.respPerM, 2);
     numTrial = length(respPer);
     curSubj  = curRes.subject;
     
     trialsUnder0 = respPer < 0;
     numTrialsLess0 = sum(trialsUnder0);
     perLess0   = round(100*numTrialsLess0/numTrial,2);
     
     allTrials = allTrials + numTrial;
     allTrialsLess0 = allTrialsLess0 + numTrialsLess0;
     
     
     respPerMin = min(respPer);
     respPerMax = max(respPer);
     lB = respPerMin - 20;
     rB = respPerMax + 20;
     
     histogram(curRes.respVarSingle{1,1}.respPer, 10)
     hold on
     plot([0 0], [-100 100], 'k--')   
     box off
     
     if ii == 13
         xlabel('Response Percentage (%)')
     end
     
     set(gca, 'FontName', fontN,...
             'FontSize', axisLSize,...
             'FontWeight','bold')
     
     if perLess0 <= 0
         titColor = 'k';
     else
         titColor = 'r';
     end
     
%      title({curSubj,['Mean: ' num2str(respPerM) '%, (n=' num2str(numTrial) ')']})
    title([curSubj ': ' num2str(perLess0) '%'], 'color', titColor)
     
     axis([lB rB 0 7])
end

AllPerLess0   = round(100*allTrialsLess0/allTrials,2);
suptitle({'Auditory Feedback Perturbation Experiment', ['Total Following Response Trials: ' num2str(AllPerLess0) '%']})

pltTitle = 'InvestRespAudPertExp.png';
plotFileName = fullfile(dirs.SavResultsDir, pltTitle);
export_fig(plotFileName)
end