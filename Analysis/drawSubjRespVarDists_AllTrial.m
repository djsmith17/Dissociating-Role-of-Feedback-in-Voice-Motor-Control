function drawSubjRespVarDists_AllTrial(dirs, pooledRunStr)

numSubj = length(pooledRunStr);

fontN = 'Times New Roman';
axisLSize = 14;

plotpos = [10 10];
plotdim = [1600 800];
IndivTrialAudResp = figure('Color', [1 1 1]);
set(IndivTrialAudResp, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

allTrials      = 0;
allTrialsLess0 = 0;
allRespPer     = [];
for ii = 1:numSubj
     
     curRes = pooledRunStr(ii);
     respPer  = curRes.respVarSingle{1,1}.respPer;
     numTrial = length(respPer);
     
     trialsUnder0 = respPer < 0;
     numTrialsLess0 = sum(trialsUnder0);
     
     allTrials = allTrials + numTrial;
     allTrialsLess0 = allTrialsLess0 + numTrialsLess0;
     allRespPer = cat(1, allRespPer, respPer);
end

histogram(allRespPer, 30)
hold on
plot([0 0], [-100 100], 'k--')

respPerMin = min(allRespPer);
respPerMax = max(allRespPer);
lB = respPerMin - 20;
rB = respPerMax + 20;

xlabel('Response Percentage (%)')
box off; axis([lB rB 0 45])
set(gca, 'FontName', fontN,...
         'FontSize', axisLSize,...
         'FontWeight','bold')

AllPerLess0   = round(100*allTrialsLess0/allTrials,2);
title({'Auditory Feedback Perturbation Experiment', ['Total Following Response Trials: ' num2str(AllPerLess0) '%']})

pltTitle = 'InvestRespAudPertExp_AllTrial.jpg';
plotFileName = fullfile(dirs.SavResultsDir, pltTitle);
export_fig(plotFileName)
end