function drawExpPressureDist(dirs, pA, pooledRunStr)

numSubj = length(pooledRunStr);

fontN = 'Arial';
axisLSize = 14;

mu = '\mu';
sigma  = '\sigma'; 

plotpos = [10 10];
plotdim = [1600 800];
IndivTrialAudResp = figure('Color', [1 1 1]);
set(IndivTrialAudResp, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

ha = tight_subplot(3, 6, [0.13 0.03],[0.12 0.18],[0.05 0.03]);

allMeanPress = [];
for ii = 1:numSubj
     axes(ha(ii))
     
     curRes  = pooledRunStr(ii);
     curSubj = curRes.studyID;
     curSubj(strfind(curSubj, '_')) = ' ';
     
     pressureOn   = curRes.sensorPOnOff(:,1);
     numTrial     = length(pressureOn);
     pressureOnM  = round(mean(pressureOn), 2);
     pressureOnSD = round(std(pressureOn), 2);
     
     pressureOnSE = pressureOnSD/sqrt(numTrial);
     
     respPerMin = min(pressureOn);
     respPerMax = max(pressureOn);
     lB = respPerMin - 0.02;
     rB = respPerMax + 0.02;
     
     allMeanPress = cat(1, allMeanPress, pressureOnM);
     
     histogram(pressureOn, 10)
     box off
     
     if ii == 13
         xlabel('Pressure (psi)')
     end
     
     set(gca, 'FontName', fontN,...
              'FontSize', axisLSize,...
              'FontWeight','bold')
     
     
     title({[curSubj ' (n=' num2str(numTrial) ')'],[mu '=' num2str(pressureOnM) ', ' sigma '=' num2str(pressureOnSD)]})

     axis([lB rB 0 15])
end

bigTit = suptitle({'Masking Study', 'Distributions of pressure readings at the appex of the pressure step function', [mu '= Mean, ' sigma '= Standard Deviation']});

set(bigTit, 'FontName', fontN,...
              'FontSize', axisLSize,...
              'FontWeight','bold')

pltTitle = 'ExpPressureDistAllSubj.jpg';
plotFileName = fullfile(dirs.SavResultsDir, pltTitle);
export_fig(plotFileName)

% Perform Standard Summary Stats
measureVar.varName = 'MeanPressureLevels';
measureVar.condition = 'AllConditions';
measureVar.units = 'psi';
summaryStat = MeasureSummaryStats(dirs, pA, measureVar, allMeanPress, 0);

summaryStat = summaryStat.testNormality(); % Describe the normality  
summaryStat.drawHistoBoxCombo()            % Visualize Normality/Outliers
end