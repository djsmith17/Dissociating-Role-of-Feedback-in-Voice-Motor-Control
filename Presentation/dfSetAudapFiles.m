function audStimP = dfSetAudapFiles(dirs, expParam, trial)
% audStimP = dfSetAudapFiles(expParam, dirs, trial, debug) sets the ost and
% pcf functions for a custom pitch-shift reflex experiment recorded in 
% Audapter. 
% This script creates two possible pitch-shifts. The first is a linear ramp
% stimulus that has a rise/fall time of 150ms, and magnitude of 100 cents.
% The second is sigmoidal ramp stimulus that matches the rise/fall time 
% and magnitude of the laryngeal perturbation experiment.

trialType = expParam.trialType(trial);
trigs     = expParam.trigs(trial,:,1);

pertName  = expParam.AudPert;
pertSw    = expParam.AudPertSw;

InflaT = expParam.InflaT;
InflaV = expParam.InflaV;

audStimP = organizeStimulus(trialType, trigs, pertSw, pertName, InflaT, InflaV);

if trialType
%     drawStimulus(audStimP, dirs)
end
end

function audStimP = organizeStimulus(trialType, trigs, pertSw, pertName, InflaT, InflaV)

audStimP.trialType = trialType; % 0: Control 1: Catch
audStimP.pertSw    = pertSw;    % 0: -100c   1: LarMag
audStimP.pertName  = pertName;  % 'Linear Standard', 'Sigmoid Matched'
audStimP.tStep     = 0.001;     % seconds

audStimP.StTime    = trigs(1);                          % Seconds
audStimP.SpTime    = trigs(2);                          % Seconds

audStimP.InflaT    = InflaT;  % seconds
audStimP.InflaV    = InflaV;  % cents

audStimP.pertSchedBegin = [0, 1.0];
audStimP.pertSchedPre   = [audStimP.StTime, 1.0];

%Define the slope for the Aud. perturbation stimulus
if pertSw == 0 %Linear Standard Stimulus
    audStimP.pertRampT1  = 0.11; %s
    audStimP.pertRampT2  = 0.15; %s
    
    if trialType == 0
        audStimP.pertMagCent = 0;
    else
        audStimP.pertMagCent = -100;
    end
    
    audStimP.pertMagMult = 2^(audStimP.pertMagCent/1200);
    audStimP.pertEndDw   = audStimP.StTime + audStimP.pertRampT1;
    audStimP.pertEndUp   = audStimP.SpTime + audStimP.pertRampT2;
    
    audStimP.pertSchedStRamp = [audStimP.pertEndDw, audStimP.pertMagMult];
    audStimP.pertSchedSpRamp = [audStimP.pertEndUp, 1.0];
elseif pertSw == 1 %Sigmoid (Laryngeal) Matched Stimulus
    audStimP.pertRampT   = round(InflaT, 3); % seconds

    if trialType == 0
        audStimP.pertMagCent = 0;
        
        audStimP.pertEndDw   = audStimP.StTime + audStimP.pertRampT;
        audStimP.pertEndUp   = audStimP.SpTime + audStimP.pertRampT;

        audStimP.pertSchedStRamp = [audStimP.pertEndDw, audStimP.pertMagMult];    
        audStimP.pertSchedSpRamp = [audStimP.pertEndUp, 1.0]; 
    else
        audStimP.pertMagCent = round(InflaV, 2);
        
        spread  = 15;
        spreadH = spread/2;
        t = audStimP.tStep:audStimP.tStep:audStimP.pertRampT;
        numSteps = length(t);
        x = linspace(0, spread, numSteps);
        sigmoid = audStimP.pertMagCent*sigmf(x, [1 spreadH]);
        
        sigmoidCent = 2.^(sigmoid/1200);
        audStimP.pertEndDw   = audStimP.StTime + t;
        audStimP.pertEndUp   = audStimP.SpTime + t;

        audStimP.pertSchedStRamp = [audStimP.pertEndDw', sigmoidCent'];    
        audStimP.pertSchedSpRamp = [audStimP.pertEndUp', fliplr(sigmoidCent)']; 
    end  
    
    audStimP.pertMagMult = 2^(audStimP.pertMagCent/1200);
end

audStimP.pertSchedSteady = [audStimP.SpTime, audStimP.pertMagMult];

audStimP.pertSched   = [audStimP.pertSchedBegin;...
                        audStimP.pertSchedPre;...
                        audStimP.pertSchedStRamp;...
                        audStimP.pertSchedSteady;...
                        audStimP.pertSchedSpRamp];
                    
audStimP.pertStim = audStimP.pertSched;
audStimP.pertStim = cat(1, audStimP.pertStim, [4 1.0]);
audStimP.pertStim(:,2) = 1200*log2(audStimP.pertStim(:,2));
end

function drawStimulus(audStimP, dirs)
% close all

time     = audStimP.pertStim(:,1);
stim     = audStimP.pertStim(:,2);

pertName = audStimP.pertName;
pertMag  = audStimP.pertMagCent;
pertTime = audStimP.pertRampT;

pertAx  = [audStimP.StTime, audStimP.SpTime];
pertAy  = [200 200];

plotpos = [10 50];
plotdim = [1200 500];
pertColor = [0.8 0.8 0.8];
AudStim = figure('Color', [1 1 1]);
set(AudStim, 'Position',[plotpos plotdim],'PaperPositionMode','auto')

pA = area(pertAx, pertAy, -200, 'FaceColor', pertColor, 'EdgeColor', pertColor);
hold on

plot(time, stim, 'LineWidth', 3)

xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('f0 Shift (cents)', 'FontSize', 16, 'FontWeight', 'bold')
title({'Auditory Feedback Perturbation Stimulus'; pertName}, 'FontSize', 18, 'FontWeight', 'bold')
axis([0 4 -120 1]); box off;

annotation('textbox',[0.62 0.20 0.40 0.1],...
           'String', {['Fall/Rise Time: ' num2str(pertTime) ' seconds'],...
                      ['Perturbation Magnitude: ' num2str(pertMag) ' cents']},...
                    'LineStyle','none',...
                    'FontWeight','bold',...
                    'FontSize',14,...
                    'FontName','Arial');

set(gca, 'FontSize', 14,...
         'FontWeight', 'bold');
     
     
plTitle = ['AFPerturbStim_' pertName '.jpg'];     
saveFileName = fullfile(dirs.SavResultsDir, plTitle);
export_fig(saveFileName) 
end