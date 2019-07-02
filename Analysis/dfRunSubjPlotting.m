function dfRunSubjPlotting()
%The function that draws the plots from analyzed results

% This calls the following drawing functions:
% -drawDAQAlignedPressure
% -drawMeanTrialMicf0
% -drawAllPertTrialMicf0
% -drawAudRespIndivTrial
% -drawAudRespMeanTrial

close all;
sPlt.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
sPlt.participants  = {'DRF1',...
                      'DRF2',...
                      'DRF4',...
                      'DRF5',...
                      'DRF6',...
                      'DRF7',...
                      'DRF8',...
                      'DRF9',...
                      'DRF10',...
                      'DRF12',...
                      'DRF13',...
                      'DRF14',...
                      'DRF15',...
                      'DRF16',...
                      'DRF17',...
                      'DRF18',...
                      'DRF19',...
                      'DRF20'}; %List of multiple participants.
sPlt.numPart       = length(sPlt.participants);
sPlt.runs          = {'AF1','AF2'}; %All runs to consider 
sPlt.numRuns       = length(sPlt.runs);
dirs               = dfDirs(sPlt.project);
ext                = '';

%Plot Toggles. This could eventually become an input variable
sv2File                      = 1;

sPlt.drawDAQAlignedPressure  = 0; % Superimposed Pressure recordings from perturbed trials
sPlt.drawPertMicResponse     = 0; % All Perturbed Trials Microphone input
sPlt.drawMicHeadResponse     = 1; % Mean Perturbed Trials. Microphone vs Headphones

presFlag = 1;
 
for ii = 1:sPlt.numPart
    participant = sPlt.participants{ii};
    for jj = 1:sPlt.numRuns 
        run = sPlt.runs{jj};
        dirs.PlotResultsDir = fullfile(dirs.Results, participant, run); % Analyzed Results Folder...Where Plots will go
        
        dirs.SavResultsFile = fullfile(dirs.Results, participant, run, [participant run ext 'ResultsDRF.mat']); % Analyzed Results File.. The analyzed data to open
        if exist(dirs.SavResultsFile, 'file') == 0
            error('File %s does not exist!\n', dirs.SavResultsFile)
        else
            load(dirs.SavResultsFile)
        end
        
        if sPlt.drawDAQAlignedPressure == 1
            drawDAQAlignedPressure(res, dirs.PlotResultsDir, sv2File, 0)
        end
        
        if sPlt.drawPertMicResponse == 1
            drawMeanTrialMicf0(res, dirs.PlotResultsDir, presFlag)
            drawAllPertTrialMicf0(res, dirs.PlotResultsDir, presFlag)
        end            
        
        if sPlt.drawMicHeadResponse == 1
            drawAudRespIndivTrial(res, dirs.PlotResultsDir)
            drawAudRespMeanTrial(res, dirs.PlotResultsDir)
        end
        
        close all
    end
end
end