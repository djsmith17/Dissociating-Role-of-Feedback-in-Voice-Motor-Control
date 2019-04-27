function dfRunSubjPlotting()
%The function that draws the plots from analyzed results

%This uses the following plot functions
%drawDAQAll
%drawDAQPresMic
%drawDAQAlignedPressure
%drawDAQMeanTrialMicf0
%drawAudRespIndivTrial
%drawAudRespMeanTrial

close all;
sPlt.project       = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
sPlt.participants  = {'DRF1'}; %List of multiple participants.
sPlt.numPart       = length(sPlt.participants);
sPlt.runs          = {'AF1', 'AF2'}; %All runs to consider 
sPlt.numRuns       = length(sPlt.runs);
dirs               = dfDirs(sPlt.project);
ext                = '';

%Plot Toggles. This could eventually become an input variable
sv2File                      = 1;
sPlt.drawDAQAll              = 0; % All signals recorded by the NIDAQ
sPlt.drawDAQPresMic          = 0; % Pressure vs Microphone Data

sPlt.drawDAQAlignedPressure  = 0; % Superimposed Pressure recordings from perturbed trials
sPlt.drawPertMicResponse     = 0; % All Perturbed Trials Microphone input
sPlt.drawMicHeadResponse     = 1; % Mean Perturbed Trials. Microphone vs Headphones

presFlag = 1;
 
for ii = 1:sPlt.numPart
    participant = sPlt.participants{ii};
    for jj = 1:sPlt.numRuns 
        run = sPlt.runs{jj};
        dirs.SavResultsFile = fullfile(dirs.Results, participant, run, [participant run ext 'ResultsDRF.mat']); %The Analyzed Results File
        
        dirs.PlotResultsDir = fullfile(dirs.Results, participant, run); % Analyzed Results Folder...Where Plots will go
        

        if exist(dirs.SavResultsFile, 'file') == 0
            fprintf('\nERROR: File %s does not exist!\n', dirs.SavResultsFile)
            return
        else
            load(dirs.SavResultsFile)
        end        
        
        if sPlt.drawDAQAll == 1
            drawDAQAll(res, dirs.PlotResultsDir, sv2File)
        end
        
        if sPlt.drawDAQPresMic == 1
            drawDAQPresMic(res, dirs.PlotResultsDir)
        end
        
        if sPlt.drawDAQAlignedPressure == 1
            drawDAQAlignedPressure(res, dirs.PlotResultsDir, sv2File)
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
close all
end