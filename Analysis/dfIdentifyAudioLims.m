function lims = dfIdentifyAudioLims(An)
% lims = dfIdentifyAudioLims(An) calculates limits of analyzed data 
% so that the limits used in plotting are dynamic and fit the data. 
%
% lims is a structure with two values, lims.audio and lims.audioMean
% Each of those limit sets are of the shape [X1 X2 Y1 Y2]
% This sub function is redundant within other functions and might
% eventually become its own function

timeTrials = An.timef0;
timeMeaned = An.secTime;
singleSt = round(timeTrials(1), 1);
singleSp = round(timeTrials(end), 1);

meanSt = round(timeMeaned(1), 1);
meanSp = round(timeMeaned(end), 1);

timeAnalysisPeriod = find(timeTrials > (singleSt + 0.5) & timeTrials < (singleSp - 0.5));

%%%%%%%%%%%lims.audio%%%%%%%%%%%
%Individual Full Trials (Perturbed): f0 Audio
if ~isempty(An.audioMf0p)
    pertTrialsM = An.audioMf0p;
    pertTrialsH = An.audioHf0p;

    uLMa = max(pertTrialsM(timeAnalysisPeriod,:));
    lLMa = min(pertTrialsM(timeAnalysisPeriod,:));

    uLM = round(max(uLMa)) + 20;
    lLM = round(min(lLMa)) - 20;
       
    uLHa = max(pertTrialsH(timeAnalysisPeriod,:));
    lLHa = min(pertTrialsH(timeAnalysisPeriod,:));
    
    uLH = round(max(uLHa)) + 20;
    lLH = round(min(lLHa)) - 20;
    
    % Consider the combination of microphone and headphone recordings
    lLMH = CompareAndChooseBounds(lLM, lLH, 'min');
    uLMH = CompareAndChooseBounds(uLM, uLH, 'max');
    
    if strcmp(An.expType, 'Auditory Perturbation_Perceptual')
        lims.audio = [singleSt singleSp lLMH uLMH];
    else
        lims.audio = [singleSt singleSp lLM uLM];
    end
else
    lims.audio     = [singleSt singleSp -20 20];
end

%%%%%%%%%%%lims.audioMean%%%%%%%%%%%
%Mean Sectioned Trials (Perturbed): f0 Audio 
if ~isempty(An.audioMf0_meanp)
    % Check the microphone recordings
    [lwBoundOn, upBoundOn] = MaxMinBounds(An.audioMf0_meanp(:,1), An.audioMf0_meanp(:,2), 10);
    [lwBoundOf, upBoundOf] = MaxMinBounds(An.audioMf0_meanp(:,3), An.audioMf0_meanp(:,4), 10);

    lwBoundM = CompareAndChooseBounds(lwBoundOn, lwBoundOf, 'min');
    upBoundM = CompareAndChooseBounds(upBoundOn, upBoundOf, 'max');

    % Check the headphone recordings
    [lwBoundOn, upBoundOn] = MaxMinBounds(An.audioHf0_meanp(:,1), An.audioHf0_meanp(:,2), 10);
    [lwBoundOf, upBoundOf] = MaxMinBounds(An.audioHf0_meanp(:,3), An.audioHf0_meanp(:,4), 10);

    lwBoundH = CompareAndChooseBounds(lwBoundOn, lwBoundOf, 'min');
    upBoundH = CompareAndChooseBounds(upBoundOn, upBoundOf, 'max');
    
    % Consider the combination of microphone and headphone recordings
    lwBoundMH = CompareAndChooseBounds(lwBoundM, lwBoundH, 'min');
    upBoundMH = CompareAndChooseBounds(upBoundM, upBoundH, 'max');

    if strcmp(An.expType, 'Auditory Perturbation_Perceptual')
        lims.audioMean = [meanSt meanSp lwBoundMH upBoundMH];
    else
        lims.audioMean = [meanSt meanSp lwBoundM upBoundM];
    end
else
    lims.audioMean     = [meanSt meanSp -50 50];
end
end

function [lwBound, upBound] = MaxMinBounds(audio, audioErr, buff)
% MaxMinBounds takes an audio trace and the error trace for that audio, and
% identifies what the max and min values of the trace are. It identifies
% how wide the bounds needs to be fully show the trace with some buff

[~, Imax] = max(audio);
upBound   = round(audio(Imax) + audioErr(Imax) + buff);
[~, Imin] = min(audio);
lwBound   = round(audio(Imin) - audioErr(Imin) - buff);
end

function bound = CompareAndChooseBounds(bound1, bound2, type)
% CompareAndChooseBounds does a simple and quick check between two proposed
% bounds and decides which is greater/lesser. Since I am often comparing
% two sets of things, this helps to keep data ranges consistent and
% symmetrical. 

allBounds = [bound1 bound2];

if strcmp(type, 'max')
    bound = max(allBounds);
elseif strcmp(type, 'min')
    bound = min(allBounds);
else
    bound = allBounds(1);
end

end