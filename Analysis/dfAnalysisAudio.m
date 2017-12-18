function dfAnalysisAudio(dirs, An)
%I found that I was mostly writing the same things twice in my scripts for
%audapter analysis and NIDAQ analysis. Ultimately what I care about is the
%audio in one form or another. Let's just put it here.

%Instatiate the variables we intend to use. 
An = initAudVar(An);

%Main script that does the Signal Frequency Analysis
dirs.audiof0AnalysisFile = fullfile(dirs.SavResultsDir, [An.subject An.run 'f0Analysis.mat']);

if exist(dirs.audiof0AnalysisFile, 'file') == 0
    [f0A.time_audio, f0A.audioMf0, f0A.fsA] = signalFrequencyAnalysis(dirs, An.time, An.audioM, An.sRate, fV, An.bTf0b, 1);
    [f0A.time_audio, f0A.audioHf0, f0A.fsA] = signalFrequencyAnalysis(dirs, An.time, An.audioH, An.sRate, fV, An.bTf0b, 1);
    save(dirs.audiof0AnalysisFile, 'f0A')
else
    load(dirs.audiof0AnalysisFile)
end

An.time_audio = f0A.time_audio;
An.fsA        = f0A.fsA;
An.audioMf0   = f0A.audioMf0; 
An.audioHf0   = f0A.audioHf0;    

%Smooth the f0 data
An.audioMf0S   = smoothf0(An.audioMf0);
An.audioHf0S   = smoothf0(An.audioHf0);

%Normalize f0 and convert to cents
prePert       = (0.5 < An.time_audio & 1.0 > An.time_audio);
An.trialf0b   = mean(An.audioMf0S(prePert,:),1);
An.f0b        = mean(An.trialf0b);

An.audioMf0_norm = normalizeDAQf0(An.audioMf0S, An.trialf0b);
An.audioHf0_norm = normalizeDAQf0(An.audioHf0S, An.trialf0b);

%Find the Perturbed Trials
An.audioMf0_p = parseTrialTypes(An.audioMf0_norm, An.pertIdx);
An.audioHf0_p = parseTrialTypes(An.audioHf0_norm, An.pertIdx);
An.audioMf0_c = parseTrialTypes(An.audioMf0_norm, An.contIdx);
An.audioHf0_c = parseTrialTypes(An.audioHf0_norm, An.contIdx);

%Find troublesome trials and remove
[An.audioMf0_pPP, An.audioHf0_pPP, An.numPertTrialsPP, An.pertTrigPP] = audioPostProcessing(An.time_audio, An.audioMf0_p, An.audioHf0_p, An.numPertTrials, An.pertTrig, An.curSess, 'Pert');
[An.audioMf0_cPP, An.audioHf0_cPP, An.numContTrialsPP, An.contTrigPP] = audioPostProcessing(An.time_audio, An.audioMf0_c, An.audioHf0_c, An.numContTrials, An.contTrig, An.curSess, 'Cont');

%Section the data around onset and offset
[An.secTime, An.audioMf0_Secp] = sectionAudioData(An.time_audio, An.audioMf0_pPP, An.fsA, An.pertTrigPP);
[An.secTime, An.audioHf0_Secp] = sectionAudioData(An.time_audio, An.audioHf0_pPP, An.fsA, An.pertTrigPP);
[An.secTime, An.audioMf0_Secc] = sectionAudioData(An.time_audio, An.audioMf0_cPP, An.fsA, An.contTrigPP);
[An.secTime, An.audioHf0_Secc] = sectionAudioData(An.time_audio, An.audioHf0_cPP, An.fsA, An.contTrigPP);

%Mean around the onset and offset
An.audioMf0_meanp = meanAudioData(An.audioMf0_Secp);
An.audioHf0_meanp = meanAudioData(An.audioHf0_Secp);
An.audioMf0_meanc = meanAudioData(An.audioMf0_Secc);
An.audioHf0_meanc = meanAudioData(An.audioHf0_Secc); 

%The Inflation Response
[An.respVar, An.respVarMean, An.respVarSD, An.InflaStimVar] = InflationResponse(An.secTime, An.audioMf0_Secp);

end

function An = initAudVar(An)
%Initialize some variables to keep track of them

An.time_audio     = []; %time vector of audio samples recorded
An.fsA            = []; %sampling rate of audio samples
An.audioMf0       = []; %Raw Microphone Audio Data
An.audioHf0       = []; %Raw Headphone Audio Data
An.trialf0b       = []; %Per Trial calculated f0
An.f0b            = []; %Average trial f0
An.audioMf0_norm  = []; %Normalized mic data
An.audioHf0_norm  = []; %Normalized head data
An.audioMf0_p     = []; %Perturbed mic data (norm)
An.audioHf0_p     = []; %Pertrubed head data (norm)
An.audioMf0_c     = []; %Control mic data (norm)
An.audioHf0_c     = []; %Control head data (norm)
An.audioMf0_pPP   = []; %
An.audioHf0_pPP   = [];
An.numPertTrialsPP = []; 
An.pertTrigPP      = [];
An.audioMf0_cPP    = []; 
An.audioHf0_cPP    = [];
An.numContTrialsPP = [];
An.contTrigPP      = [];

An.secTime        = [];
An.audioMf0_Secp  = []; 
An.audioHf0_Secp  = [];
An.audioMf0_Secc  = []; 
An.audioHf0_Secc  = [];

An.audioMf0_meanp = []; 
An.audioHf0_meanp = [];
An.audioMf0_meanc = []; 
An.audioHf0_meanc = [];

An.respVar        = []; 
An.respVarMean    = [];
An.respVarSD      = [];
An.InflaStimVar   = [];

end

function [timef0, audiof0, fsA] = signalFrequencyAnalysis(dirs, time, audio, fs, fV, bTf0b, flag)
[~, numTrial] = size(audio);

if flag == 1
    [timef0, audiof0, fsA] = dfCalcf0Praat(dirs, audio, fs, bTf0b);
else
    %Low-Pass filter for the given cut off frequency
    [B,A]    = butter(4,(fV.freqCutOff)/(fs/2));

    audiof0 = zeros(fV.numWin, numTrial);
    for j = 1:numTrial %Trial by Trial
        sensorHP = filtfilt(B,A,audio(:,j));
        
        timef0 = [];
        for i = 1:fV.numWin
            winIdx = fV.winSts(i):fV.winSts(i)+ fV.winP - 1;
            timef0 = cat(1, timef0, mean(time(winIdx)));
            f0 = dfCalcf0Chile(sensorHP(winIdx), fs);
            if isnan(f0)
%                 disp('hit')
                f0 = 100;
            end
            audiof0(i,j) = f0;
        end
    end
    fsA = fV.fsA;
end
end