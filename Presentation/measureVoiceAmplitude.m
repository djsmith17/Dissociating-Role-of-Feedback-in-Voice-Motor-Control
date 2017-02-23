function [allrmsMean, finalrmsMean] = measureVoiceAmplitude()
%This takes a few sample recordings of the participant's voice and returns
%the average RMS value of the recorded voice. This should be used at the
%beginning of recording session to determine a baseline voice amplitude. 
%It is assumed that the participant speaks at a comfortable and typical
%speaking volume, the microphone is placed at a consistent distance for
%each participant, and microphone gain levels are constant for each
%participant and through the trial. 

%Paradigm Configurations
expParam.project    = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
expParam.expType    = 'Somatosensory Perturbation_Perceptual';
expParam.numTrial   = 3;
expParam.gender     = 'male';
expParam.masking    = 0;
expParam.trialLen   = 4; %Seconds

dirs = sfDirs(expParam.project, expParam.expType);

expParam.sRate              = 48000;  % Hardware sampling rate (before downsampling)
expParam.downFact           = 3;
expParam.sRateAnal          = expParam.sRate/expParam.downFact;
expParam.frameLen           = 96;  % Before downsampling
expParam.audioInterfaceName = 'MOTU MicroBook'; %'ASIO4ALL' 'Komplete'

%Set up Audapter
Audapter('deviceName', expParam.audioInterfaceName);
Audapter('setParam', 'downFact', expParam.downFact, 0);
Audapter('setParam', 'sRate', expParam.sRateAnal, 0);
Audapter('setParam', 'frameLen', expParam.frameLen / expParam.downFact, 0);
p = getAudapterDefaultParams(expParam.gender);

%Set up OST and PCF Files. Just for the sake of having them
expParam.ostFN = fullfile(dirs.Prelim, 'SFPerturbOST.ost'); check_file(expParam.ostFN);
expParam.pcfFN = fullfile(dirs.Prelim, 'SFPerturbPCF.pcf'); check_file(expParam.pcfFN);

[expParam, p]      = setAudFeedType(expParam, dirs, p); %Sets some p params

refSPL  = 0.00002; %20 micropascals

allrmsMean = [];
for ii = 1:expParam.numTrial
    
    fprintf('Ready to Record?\n')
    pause()
    %Set the OST and PCF functions
    Audapter('ost', expParam.ostFN, 0);
    Audapter('pcf', expParam.pcfFN, 0);
    
    fprintf('Trial %d\n', ii)
    AudapterIO('init', p);
    Audapter('reset');
    Audapter('start');

    pause(expParam.trialLen);
    
    Audapter('stop');   
    
    data = AudapterIO('getData');
    
    rms   = data.rms(:,1);
    rmsdB = 20*log10(rms/refSPL);
    rmsMean = mean(rmsdB);
    
    allrmsMean = [allrmsMean, rmsMean]; 
end

finalrmsMean = mean(allrmsMean);

fprintf('\nThe mean amplitude from each of the three voice recordings were %4.2f dB, %4.2f dB, and %4.2f dB\n', allrmsMean)
fprintf('\nThe mean amplitude from all three voice recordings is %4.2f dB\n', finalrmsMean)
end