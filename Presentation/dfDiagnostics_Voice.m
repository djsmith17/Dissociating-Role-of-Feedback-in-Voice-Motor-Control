function [allrmsMean, finalrmsMean] = dfDiagnostics_Voice()
%This takes a few sample recordings of the participant's voice and returns
%the average RMS value of the recorded voice. This should be used at the
%beginning of recording session to determine a baseline voice amplitude. 
%This script assumes:
%1: The participant speaks at a comfortable and typical speaking volume
%2: The microphone is placed at a consistent distance for each participant
%3: The microphone gain levels are constant for each participant and through the trials
%4: The participant phonates a steady-state vowel sound through these recordings

%Paradigm Configurations
expParam.project    = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
expParam.expType    = 'Somatosensory Perturbation_Perceptual';
expParam.numTrial   = 3;
expParam.gender     = 'male';
expParam.masking    = 0;
expParam.trialLen   = 4; %Seconds

dirs = dfDirs(expParam.project);

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

[expParam, p]      = dfSetAudFB(expParam, dirs, p); %Sets some p params

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
    
    data    = AudapterIO('getData');
    rmsMean = calcMeanRMS(data, refSPL);
    
    allrmsMean = cat(1, allrmsMean, rmsMean); 
end

finalrmsMean = mean(allrmsMean);

fprintf('\nThe mean amplitude from each of the three voice recordings were %4.2f dB, %4.2f dB, and %4.2f dB\n', allrmsMean)
fprintf('\nThe mean amplitude from all three voice recordings is %4.2f dB\n', finalrmsMean)
end

function rmsMean = calcMeanRMS(data, refSPL)
rms   = data.rms(:,1);
rmsdB = 20*log10(rms/refSPL);
rmsMean = mean(rmsdB);
end

function quikFFT(data)
x = data.signalIn;
fs = data.params.sRate;
Y = fft(x);
L = length(x);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;

figure
plot(f,P1)
end