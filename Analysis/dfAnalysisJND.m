function dfAnalysisJND()
% dfAnalysisJND() loads a f0 JND data set (fAX#DRF.mat) and analyzed the raw
% data into interpretable results. After the results have been computed and
% organized, they are immediately plotted. Multiple JND results can be
% analyzed and plotted at the same time. By following the saved file
% structure of fAX1, fAX2,... etc, you can select the run numbers to
% analyze and view the run results next to each other in the plots.
%
% This function calls the following external functions
% -dfDirs
% -dfAnalyzeThresholdJND
% -drawJNDResults
%
% This script has the following subfunctions:
% -initJNDAnalysis()
% -AnalyzeRawJNDData()
tic
close all

% Initalize the analysis structure
JNDa.project      = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
JNDa.participants = {'DRF1',...
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
                     'DRF20'}; % List of multiple participants.
JNDa.numPart      = length(JNDa.participants);
JNDa.runs         = {'fAX1','fAX2','fAX3','fAX4'}; %List of multiple runs.
JNDa.numRuns      = length(JNDa.runs);

dirs = dfDirs(JNDa.project);

for jj = 1:JNDa.numPart
    curPart = JNDa.participants{jj}; % Current Participant
    dirs.baselineData  = fullfile(dirs.SavData, curPart, 'BV1', [curPart 'BV1' 'DRF.mat']); % Where to find data  
    dirs.GTRecordings  = fullfile(dirs.SavData, curPart, 'GT1', [curPart 'GT1' 'DRF.mat']); % Where to find data  
    dirs.SavResultsDir = fullfile(dirs.Results, curPart, 'JND'); %Where to save results

    % Look for the baseline voice info for this participant, then load it
    if exist(dirs.baselineData, 'file') == 0
        fprintf('ERROR: Could not find baseline data set at %s\n', dirs.baselineData)
        return
    else
        fprintf('Loading baseline data set for %s %s\n', curPart, 'BV1')
        load(dirs.baselineData) % Returns DRF
        bV = DRF;
        load(dirs.GTRecordings)
    end
    
    bVaudioRMS = assessTokenQuality(GT, bV);
    
    if exist(dirs.SavResultsDir, 'dir') == 0
        mkdir(dirs.SavResultsDir) %If the folder we are saving does not exist, let's make it
    end

    for ii = 1:JNDa.numRuns
        curRun     = JNDa.runs{ii}; % Current Run
        dirs.SavFileDir = fullfile(dirs.SavData, curPart, curRun, [curPart curRun 'DRF.mat']); %Where to find raw data

        fprintf('****%s  %s****\n', curPart, curRun)
        fprintf('Loading Raw JND Data\n')
        load(dirs.SavFileDir) % Returns UD
        
        % Organize the raw data into a result structure named 'resJND'
        resJND = AnalyzeRawJNDData(UD);
        resJND.baseAudioRMS = bVaudioRMS; % dB
        dirs.SavResultsFile = fullfile(dirs.SavResultsDir, [curPart curRun 'ResultsDRF.mat']); %What to name the organized results
        fprintf('Saving Individual JND Results\n\n')
        save(dirs.SavResultsFile, 'resJND')    
    end
end
fprintf('Elapsed time was %f min\n', toc/60)
end

function bVLoudDB = assessTokenQuality(GT, bV)

GTTime = linspace(0, GT.tokenLen, GT.tokenLenP);
GTBase = GT.BaseToken;

frameLen = bV.expParam.frameLenDown;
rmsB     = bV.expParam.rmsB;
bVfs   = bV.expParam.sRateAnal;
bVBase = bV.rawData(GT.baseTrial).signalIn;
bVBase_seg = bVBase(2*bVfs:2.5*bVfs);
bvTime = linspace(2, 2.5, length(bVBase_seg));

bVaudioRMS = bV.qRes.audioRMS(GT.baseTrial);

bVStruct.rms   = reCalcRMS(bVBase_seg(800:7200), frameLen);
bVStruct.rmsM  = mean(bVStruct.rms);
GTStruct.rms   = reCalcRMS(GTBase(2205:19845), frameLen);
GTStruct.rmsM  = mean(GTStruct.rms);

bVLoudDB = dfCalcMeanRMS(bVStruct, rmsB);
GTLoudDB = dfCalcMeanRMS(GTStruct, rmsB);

fprintf('Baseline Voice Loudness recorded as %0.2f dB\n', bVaudioRMS)
fprintf('Baseline Voice Loudness reanalyzed as %0.2f dB\n', bVLoudDB)
fprintf('Baseline Token Loudness reanalyzed as %0.2f dB\n', GTLoudDB)

% fprintf('Difference between recorded and reanalyzed baseline voice = %f\n', bVLoudDB - bVaudioRMS)
% fprintf('Difference between baseline Token and baseline voice = %f\n', bVLoudDB - GTLoudDB)

% figure('Color', [1 1 1])
% subplot(1,2,1)
% plot(bvTime, bVBase_seg)
% xlabel('Time (s)')
% box off
% title({'Baseline Voice', ['RMS = ' num2str(round(bVStruct.rmsM, 4))]});
% 
% subplot(1,2,2)
% plot(GTTime, GTBase)
% xlabel('Time (s)')
% box off
% title({'Baseline Token', ['RMS = ' num2str(round(GTStruct.rmsM, 4))]});
% suptitle(bV.expParam.subject)

allTokenRMSdB = [];
for ii = 1:GT.numPertToken
    tokenStruct.rms = reCalcRMS(GT.PertTokens(ii,2205:19845), frameLen);
    allTokenRMSdB = cat(1, allTokenRMSdB, dfCalcMeanRMS(tokenStruct, rmsB));
end
allTokenRMSdB_Mean = round(mean(allTokenRMSdB), 2);
allTokenRMSdB_SD   = round(std(allTokenRMSdB), 3);

fprintf('Mean Token Loudness = %0.2f db (SD = %0.3f dB)\n\n', allTokenRMSdB_Mean, allTokenRMSdB_SD)
end

function sigRMS = reCalcRMS(sig, frameLen)

lenSig = length(sig);
numFrame = floor(lenSig/frameLen);

sigRMS = zeros(numFrame, 1);
for ii = 1:numFrame
    sigRMS(ii) = rms(sig([1:frameLen]*ii));
end

end

function resJND = initJNDAnalysis()

resJND.participant = [];
resJND.gender      = [];
resJND.f0          = [];
resJND.baseAudioRMS = [];
resJND.run         = [];

resJND.instructions     = {};
resJND.selectOpt        = {};

resJND.reversalsReached = [];
resJND.trialsCompleted  = [];
resJND.timeElapsed      = [];

resJND.distProgression = [];

resJND.trialsAtReversals = [];
resJND.distAtReversals   = [];

resJND.trialsAtCorrectOpt1   = [];
resJND.distAtCorrectOpt1     = [];
resJND.trialsAtIncorrectOpt1 = [];
resJND.distAtIncorrectOpt1   = [];
resJND.trialsAtCorrectOpt2   = [];
resJND.distAtCorrectOpt2     = [];
resJND.trialsAtIncorrectOpt2 = [];
resJND.distAtIncorrectOpt2   = [];

resJND.JNDScore        = [];
resJND.LastSetAccuracy = [];
resJND.catchAccuracy   = [];
end

function resJND = AnalyzeRawJNDData(UD)
% This calls the following subfunctions
% -initJNDAnalysis()
% -dfAnalyseThresholdJND()

fprintf('Starting Individual JND Analysis\n')

resJND = initJNDAnalysis();

resJND.participant = UD.subject;
resJND.gender      = UD.gender;
resJND.f0          = round(UD.subjf0, 1);
resJND.run         = UD.run;
    
resJND.instructions = UD.inst;    
resJND.selectOpt = {'First'; 'Last'};

resJND.reversalsReached = UD.reversals;
resJND.trialsCompleted  = UD.performedTrials;
resJND.timeElapsed      = UD.elapsedTime;

resJND.distProgression = UD.x;

resJND.trialsAtReversals = find(UD.reversal~=0);
resJND.distAtReversals   = UD.x(resJND.trialsAtReversals);

resJND.trialsAtCorrectOpt1   = find(UD.allTrialTypes == 1);
resJND.trialsAtIncorrectOpt1 = find(UD.allTrialTypes == 2);
resJND.trialsAtCorrectOpt2   = find(UD.allTrialTypes == 3);
resJND.trialsAtIncorrectOpt2 = find(UD.allTrialTypes == 4);

resJND.distAtCorrectOpt1   = UD.x(resJND.trialsAtCorrectOpt1)';
resJND.distAtIncorrectOpt1 = UD.x(resJND.trialsAtIncorrectOpt1)';
resJND.distAtCorrectOpt2   = UD.x(resJND.trialsAtCorrectOpt2)';
resJND.distAtIncorrectOpt2 = UD.x(resJND.trialsAtIncorrectOpt2)';

% Determine the JND Score and Accuracy of the last set of trials
[JNDScore, LastSetAccuracy] = dfAnalyzeThresholdJND(UD, 'reversals', 4); %Cents
resJND.JNDScore        = JNDScore;
resJND.LastSetAccuracy = round(LastSetAccuracy, 1);
resJND.catchAccuracy   = round(0); %Currently N/A, but maybe used again later
end