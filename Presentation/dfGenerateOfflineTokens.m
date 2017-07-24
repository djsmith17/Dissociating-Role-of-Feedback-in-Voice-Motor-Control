function modifiedToken = dfGenerateOfflineTokens(baseToken, dirs, gender, level)
%This scripts loads a previously recorded audio signal and provides a
%Pitch-shift to it in a similiar fashion that happens during online
%testing.

%Paradigm Configurations
expParam.baseToken          = baseToken;
expParam.gender             = gender;
expParam.level              = level; %In semitones
expParam.trialLen           = 2; %seconds

expParam.sRate              = 48000;  % Hardware sampling rate (before downsampling)
expParam.downFact           = 3;
expParam.sRateAnal          = expParam.sRate/expParam.downFact; %Everything get automatically downsampled! So annoying
expParam.frameLen           = 96;  % Before downsampling
expParam.audioInterfaceName = 'MOTU MicroBook'; %'ASIO4ALL' 'Komplete'

%Resample at 48000Hz
expParam.baseToken_reSamp = resample(expParam.baseToken, 48000, 16000);
%Split the signal into frames
expParam.baseToken_frames = makecell(expParam.baseToken_reSamp, expParam.frameLen);

%Set up Audapter
Audapter('deviceName', expParam.audioInterfaceName);
Audapter('setParam', 'downFact', expParam.downFact, 0);
Audapter('setParam', 'sRate', expParam.sRateAnal, 0);
Audapter('setParam', 'frameLen', expParam.frameLen / expParam.downFact, 0);
p = getAudapterDefaultParams(expParam.gender);
p.bPitchShift = 1;

%Set up OST and PCF Files
expParam.ostFN = fullfile(dirs.Prelim, 'AFPerturbOST.ost'); check_file(expParam.ostFN);
expParam.pcfFN = fullfile(dirs.Prelim, 'AFPerturbPCF.pcf'); check_file(expParam.pcfFN);

setAudapterFilesForTokens(expParam.ostFN, expParam.pcfFN, expParam.trialLen, expParam.level)

%Set the OST and PCF functions
Audapter('ost', expParam.ostFN, 0);
Audapter('pcf', expParam.pcfFN, 0);

% fprintf('Running Trial %d\n', ii)
AudapterIO('init', p);
Audapter('reset');

for n = 1:length(expParam.baseToken_frames)
    Audapter('runFrame', expParam.baseToken_frames{n});
end

try
    data = AudapterIO('getData');
    modifiedToken = data.signalOut;
catch
    sprintf('\nAudapter decided not to show up today')
    data = [];
    modifiedToken = [];
    return
end
end

function setAudapterFilesForTokens(ost, pcf, trialLen, level)

OST_tline = writeOSTportions(trialLen);
PCF_tline = writePCFportions(level);

svPSRLevels(ost, OST_tline);
svPSRLevels(pcf, PCF_tline);
end

function OST_tline = writeOSTportions(trialLen)
%The Online Status Tracking file regulates the timing of when actions or
%changes to the speech occur. The steps between timed actions followed
%rules outlined in the Audapter Manuel. This has been specifically
%organized for a simple full recording perturbation.

OST_tline{1} = '# Online status tracking (OST) configuration file';
OST_tline{2} = 'rmsSlopeWin = 0.030000';
OST_tline{3} = ' ';
OST_tline{4} = '# Main section: Heuristic rules for tracking';
OST_tline{5} = 'n = 3';
OST_tline{6} = '0 INTENSITY_RISE_HOLD 0.01 0.05 {} # Detect voicing onset';
OST_tline{7} = ['2 ELAPSED_TIME ' num2str(trialLen) ' NaN {} #The amount of time pre-perturbation']; %Random start between 1.7 and 2.1s
OST_tline{8} = '3 OST_END NaN NaN {} #End the dang thing';
end

function PCF_tline = writePCFportions(level)
%The Pertrubation Configuration file defines the levels for acoustic 
%variables at each action step defined in the OST. This have been
%specifically organized for customized Pitch-Shift Reflex experiments. 
%The PCF expects f0 in units of semitones. Semitones is expected as input.

%Simple pert through whole recording 
PCF_tline{1} = '# Section 1 (Time warping): tBegin, rate1, dur1, durHold, rate2';
PCF_tline{2} = '0';
PCF_tline{3} = ' ';
PCF_tline{4} = '# Section 2: stat pitchShift(st) gainShift(dB) fmtPertAmp fmtPertPhi(rad)';
PCF_tline{5} = '4';
PCF_tline{6} = ['0, ' num2str(level) ', 0.0, 0, 0'];
PCF_tline{7} = ['1, ' num2str(level) ', 0.0, 0, 0'];
PCF_tline{8} = ['2, ' num2str(level) ', 0.0, 0, 0'];
PCF_tline{9} = ['3, ' num2str(level) ', 0.0, 0, 0'];
end

function tline = svPSRLevels (file, tline)
fid = fopen(file, 'w');
for i = 1:length(tline)    
    if i ~= length(tline)
        fprintf(fid,'%s\n',tline{i});
    else
        fprintf(fid,'%s',tline{i});
    end
end 
fclose(fid);
end