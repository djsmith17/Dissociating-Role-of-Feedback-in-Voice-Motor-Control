function dfGeneratef0AcuityTokens()
%dfGeneratef0AcuityTokens() generates both a baseline token and
%pitch-shifted tokens for JND type pitch perception tasks. This looks for a
%baseline voice recording made earlier, and saves both wav files and a
%single MATLAB data structure containing the sound tokens needed for a JND
%task.
%
%This script makes use of the following functions
%dfDirs.m
%dfGenerateBF.m
%dfcalcf0Praat.m
%dfGeneratePF.m

close all;
prompt = {'Subject ID:',...
          'Session ID:',...
          'Baseline Run:',...
          'Baseline Trial:'};
name = 'Subject Information';
numlines = 1;
defaultanswer = {'null', 'GT1', 'BV1', '3'};
answer = inputdlg(prompt, name, numlines, defaultanswer);

if isempty(answer)
    return
end

GT.project = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
GT.subject = answer{1};
GT.run     = answer{2};
GT.baseRec = answer{3};
GT.baseTrial = str2double(answer{4});

dirs = dfDirs(GT.project);
% Folder paths to save data files
dirs.RecFileDir = fullfile(dirs.RecData, GT.subject, GT.run);
dirs.SavFileDir = fullfile(dirs.RecData, GT.subject, GT.baseRec, [GT.subject GT.baseRec 'DRF.mat']);

dirs.tokenDir = fullfile(dirs.RecFileDir, 'speechTokens');
dirs.baseTokenFile = fullfile(dirs.tokenDir,[GT.subject GT.run 'BaseToken.wav']);

if ~exist(dirs.RecFileDir, 'dir')
    mkdir(dirs.RecFileDir);
end

if ~exist(dirs.SavFileDir, 'file')
    disp('ERROR: No voice file at this location!')
    return
end

if ~exist(dirs.tokenDir, 'dir')
    mkdir(dirs.tokenDir);
end

GT.xMax = 200; %max difference between speaker's fo and fo of stimulus in headphones
GT.xMin = 1; %min difference between speaker's fo and fo of stimulus in headphones
GT.xAll = -100:0.5:100;
GT.xAll = GT.xAll(~logical(GT.xAll == 0));
GT.xLen = length(GT.xAll);

% Generate audio tokens
[BaseToken, fs] = dfGenerateBT(dirs, GT.baseTrial); %Extract a Speech Token. Located in JND Folder
subjf0 = dfcalcf0Praat(dirs);                       %Calculate f0 using praat. Located in JND Folder
PertFreqs  = targetf0calc(subjf0, GT.xAll, GT.xLen); %Located Below
PertTokens = dfGeneratePT(dirs, GT, PertFreqs); %Generate Pert Tokens. Located in JND Folder

GT.subjf0     = subjf0;
GT.pertFreqs  = PertFreqs;
GT.fs         = fs;
GT.BaseToken  = BaseToken;
GT.PertTokens = PertTokens;

GTFiles = fullfile(dirs.RecFileDir, [GT.subject GT.run 'DRF.mat']);
save(GTFiles, 'GT');
end

function freqs = targetf0calc(f0, AllFreq, FreqLen)
%calculates all possible freq vals spaced 0.5 cent apart. 

for i = 1: FreqLen
    if i ~= 0 %I dont want the case of pure baseline
        freqs(i) = f0*2^(AllFreq(i)/1200);
    end
end
end