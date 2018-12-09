function dfCreatef0AcuityTokens()
% dfCreatef0AcuityTokens() generates both a baseline token and
% pitch-shifted tokens for JND type pitch perception tasks. This loads a
% baseline voice recording made earlier, and saves both .wav files and a
% single MATLAB data structure containing the sound tokens needed for a JND
% task.
%
% This script is dependent on the following external functions:
% -dfDirs.m
% -dfGenerateBT.m
% -dfcalcf0PraatSingle.m
% -dfGeneratePT.m
%
% This script has the following subfunctions:
% -calcShiftedf0

close all;
prompt = {'Subject ID:',...
          'Session ID:',...
          'Baseline Run:',...
          'Baseline Trial:'};
name = 'Subject Information';
numlines = 1;
defaultanswer = {'null', 'GT1', 'BV1', '3'};
answer = inputdlg(prompt, name, numlines, defaultanswer);

ET = tic;
if isempty(answer)
    return
end

GT.project   = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
GT.subject   = answer{1};             % Participant Identifier
GT.gender    = [];
GT.age       = [];
GT.run       = answer{2};             % What is the name of this set of tokens?
GT.baseRec   = answer{3};             % Which baseline recording
GT.baseTrial = str2double(answer{4}); % Which trial to use from baseline

dirs = dfDirs(GT.project);
% Folder paths to save data files
dirs.RecFileDir    = fullfile(dirs.RecData, GT.subject, GT.run);  % Dir to save any files at all
dirs.TokenDir      = fullfile(dirs.RecFileDir, 'speechTokens');   % Dir to save the tokens
dirs.baseTokenFile = fullfile(dirs.TokenDir, [GT.subject GT.run 'BaseToken.wav']); % File name of the Baseline toke

if ~exist(dirs.RecFileDir, 'dir')
    mkdir(dirs.RecFileDir);
end

if ~exist(dirs.TokenDir, 'dir')
    mkdir(dirs.TokenDir);
end

% Where to find the baseline recordings
dirs.BaseFile = fullfile(dirs.RecData, GT.subject, GT.baseRec, [GT.subject GT.baseRec 'DRF.mat']);

[GT.BVDataMic,...
 GT.BVfs,...
 GT.gender,...
 GT.age]      = loadBaselineRecording(dirs, GT.baseTrial);

GT.xMax = 200;  %max difference between speaker's f0 and f0 of stimulus in headphones
GT.xMin = 1;    %min difference between speaker's f0 and f0 of stimulus in headphones

GT.fs            = 44100;             % sampling rate to play tokens at
GT.tokenLen      = 0.5;               % seconds
GT.tokenLenP     = GT.tokenLen*GT.fs; % length of speech tokens in points 
GT.allCentShifts = -100:0.5:100;      % Change this to change range of shifted cents
GT.allCentShifts = GT.allCentShifts(GT.allCentShifts ~= 0); % Excluding a shift of 0 cents
GT.numPertToken  = length(GT.allCentShifts);

% Generate audio tokens
GT.BaseToken  = dfGenerateBT(dirs, GT);    % Extract a Speech Token.
GT.subjf0     = dfCalcf0PraatSingle(dirs); % Calculate f0 using praat.
GT.PertFreqs  = calcShiftedf0(GT.subjf0, GT.allCentShifts); % Located Below
GT.PertTokens = dfGeneratePT(dirs, GT);    % Generate Pert Tokens using praat. 

GTFiles = fullfile(dirs.RecFileDir, [GT.subject GT.run 'DRF.mat']);
save(GTFiles, 'GT');
fprintf('Completed creating tokens for participant %s with f0 of %0.2f Hz\n', GT.subject, GT.subjf0)

elapsed_time = toc(ET);
fprintf('Elapsed Time: %f (s)\n', elapsed_time)
end

function [BVMic, BVfs, gender, age] = loadBaselineRecording(dirs, baseTrial)

baseFile = dirs.BaseFile;

if ~exist(baseFile, 'file')
    error('ERROR: No baseline voice file at the designated location!')
else
    fprintf('Loading Baseline recording from...\n%s\n\n', dirs.BaseFile)
    load(baseFile) % Returns a structure called 'DRF'
end

rawData  = DRF.rawData;
BVData   = rawData(baseTrial); % Grab a specific trial
BVMic    = BVData.signalIn;    % Grab the microphone channel.

expParam = DRF.expParam;
BVfs     = expParam.sRateAnal; % Sampling Rate
gender   = expParam.gender;
age      = expParam.age;
end

function shiftedFreqs = calcShiftedf0(f0, allCentShifts)
% calcShiftedf0(f0, allCentShifts) calculates the values of shifted 
% fundamental frequency, used to generate pitch-shifted audio tokens.

numToken     = length(allCentShifts);
shiftedFreqs = zeros(numToken, 1);
for ii = 1:numToken
    shiftedFreqs(ii) = f0*2^(allCentShifts(ii)/1200);
end
end