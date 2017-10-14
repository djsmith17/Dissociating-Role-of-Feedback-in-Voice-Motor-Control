function dfGeneratef0AcuityTokens()
close all;

prompt = {'Subject ID:',...
          'Baseline Run:',...
          'Baseline Trial:',...
          'Gender ("male" or "female")'};
name = 'Subject Information';
numlines = 1;
defaultanswer = {'null', 'BV1', '3', 'female'};
answer = inputdlg(prompt, name, numlines, defaultanswer);

if isempty(answer)
    return
end

GT.project = 'Dissociating-Role-of-Feedback-in-Voice-Motor-Control';
GT.subject = answer{1};
GT.baseRec = answer{2};
GT.baseTrial = str2double(answer{3});
GT.gender  = answer{4};

dirs = dfDirs(GT.project);
% Folder paths to save data files
dirs.RecFileDir = fullfile(dirs.RecData, GT.subject);
dirs.SavFileDir = fullfile(dirs.RecData, GT.subject, GT.baseRec, [GT.subject GT.baseRec 'DRF.mat']);

dirs.tokenDir = fullfile(dirs.RecFileDir, 'speechTokens');
dirs.baseTokenFile = fullfile(dirs.tokenDir,[GT.subject 'BaseToken.wav']);

if ~exist(dirs.RecFileDir, 'dir')
    mkdir(dirs.RecFileDir);
end

if ~exist(dirs.SavFileDir, 'file')
    disp('ERROR: No voice file at this location!')
    return
end

if exist(dirs.tokenDir, 'dir')
    rmdir(dirs.tokenDir, 's')  
end
mkdir(dirs.tokenDir);

GT.xMax = 200; %max difference between speaker's fo and fo of stimulus in headphones
GT.xMin = 1; %min difference between speaker's fo and fo of stimulus in headphones
GT.xAll = -100:0.5:100;
GT.xAll = GT.xAll(~logical(GT.xAll == 0));
GT.xLen = length(GT.xAll);

% Generate audio tokens
[BaseToken, fs] = dfGenerateBT(dirs, GT.baseTrial); %Extract a Speech Token. Located in JND Folder
subjf0 = dfcalcf0Praat(dirs);                       %Calculate f0 using praat. Located in JND Folder
PertFreqs = targetf0calc(subjf0, GT.xAll, GT.xLen); %Located Below
numPertFreqs = length(PertFreqs);
PertTokens = dfGeneratePT(dirs, numPertFreqs, PertFreqs, GT); %Generate Pert Tokens. Located in JND Folder

GT.subjf0     = subjf0;
GT.pertFreqs  = PertFreqs;
GT.fs         = fs;
GT.BaseToken  = BaseToken;
GT.PertTokens = PertTokens;
end